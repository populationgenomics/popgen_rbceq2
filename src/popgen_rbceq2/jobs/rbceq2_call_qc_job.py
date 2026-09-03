#!/usr/bin/env python3
"""Flag blood-group systems whose call rests on a low-quality allele-defining site.

Joins the DP/GQ extracted at every allele-defining coordinate by
`FilterAndConvertGvcfsForRbceq2` to the committed site -> blood-group-system map, and
writes a `<sg>.qc.tsv` in the same wide layout as RBCeq2's genotype and phenotype TSVs
(one row per sample, one column per system).

Cell values are `PASS`, `NA` for a system with no assessable defining site in the map, or
the semicolon-joined flags of the system's defining sites:

    NOCOV:1:3774964(A>G)                                   no record covers the site
    LOWQ:1:3774964(A>G,DP=19,GQ=15)                        a call at the site itself
    DEL:1:3774964(A>G,del=CATGA>C,GT=0/1,DP=30,GQ=50)      a deletion removed the base
    LOWQ:1:3774964(A>G,block=101bp,DP=26,MIN_DP=19,GQ=15)  a reference block covers it

A flag name states two independent findings about the site, joined by `+`: its **severity**
(NOCOV, DEL or LOWQ, or absent when the site is fine) and its **provenance** (POSTHOC when
the post-hoc caller supplied the record, absent when the primary one did).

    POSTHOC:1:159204893(T>C,src=gatk-hc-4.6.2.0,DP=42,GQ=99)
    LOWQ+POSTHOC:1:3774964(A>G,src=gatk-hc-4.6.2.0,block=91bp,DP=1,MIN_DP=1,GQ=3)

The halves are joined rather than ranked, because ranking makes one displace the other and a
recovered site that is also poor has to report both. So grepping `POSTHOC` finds every call
that rests on recovered data, and equivalently every system that was not typable from the
primary caller alone; `rests_on_posthoc` answers the same question on a whole cell. A site
that is neither poor nor recovered is not listed at all, which is what makes a bare `PASS`
cell mean "nothing to report".

Within the parentheses the first field is always the allele the db defines the antigen on,
and every later field is `key=value`. Which keys appear says what the numbers describe,
because a flag names both the site and what the caller actually reported there, and the two
can differ. The field count varies deliberately: NOCOV carries no metrics because no record
was found, which is not the same as a record that omitted a field (rendered `.`).
"""

import csv
import logging
from collections import defaultdict
from dataclasses import dataclass
from typing import Literal

import click
from cpg_utils import to_path

from popgen_rbceq2.logging_setup import setup_logging
from popgen_rbceq2.scripts.bg_db import SITE_SYSTEM_COLUMNS, DefiningSite, chrom_key

logger = logging.getLogger(__name__)

# Flag prefixes appearing in the QC TSV cells.
PASS = 'PASS'  # noqa: S105 — a QC flag, not a credential
LOWQ = 'LOWQ'
NOCOV = 'NOCOV'
# The defining coordinate falls inside the REF span of a deletion the sample carries, so
# the base the antigen is defined on is not present. RBCeq2 sees no variant at the
# coordinate and calls the system reference; that call rests on a base the sample deleted,
# which is neither a quality problem nor an absence of data.
DELETED = 'DEL'
# A system RBCeq2 called but that has no assessable defining site in the map: either the db
# gives it no coordinate, or its only defining alleles are structural variants that one
# base of DP/GQ cannot speak to (ABCC1, ATP11C and CD99 at rbceq2 2.4.3). Never assessed,
# so neither PASS nor LOWQ is honest.
NOT_ASSESSED = 'NA'
# The primary caller had no record at this site and the value came from the post-hoc caller
# instead (see PosthocGenotypeOffTargetSites and the merge in FilterAndConvertGvcfsForRbceq2).
#
# This is provenance, and like severity it is a property of the site, so it lives in the site's
# own flag name. A recovered site must never read as though the primary caller supported the
# antigen: the supplement is a different caller, without the sample's DRAGstr model, on reads
# the capture design did not target. It says the same thing from the other side too — that the
# system was not typable without the recall — which is worth reporting even where the recovered
# site is perfectly good.
#
# Provenance is joined to severity rather than ranked against it. The two are independent
# findings, and ranking them makes one displace the other: a recovered site that is also
# sub-threshold would report only LOWQ, which on the validation cohorts hid the recovery for
# 234 of 681 reliant systems.
POSTHOC = 'POSTHOC'

# Joins the severity and provenance halves of a flag name, as in `LOWQ+POSTHOC`. Neither `:`
# nor `;` nor `,` can serve: those already separate the name from the coordinate, one site
# flag from the next, and one metric from the next.
FLAG_JOIN = '+'

# Alleles longer than this are rendered as a length in the flag; the FY and XK large
# indels carry ~200-base REFs that would otherwise dominate the cell.
MAX_ALLELE_CHARS = 12

# Reference-block spans at or above this are rendered in kb rather than bases.
KB = 1000

# The symbolic allele on GVCF reference blocks and on the twin record that
# `bcftools norm -m -any` leaves beside every split variant.
NON_REF = '<NON_REF>'

EXTRACT_COLUMNS = ('chrom', 'pos', 'ref', 'alt', 'end', 'gt', 'dp', 'gq', 'min_dp', 'posthoc')

# INFO/POSTHOC as bcftools renders it for a record that does not carry the tag, which is every
# record on a genome run and every primary record on an exome one.
NO_POSTHOC_TAG = '.'

# Where a site's DP and GQ came from, which decides how the flag reads:
#   site      a call at the coordinate, so the values describe the site itself
#   block     a reference block covering it, whether it starts there or reaches it from an
#             earlier POS, so the values describe the whole band rather than this one base
#   deletion  a deletion whose REF span swallows it, so the defining base is not present
#
# Note that `site` is not simply "a record starting here". A reference block can begin on any
# coordinate, including a defining one, and that is still a band. resolve_coverage reads this
# off the record via `is_block`, not off the branch that found the record.
CoverageSource = Literal['site', 'block', 'deletion']


@dataclass(frozen=True, slots=True)
class GvcfRecord:
    """One GVCF record from the extract, either a variant or a reference block.

    DP, GQ, MIN_DP and END are None when the record does not carry that field.

    `posthoc` names the caller that supplied the record, or None for one from the primary
    gVCF. It is set from INFO/POSTHOC, which the merge in `FilterAndConvertGvcfsForRbceq2`
    stamps on every record it takes from the post-hoc caller.
    """

    chrom: str
    pos: int
    ref: str
    alt: str
    end: int | None
    gt: str
    dp: int | None
    gq: int | None
    min_dp: int | None
    posthoc: str | None = None

    @property
    def is_primary(self) -> bool:
        """True if this record came from the primary caller rather than the post-hoc one."""
        return self.posthoc is None

    @property
    def is_block(self) -> bool:
        """True if this record is a reference block, meaning it spans past its own POS.

        A one-base block (END == POS) is False, which is what `resolve_coverage` needs: a
        record covering only its own POS is a call at that site, not a span reaching it.
        """
        return self.end is not None and self.end > self.pos

    def covers(self, pos: int) -> bool:
        """Test whether this record spans a position.

        A record covers the position it starts at, any position inside its INFO/END span,
        and any position inside its REF allele.

        Args:
            pos: 1-based position to test.

        Returns:
            True if the record spans `pos`.
        """
        return self.pos <= pos <= max(self.end or self.pos, self.pos + len(self.ref) - 1)

    @property
    def span(self) -> int:
        """The number of bases this record covers, counting its own POS."""
        return max(self.end or self.pos, self.pos + len(self.ref) - 1) - self.pos + 1


@dataclass(frozen=True, slots=True)
class Coverage:
    """The DP and GQ attributed to one defining site, and the record they came from.

    Either value may be None if the record omitted the field. `record` is kept so a flag
    can say what the caller actually reported at the coordinate, which is not always the
    allele the site is defined on, and so a block's MIN_DP can be reported beside its DP.
    """

    dp: int | None
    gq: int | None
    source: CoverageSource
    record: GvcfRecord


def _optional_int(value: str) -> int | None:
    """Parse a bcftools field, mapping its missing marker to None.

    Args:
        value: Raw field text.

    Returns:
        The integer value, or None for `.` and the empty string.
    """
    return None if value in ('.', '') else int(value)


def parse_extract(text: str) -> list[GvcfRecord]:
    """Parse the headerless `bcftools query` DP/GQ extract.

    Args:
        text: Extract contents, one tab-separated record per line in EXTRACT_COLUMNS order.

    Returns:
        The parsed records, in file order. Blank lines are skipped.

    Raises:
        ValueError: A line does not have exactly len(EXTRACT_COLUMNS) fields.
    """
    records = []
    for line in text.splitlines():
        if not line.strip():
            continue
        fields = line.split('\t')
        if len(fields) != len(EXTRACT_COLUMNS):
            raise ValueError(
                f'Expected {len(EXTRACT_COLUMNS)} columns {EXTRACT_COLUMNS} in the DP/GQ extract, '
                f'got {len(fields)}: {line!r}. An extract written before the columns last '
                'changed will not parse; delete it so the conversion stage rewrites it.'
            )
        chrom, pos, ref, alt, end, gt, dp, gq, min_dp, posthoc = fields
        records.append(
            GvcfRecord(
                chrom=chrom,
                pos=int(pos),
                ref=ref,
                alt=alt,
                end=_optional_int(end),
                gt=gt,
                dp=_optional_int(dp),
                gq=_optional_int(gq),
                min_dp=_optional_int(min_dp),
                posthoc=None if posthoc == NO_POSTHOC_TAG else posthoc,
            )
        )
    return records


def load_site_systems(text: str) -> list[DefiningSite]:
    """Read the committed site -> blood-group-system map.

    Args:
        text: Contents of `bg_site_systems.<genome>.tsv`, including its header row.

    Returns:
        One DefiningSite per data row, in file order.

    Raises:
        ValueError: The header does not match SITE_SYSTEM_COLUMNS.
    """
    reader = csv.DictReader(text.splitlines(), delimiter='\t')
    if tuple(reader.fieldnames or ()) != SITE_SYSTEM_COLUMNS:
        raise ValueError(f'Site-system map must have columns {SITE_SYSTEM_COLUMNS}, got {reader.fieldnames}')
    sites = [
        DefiningSite(
            chrom=row['chrom'],
            pos=int(row['pos']),
            ref=row['ref'],
            alt=row['alt'],
            kind=row['kind'],  # type: ignore[arg-type]
            system=row['system'],
        )
        for row in reader
    ]
    # gen_bg_resources.py never writes an SV site, because one base of DP/GQ cannot assess a
    # multi-kb allele. Reject one rather than reporting a quality for it.
    if svs := [s for s in sites if s.kind == 'sv']:
        raise ValueError(
            f'Site-system map carries {len(svs)} structural-variant site(s), which cannot be assessed on '
            f'one base of DP/GQ, e.g. {svs[0].chrom}:{svs[0].pos} ({svs[0].system}). Regenerate the map '
            'with gen_bg_resources.py rather than editing it.'
        )
    return sites


def _coverage_from(record: GvcfRecord, source: CoverageSource) -> Coverage:
    """Read the DP and GQ off a record covering a defining site.

    Args:
        record: The record chosen as covering a defining site.
        source: What that record is, relative to the site.

    Returns:
        The site's coverage.
    """
    # A block reports one DP, one MIN_DP and one GQ for every base it covers, and nothing
    # per-base. Expanding it to per-base records would repeat those three numbers, not
    # recover anything.
    #
    # Illumina documents a hom-ref block's DP as the *median* depth across the band and
    # MIN_DP as its minimum. (DRAGEN's own v3.7 page calls FORMAT/DP the band minimum
    # instead, which our data refutes: DP exceeded MIN_DP in 69% of 608,745 blocks, so the
    # two cannot both be minima.)
    #
    # A defining site's position inside the block is unknowable, so DP, the median, is what
    # it is judged on, being the best available estimate of the site's own depth. MIN_DP is a
    # floor that one shallow base anywhere in the span can set; over those same blocks it fell
    # below min_depth while DP did not for 0.2% of them, which is enough to flag a system off
    # a single read. It is reported beside DP in the flag rather than tested, so a borderline
    # site stays visible to whoever reads the column.
    #
    # GQ needs no equivalent care, because blocks are banded on GQ and not on depth. Our
    # gVCFs band at 10/20/30/40 and declare it per file:
    #
    #   ##GVCFBlock=minGQ=20(inclusive),maxGQ=30(exclusive)
    #
    # so a block's single GQ answers a threshold on a band edge exactly: one reporting
    # GQ >= 20 cannot contain a base below 20. See min_gq in the config.
    return Coverage(dp=record.dp, gq=record.gq, source=source, record=record)


def resolve_coverage(records: list[GvcfRecord], chrom: str, pos: int) -> Coverage | None:
    """Find the DP and GQ the caller reported at one defining site.

    Prefers a record starting at the site with a real ALT, then any record starting at the
    site, then the innermost record spanning it.

    Args:
        records: Candidate records; those on other contigs are ignored.
        chrom: Contig of the defining site, in GVCF naming.
        pos: 1-based position of the defining site.

    Returns:
        The site's coverage tagged with where it came from, or None if no record covers it.
    """
    covering = [r for r in records if r.chrom == chrom and r.covers(pos)]
    if not covering:
        return None
    # The primary caller wins wherever both speak, matching the merge's own fill rule. The
    # merge keeps a post-hoc reference block whole, so one straddling a capture edge can reach
    # a site the primary caller did call; without this the innermost-record rule below would
    # prefer that block purely because it starts later, and report a recovered site where
    # there was never a hole.
    covering = [r for r in covering if r.is_primary] or covering
    at_site = [r for r in covering if r.pos == pos]
    if at_site:
        # `bcftools norm -m -any` splits every GVCF variant into the real ALT plus a
        # <NON_REF> twin at the same POS, carrying identical DP/GQ. Prefer the real ALT,
        # but do not count on it reaching the extract: at a site that is a deletion's own
        # anchor base, `--targets-overlap 2` drops the real-ALT record — an indel's variant
        # span starts at POS+1 — and only the twin arrives, so the fallback is the normal
        # path there. Either serves for DP/GQ; the gt/alt of an at-site record are NOT
        # reliable evidence of what rbceq2 read.
        real_alt = [r for r in at_site if r.alt != NON_REF]
        chosen = (real_alt or at_site)[0]
        # Starting at the site does not make a record a call at it. Reference blocks are
        # banded on GQ, so one can begin on any coordinate, including a defining one; that
        # is still a span whose DP is a median over the band and whose MIN_DP the flag has
        # to report. Read the source off the record, not off the branch that found it.
        return _coverage_from(chosen, 'block' if chosen.is_block else 'site')
    # Nothing starts at the site, so the chosen record reaches it from an earlier POS: a
    # reference block spanning it, or a deletion whose REF swallows it. A deletion starting
    # at the site itself is not this case. Its first REF base is the retained anchor, and
    # it is handled above.
    spanning = max(covering, key=lambda r: r.pos)
    return _coverage_from(spanning, 'block' if spanning.is_block else 'deletion')


def _render_allele(allele: str) -> str:
    """Render an allele, abbreviating one longer than MAX_ALLELE_CHARS to its length.

    Args:
        allele: REF or ALT sequence.

    Returns:
        The allele, or `<n>bp` if it is too long to show.
    """
    return f'{len(allele)}bp' if len(allele) > MAX_ALLELE_CHARS else allele


def describe_site(site: DefiningSite) -> tuple[str, str]:
    """Render the coordinate and allele a flag identifies a site by.

    This describes the site as the db defines it, not what the caller observed there; a
    flag adds the latter from `Coverage.record`.

    Args:
        site: The defining site.

    Returns:
        The `<chrom>:<pos>` coordinate, with the `chr` prefix dropped to match RBCeq2's
        internal loci naming, and the allele as `<ref>><alt>`, or `ref` for a lane site,
        which has no ref/alt of its own.
    """
    coordinate = f'{site.chrom.removeprefix("chr")}:{site.pos}'
    allele = 'ref' if site.kind == 'ref' else f'{_render_allele(site.ref)}>{_render_allele(site.alt)}'
    return coordinate, allele


def _render_span(record: GvcfRecord) -> str:
    """Render a reference block's span, in kb once it reaches KB bases.

    Args:
        record: The reference block.

    Returns:
        The span as `<n>bp` or `<n.n>kb`.
    """
    span = record.span
    return f'{span / KB:.1f}kb' if span >= KB else f'{span}bp'


def _render_int(value: int | None) -> str:
    """Render an optional metric for a flag.

    Args:
        value: Metric value, or None if the record omitted the field.

    Returns:
        The value, or `.` for None so absence reads apart from a reported zero.
    """
    return '.' if value is None else str(value)


def _render_metrics(coverage: Coverage) -> str:
    """Render the depth and GQ a site was judged on, and where they came from.

    A reference block also reports its MIN_DP, which is not what the site was judged on but
    bounds how low any single base in the span went, so the gap between the two is visible
    at the point of reading rather than having to be assumed.

    Args:
        coverage: The site's coverage.

    Returns:
        `DP=<n>,GQ=<n>`, prefixed for a record that is not a call at the site: `del=`/`GT=`
        for a deletion, or `block=`/`MIN_DP=` for a reference block. Every field is
        `key=value`, so which keys appear says what the numbers describe. A missing metric
        renders as `.`, to distinguish it from a reported zero.

        A record from the post-hoc caller leads with `src=<caller>-<version>`. The flag name
        says *that* the site was recovered; `src=` says by which caller and version, so a
        version bump is visible in the output rather than only in the code.
    """
    record = coverage.record
    depth = f'DP={_render_int(coverage.dp)}'
    source = '' if record.is_primary else f'src={record.posthoc},'
    match coverage.source:
        case 'deletion':
            observed = f'del={_render_allele(record.ref)}>{_render_allele(record.alt)},GT={record.gt},'
        case 'block':
            observed = f'block={_render_span(record)},'
            depth += f',MIN_DP={_render_int(record.min_dp)}'
        case _:
            observed = ''
    return f'{source}{observed}{depth},GQ={_render_int(coverage.gq)}'


def flag_name(severity: str | None, *, recovered: bool) -> str | None:
    """Join a site's severity and provenance into one flag name.

    Args:
        severity: NOCOV, DEL or LOWQ, or None when the site clears both thresholds.
        recovered: Whether the post-hoc caller supplied the record.

    Returns:
        The two joined by `+` when both apply, whichever applies alone, or None when neither
        does — a primary site at adequate quality, which is not listed in the cell.
    """
    parts = [part for part in (severity, POSTHOC if recovered else None) if part is not None]
    return FLAG_JOIN.join(parts) if parts else None


def flag_severity(flag: str) -> str | None:
    """The severity half of one site flag, or None if the site only carries provenance.

    Args:
        flag: One site flag from a QC cell, e.g. `LOWQ+POSTHOC:1:3774964(...)`.

    Returns:
        NOCOV, DEL or LOWQ, or None for a bare POSTHOC flag.
    """
    names = [name for name in flag.split(':', 1)[0].split(FLAG_JOIN) if name != POSTHOC]
    return names[0] if names else None


def rests_on_posthoc(cell: str) -> bool:
    """Whether a system's call rests on a site the post-hoc caller supplied.

    True is the annotation this pipeline exists to make visible: the antigen rests on a
    different caller, without the sample's DRAGstr model, over reads the capture design did
    not target. It says the same thing from the other side — that the system was not typable
    from the primary caller alone, because without the recall that site would have been NOCOV.

    Args:
        cell: One cell of the QC TSV.

    Returns:
        True if any of the system's site flags carries POSTHOC.
    """
    if cell in (PASS, NOT_ASSESSED):
        return False
    return any(POSTHOC in flag.split(':', 1)[0].split(FLAG_JOIN) for flag in cell.split(';'))


def flag_site(site: DefiningSite, coverage: Coverage | None, min_depth: int, min_gq: int) -> str | None:
    """Flag one defining site on its coverage and on which caller supplied it.

    Args:
        site: The defining site.
        coverage: Its coverage, or None if no GVCF record covers it.
        min_depth: DP below which a site is flagged.
        min_gq: GQ below which a site is flagged.

    Returns:
        `<name>:<coordinate>(<allele>,<metrics>)`, or None for a site with nothing to report,
        meaning one the primary caller covered at adequate quality.

        The name joins two independent findings about the site with `+`:

        - **severity** — NOCOV, DEL or LOWQ, in that order of precedence, or absent when the
          site clears both thresholds;
        - **provenance** — POSTHOC when the record came from the post-hoc caller, absent when
          it came from the primary one.

        So `LOWQ` is a primary site below threshold, `POSTHOC` a recovered site that passes,
        and `LOWQ+POSTHOC` a recovered site that is also below threshold. NOCOV never carries
        POSTHOC: there is no record from either caller, so there is no provenance to report.

        Both findings are properties of the site, so both belong in the site's own flag. They
        are joined rather than ranked because ranking makes one displace the other, and a
        recovered site that is also poor has to report both — grepping POSTHOC would otherwise
        miss it, which on the validation cohorts was 234 of 681 reliant systems.
    """
    coordinate, allele = describe_site(site)
    if coverage is None:
        return f'{NOCOV}:{coordinate}({allele})'
    # DEL outranks LOWQ, and is reported however good the deletion's own DP and GQ are: a
    # confidently-called deletion of the defining base is the finding, not a quality one.
    if coverage.source == 'deletion':
        severity = DELETED
    elif (dp := coverage.dp) is None or (gq := coverage.gq) is None or dp < min_depth or gq < min_gq:
        # A missing DP or GQ flags like a low one: the caller reported no quality for a site
        # that defines a blood-group antigen, which is not evidence of a good call.
        severity = LOWQ
    else:
        severity = None
    name = flag_name(severity, recovered=not coverage.record.is_primary)
    if name is None:
        return None
    return f'{name}:{coordinate}({allele},{_render_metrics(coverage)})'


def load_fillable_sites(text: str) -> frozenset[tuple[str, int]]:
    """Read the defining sites the post-hoc caller was allowed to fill.

    Args:
        text: Contents of the off-design defining-sites BED `SelectOffDesignDefiningSites`
            writes: one 0-based half-open single-base interval per site, no header.

    Returns:
        The sites as (contig, 1-based position), the coordinates `DefiningSite` carries.

    Raises:
        ValueError: A row is not a single-base interval, so it names no defining site.
    """
    sites = set()
    for line_no, line in enumerate(text.splitlines(), start=1):
        if not line.strip():
            continue
        chrom, start, end = line.split('\t')[:3]
        if int(end) != int(start) + 1:
            raise ValueError(f'Fillable-sites BED line {line_no} is not a single-base interval: {line!r}')
        sites.add((chrom, int(end)))
    return frozenset(sites)


def flags_by_system(
    sites: list[DefiningSite],
    records: list[GvcfRecord],
    min_depth: int,
    min_gq: int,
    fillable_sites: frozenset[tuple[str, int]] | None = None,
) -> tuple[dict[str, str], list[DefiningSite]]:
    """Aggregate per-site flags into one cell per blood-group system.

    Args:
        sites: Every defining site to assess.
        records: Extracted GVCF records to resolve coverage from.
        min_depth: DP below which a site is flagged.
        min_gq: GQ below which a site is flagged.
        fillable_sites: The (contig, position) set the merge was allowed to fill from the
            post-hoc caller, or None on a run that merged nothing. A post-hoc record counts as
            covering a site only if the site is in this set. The merge keeps a post-hoc
            reference block whole, so one selected for an off-design hole can reach an
            in-design hole beside it, where it is the only record; without this the QC would
            report that hole as recovered when the rule says it stays NOCOV.

    Returns:
        A `{system: cell}` map covering every system in `sites`, where a system with no
        flagged site is PASS and several flagged sites are semicolon-joined in coordinate
        order; and the sites with no covering record, for the caller to log.

        Every finding is carried by the site flags themselves, so this only groups and orders
        them. `rests_on_posthoc` answers whether a cell's system needed the recall.

    Raises:
        ValueError: `fillable_sites` is None but a record carries INFO/POSTHOC. The merge
            that produced it read a fillable set, so the QC has to be handed the same one; a
            run that has none cannot say which post-hoc records to believe.
    """
    if fillable_sites is None and any(not record.is_primary for record in records):
        raise ValueError(
            'The extract carries post-hoc records but no fillable-sites BED was given. The QC '
            'needs the off-design defining sites the merge filled from, to disregard a post-hoc '
            'record at any other site.'
        )

    by_chrom: dict[str, list[GvcfRecord]] = defaultdict(list)
    for record in records:
        by_chrom[record.chrom].append(record)

    per_system: dict[str, list[tuple[tuple[int, str], int, str]]] = defaultdict(list)
    uncovered = []
    for site in sites:
        candidates = by_chrom[site.chrom]
        if fillable_sites is not None and (site.chrom, site.pos) not in fillable_sites:
            candidates = [record for record in candidates if record.is_primary]
        coverage = resolve_coverage(candidates, site.chrom, site.pos)
        if coverage is None:
            uncovered.append(site)
        flag = flag_site(site, coverage, min_depth, min_gq)
        if flag is not None:
            per_system[site.system].append((chrom_key(site.chrom), site.pos, flag))
    systems = {site.system for site in sites}
    return (
        {
            system: ';'.join(flag for *_, flag in sorted(per_system[system])) if system in per_system else PASS
            for system in sorted(systems)
        },
        uncovered,
    )


def is_posthoc_only(cell: str) -> bool:
    """Whether a system's cell reports a recovery and no quality problem.

    Such a system was recovered rather than found wanting: every site it lists was supplied by
    the post-hoc caller at adequate quality. Counted apart from the quality flags so a run
    summary does not read as though recovering sites had made the cohort worse.

    A cell carrying both a recovery and a LOWQ, DEL or NOCOV site is False here. It still
    `rests_on_posthoc`, and its flags still say so; it is simply not a clean recovery, so it
    belongs in the quality count.

    Args:
        cell: One cell of the QC TSV.

    Returns:
        True if the cell lists at least one site and every one of them is a bare POSTHOC flag.
    """
    if cell in (PASS, NOT_ASSESSED):
        return False
    return all(flag_severity(flag) is None for flag in cell.split(';'))


def build_qc_tsv(geno_tsv: str, system_flags: dict[str, str]) -> str:
    """Render the QC TSV against the columns of RBCeq2's genotype TSV.

    The first header cell (RBCeq2's `UUID: <uuid>`) and the sample id are copied verbatim,
    so the QC TSV joins and concatenates like the genotype and phenotype TSVs.

    Args:
        geno_tsv: Contents of `<sg>.geno.tsv`, holding a header row and a single sample row.
        system_flags: Flag per blood-group system, as returned by `flags_by_system`.

    Returns:
        The QC TSV, with a column per system in `geno_tsv`'s order. A system absent from
        `system_flags` is NOT_ASSESSED.

    Raises:
        ValueError: `geno_tsv` has no sample row, or its rows are ragged.
    """
    rows = [line.split('\t') for line in geno_tsv.splitlines() if line.strip()]
    if len(rows) < 2:
        raise ValueError('Expected a header row and a single sample row in the rbceq2 geno TSV')
    header, values = rows[0], rows[1]
    if len(header) != len(values):
        raise ValueError(f'rbceq2 geno TSV header has {len(header)} cells but the sample row has {len(values)}')
    systems = header[1:]
    flags = [system_flags.get(system, NOT_ASSESSED) for system in systems]
    return '\n'.join(['\t'.join([header[0], *systems]), '\t'.join([values[0], *flags])]) + '\n'


@click.command()
@click.option('--geno-tsv', required=True, help="RBCeq2's <sg>.geno.tsv; fixes the QC TSV's columns")
@click.option('--site-systems', required=True, help='Committed bg_site_systems.<genome>.tsv')
@click.option('--defining-sites-extract', required=True, help='Per-site DP/GQ extract from the conversion stage')
@click.option('--output', required=True, help='<sg>.qc.tsv to write')
@click.option('--min-depth', type=int, required=True, help='Flag a defining site with DP below this')
@click.option('--min-gq', type=int, required=True, help='Flag a defining site with GQ below this')
@click.option(
    '--fillable-sites',
    default=None,
    help='Off-design defining-sites BED the merge filled from; required when the extract carries post-hoc records',
)
def main(
    geno_tsv: str,
    site_systems: str,
    defining_sites_extract: str,
    output: str,
    min_depth: int,
    min_gq: int,
    fillable_sites: str | None,
) -> None:
    """Write the per-sequencing-group blood-group QC TSV.

    Args:
        geno_tsv: RBCeq2's `<sg>.geno.tsv`, which fixes the QC TSV's columns.
        site_systems: Committed `bg_site_systems.<genome>.tsv`.
        defining_sites_extract: Per-site DP/GQ extract from the conversion stage.
        output: Path to write the QC TSV to.
        min_depth: DP below which a defining site is flagged.
        min_gq: GQ below which a defining site is flagged.
        fillable_sites: The off-design defining-sites BED the merge was allowed to fill from,
            or None on a run that merged nothing. See `flags_by_system`.
    """
    setup_logging(force=True)
    sites = load_site_systems(to_path(site_systems).read_text())
    records = parse_extract(to_path(defining_sites_extract).read_text())
    fillable = None if fillable_sites is None else load_fillable_sites(to_path(fillable_sites).read_text())
    logger.info(f'Assessing {len(sites)} defining sites against {len(records)} extracted GVCF records')
    if fillable is not None:
        logger.info(f'{len(fillable)} defining site(s) were fillable by the post-hoc caller')

    system_flags, uncovered = flags_by_system(sites, records, min_depth, min_gq, fillable)
    if uncovered:
        listed = ', '.join(f'{s.chrom}:{s.pos} ({s.system})' for s in uncovered[:20])
        ellipsis = ' ...' if len(uncovered) > 20 else ''
        logger.warning(
            f'{len(uncovered)} defining site(s) have no covering GVCF record '
            f'and are flagged {NOCOV}: {listed}{ellipsis}'
        )
    qc_tsv = build_qc_tsv(to_path(geno_tsv).read_text(), system_flags)

    # Two populations, counted separately: system_flags covers every system in the site map,
    # while the QC TSV only has columns for the systems rbceq2 put in the geno TSV. A flagged
    # system rbceq2 did not emit is dropped from the TSV, so this log is the only place it
    # shows up.
    emitted = set(qc_tsv.splitlines()[0].split('\t')[1:])
    flagged = sorted(system for system, flag in system_flags.items() if flag != PASS)
    emitted_flagged = [system for system in flagged if system in emitted]
    # Two populations again: a system flagged only because the post-hoc caller stood in for
    # the primary one at a passing site is a recovery, not a quality problem, and lumping the
    # two together would make a successful exome run look like a degraded one.
    recovered = sorted(system for system in flagged if is_posthoc_only(system_flags[system]))
    quality_flagged = [system for system in flagged if system not in set(recovered)]
    named = f': {", ".join(quality_flagged)}' if quality_flagged else ''
    logger.info(
        f'{len(quality_flagged)}/{len(system_flags)} site-map systems quality-flagged, on DP<{min_depth}, '
        f'GQ<{min_gq}, a deleted defining base or no coverage{named}'
    )
    if recovered:
        logger.info(
            f'{len(recovered)} further system(s) rest on a site recovered by the post-hoc caller and '
            f'are flagged {POSTHOC} rather than {PASS}: {", ".join(recovered)}'
        )
    # Counted over every system, not just the cleanly-recovered ones: a system whose recovered
    # site is also sub-threshold appears in the quality count above, and this is the only line
    # that says the recall is what let it be typed at all.
    reliant = sorted(system for system, cell in system_flags.items() if rests_on_posthoc(cell))
    if reliant:
        logger.info(
            f'{len(reliant)}/{len(system_flags)} site-map system(s) rest on a site the post-hoc caller '
            f'supplied, and were not typable from the primary caller alone: {", ".join(reliant)}'
        )
    logger.info(
        f'{len(emitted_flagged)}/{len(emitted)} QC TSV systems carry a flag; '
        f'{len(flagged) - len(emitted_flagged)} flagged system(s) have no geno TSV column and are not in the TSV'
    )

    to_path(output).write_text(qc_tsv)
    logger.info(f'Wrote blood-group QC TSV -> {output}')


if __name__ == '__main__':
    main()
