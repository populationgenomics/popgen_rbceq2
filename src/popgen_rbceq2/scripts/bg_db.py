"""Parse RBCeq2's allele database into the blood-group site resources we ship.

Use this to build the three committed resources (see `gen_bg_resources.py`): the regions
BED that restricts the GVCF conversion, the defining-sites BED that the QC extraction
targets, and the site -> blood-group-system map that turns extracted DP/GQ into a
per-system flag. Deriving all three here keeps them describing one site set.

`parse_positions` and `build_intervals` mirror `rbceq2.IO.vcf` (v2.4.1), so the regions
BED is a superset of the coordinates RBCeq2 itself reads.
"""

import re
from collections import defaultdict
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from typing import Literal

# RBCeq2 keys most karyotype-1 alleles as KLF1 but emits the system as KLF; every place
# it derives a system name applies this (core_logic/constants.py from_string,
# core_logic/alleles.py blood_group, db/db.py prepare_db).
_KLF = 'KLF'

# GRCh38 cells that carry no usable coordinate. RBCeq2 skips the row entirely.
_NO_COORDINATE = frozenset({'', '.', 'na'})

# Leading word of a db token naming a structural variant rather than an allele, in either
# case: `del`/`ins` with a size (`95018451_del_21kb`) and RHD/RHCE's `DEL`/`INS` with a
# base count (`25272546_DEL_148`).
_SV_WORDS = frozenset({'del', 'ins'})

# A REF or ALT sequence. Anything else at a coordinate that is not a lane site or an SV
# notation is a db form this module has not been taught, and raises rather than being
# guessed at.
_SEQUENCE = re.compile(r'^[ACGTN]+$')

SiteKind = Literal['var', 'ref', 'sv']

# Column order of the committed site-system map, shared by the generator that writes it
# and the QC job that reads it back.
SITE_SYSTEM_COLUMNS = ('chrom', 'pos', 'ref', 'alt', 'kind', 'system')


@dataclass(frozen=True, slots=True)
class DefiningSite:
    """One allele-defining coordinate from the db, tagged with its blood-group system.

    `kind` is one of:

    - `var`, a site defined by a variant (`<pos>_<ref>_<alt>` in the db).
    - `ref`, a site defined by the reference base being present at that position
      (`<pos>_ref`, which RBCeq2 calls a lane variant). Carries `.` for both `ref` and
      `alt`, having no alleles of its own.
    - `sv`, a site defined by a large structural variant (`<pos>_del_21kb`), carrying the
      notation's own two words rather than sequences. These are excluded from the shipped
      resources. See `site_system_map`.
    """

    chrom: str
    pos: int
    ref: str
    alt: str
    kind: SiteKind
    system: str


def norm_chrom(chrom: str) -> str:
    """Normalise a db `Chrom` value to GVCF contig form.

    Args:
        chrom: Value of the db's `Chrom` column, e.g. `chr1` or `chrx`.

    Returns:
        The contig name with a `chr` prefix and an upper-case body, e.g. `chrX`.
    """
    return 'chr' + chrom.removeprefix('chr').upper()


def chrom_key(chrom: str) -> tuple[int, str]:
    """Build a sort key ordering autosomes numerically ahead of the sex chromosomes.

    Args:
        chrom: Contig name, with or without a `chr` prefix.

    Returns:
        A key sorting chr1..chr22 by number, then the rest alphabetically.
    """
    bare = chrom.removeprefix('chr')
    return (int(bare), '') if bare.isdigit() else (99, bare)


def blood_group_system(genotype: str, genotype_alt: str) -> str:
    """Derive the blood-group system of an allele, as RBCeq2 names it.

    Mirrors `rbceq2.core_logic.alleles.Allele.blood_group`.

    Args:
        genotype: The db's `Genotype` column, e.g. `VEL*01.01`.
        genotype_alt: The db's `Genotype_alt` column, used when `genotype` is absent.

    Returns:
        The prefix before `*`, with KLF1 normalised to KLF.
    """
    name = genotype.strip()
    if name in ('', '.'):
        name = genotype_alt.strip()
    system = name.split('*')[0]
    return _KLF if _KLF in system.upper() else system


def parse_positions(cell: str | None) -> list[int]:
    """Extract positions from a db genome column, as `rbceq2.IO.vcf.parse_positions` does.

    Args:
        cell: A db `GRCh37`/`GRCh38` cell, comma-separated coordinate tokens.

    Returns:
        The position of every token containing `_` and starting with digits, in cell order.
        Tokens without a numeric position (`unknown_del`, `na`) are skipped.
    """
    if cell is None:
        return []
    positions = []
    for tok in str(cell).split(','):
        if '_' in tok:
            pos = tok.split('_', 1)[0]
            if pos.isdigit():
                positions.append(int(pos))
    return positions


def parse_defining_token(token: str) -> tuple[int, SiteKind, str, str] | None:
    """Parse one db coordinate token.

    Recognises the forms occurring in v2.4.1's db: `3774964_A_G` (a variant),
    `207331122_ref` (a lane site), `95018451_del_21kb` and `25272546_DEL_148` (structural
    variants), any of those with a trailing `_no_phenotype` note, and
    `159205730_TGTCC...>T` (a large indel, written with `>` rather than `_`).

    Args:
        token: One comma-separated token from a db genome column.

    Returns:
        `(pos, kind, ref, alt)`, or None for a token with no numeric position, which covers
        unsupported notations (`unknown_del`, `unknown_del_27kb`, `na`), which
        `parse_positions` also skips. An SV notation yields `kind='sv'` carrying its own
        two words, e.g. `('del', '21kb')`, which are not sequences.

    Raises:
        ValueError: The token has a numeric position but no parseable allele, meaning the
            db has grown a coordinate form this does not understand.
    """
    token = token.strip()
    pos, _, rest = token.partition('_')
    if not pos.isdigit():
        return None
    parts = rest.split('_')
    if parts[0] == 'ref':
        return (int(pos), 'ref', '.', '.')
    if parts[0].lower() in _SV_WORDS:
        return (int(pos), 'sv', parts[0], parts[1] if len(parts) > 1 else '.')
    if len(parts) >= 2:
        ref, alt = parts[0], parts[1]
    elif '>' in parts[0]:
        ref, _, alt = parts[0].partition('>')
    else:
        raise ValueError(f'Unparseable allele in db coordinate token: {token!r}')
    if not (_SEQUENCE.match(ref) and _SEQUENCE.match(alt)):
        raise ValueError(
            f'Unrecognised allele {ref}/{alt} in db coordinate token {token!r}. '
            'Teach parse_defining_token the new notation rather than shipping it as a variant.'
        )
    return (int(pos), 'var', ref, alt)


def iter_defining_sites(rows: Iterable[dict[str, str]], genome: str) -> Iterator[DefiningSite]:
    """Yield every allele-defining site in the db, with duplicates.

    Args:
        rows: Parsed db rows.
        genome: Coordinate column to read, `GRCh37` or `GRCh38`.

    Yields:
        One DefiningSite per coordinate token per row, so a multi-site allele such as
        `JK*02N.17` yields one site per defining SNP, and a site shared by several alleles
        of one system is yielded once per row. Rows with no coordinate are skipped. SV
        notations are yielded with `kind='sv'`, for the caller to count or exclude.
    """
    for row in rows:
        cell = (row[genome] or '').strip()
        if cell in _NO_COORDINATE:
            continue
        chrom = norm_chrom(row['Chrom'])
        system = blood_group_system(row.get('Genotype') or '', row.get('Genotype_alt') or '')
        for token in cell.split(','):
            parsed = parse_defining_token(token)
            if parsed is None:
                continue
            pos, kind, ref, alt = parsed
            yield DefiningSite(chrom=chrom, pos=pos, ref=ref, alt=alt, kind=kind, system=system)


def site_system_map(rows: Iterable[dict[str, str]], genome: str) -> list[DefiningSite]:
    """Build the site -> blood-group-system map.

    SV sites are excluded. One base's DP and GQ cannot say whether a sample carries a 21kb
    deletion, so assessing an SV-defined allele at its start coordinate would report a
    quality for something the check never looked at. A system whose only defining alleles
    are SVs therefore has no row here and is reported as not assessed. At v2.4.1 that is
    ABCC1, ATP11C and CD99.

    Args:
        rows: Parsed db rows.
        genome: Coordinate column to read, `GRCh37` or `GRCh38`.

    Returns:
        The deduplicated non-SV sites, sorted by coordinate then system then allele. A site
        can map to more than one system, because the RHD and RHCE coordinates overlap, so this is
        a list of rows rather than a `{site: system}` dict.
    """
    return sorted(
        {s for s in iter_defining_sites(rows, genome) if s.kind != 'sv'},
        key=lambda s: (chrom_key(s.chrom), s.pos, s.system, s.kind, s.ref, s.alt),
    )


def defining_site_positions(rows: Iterable[dict[str, str]], genome: str) -> list[tuple[str, int]]:
    """Collect the coordinates to extract DP/GQ at.

    Excludes SV sites, for the reason `site_system_map` gives: extracting a coordinate the
    map has no row for would only cost bandwidth.

    Args:
        rows: Parsed db rows.
        genome: Coordinate column to read, `GRCh37` or `GRCh38`.

    Returns:
        The distinct `(chrom, pos)` pairs, sorted.
    """
    return sorted(
        {(s.chrom, s.pos) for s in iter_defining_sites(rows, genome) if s.kind != 'sv'},
        key=lambda cp: (chrom_key(cp[0]), cp[1]),
    )


def build_intervals(rows: Iterable[dict[str, str]], genome: str, flank: int) -> dict[str, list[list[int]]]:
    """Build merged ±flank intervals, as `rbceq2.IO.vcf.build_intervals` builds them.

    Args:
        rows: Parsed db rows.
        genome: Coordinate column to read, `GRCh37` or `GRCh38`.
        flank: Bases to extend each defining position by on either side.

    Returns:
        Per-contig merged `[start, end]` intervals, sorted by start within each contig.
    """
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for row in rows:
        chrom = norm_chrom(row['Chrom'])
        for pos in parse_positions(row[genome]):
            intervals[chrom].append((max(0, pos - flank), pos + flank))
    merged: dict[str, list[list[int]]] = {}
    for chrom, spans in intervals.items():
        spans.sort()
        out: list[list[int]] = []
        for start, end in spans:
            if not out or start > out[-1][1]:
                out.append([start, end])
            else:
                out[-1][1] = max(out[-1][1], end)
        merged[chrom] = out
    return merged
