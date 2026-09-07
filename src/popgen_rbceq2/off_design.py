"""The defining sites an exome capture design never targeted, as committed resources.

Post-hoc calls may fill a hole only at a defining site outside the cohort's capture design. That
set depends on two fixed inputs, the committed defining sites and the design's BED, and on
nothing about any sample, so it is computed once per design by `scripts/gen_off_design_sites.py`
and committed under `resources/` beside the sites it was subtracted from. A run resolves the
file for its configured design at graph build (`resource_path`) and localises it like any other
resource; nothing in the pipeline reads the vendor BED.

The manifest beside the BEDs records what each was built from. Its `sites_md5` is what ties a
committed subtraction to the sites BED it is a subset of: regenerating the sites without the
subtractions would otherwise leave a stale set in place, and the fill would quietly narrow or
widen. `resource_path` refuses a design whose row does not match the shipped sites, and the test
suite checks every committed row the same way.
"""

import csv
import dataclasses
import hashlib
from pathlib import Path

import cpg_utils.config

from popgen_rbceq2 import stage_support

GENERATOR = 'scripts/gen_off_design_sites.py'

# One manifest for every design and build, beside the BEDs it describes.
MANIFEST_NAME = 'bg_off_design_sites.manifest.tsv'


@dataclasses.dataclass(frozen=True)
class ManifestRow:
    """What one committed off-design BED was built from.

    Attributes:
        design_key: The `[references]` key naming the design, as configured in
            `workflow.exome_design_bed`.
        genome: Reference build of the sites BED subtracted from.
        resource: Filename of the off-design BED under `resources/`.
        design_path: Where the design BED was read from when the resource was generated.
        design_md5: MD5 of that design BED's bytes.
        sites_md5: MD5 of the `bg_defining_sites.<genome>.bed` subtracted from.
        n_sites: Rows in that sites BED.
        n_off_design: Rows in the resource: the sites outside the design.
        bedtools_version: The bedtools that did the subtraction.
        generated: ISO date the resource was written.
    """

    design_key: str
    genome: str
    resource: str
    design_path: str
    design_md5: str
    sites_md5: str
    n_sites: int
    n_off_design: int
    bedtools_version: str
    generated: str


MANIFEST_COLUMNS = tuple(field.name for field in dataclasses.fields(ManifestRow))


def resource_name(design_key: str, genome: str) -> str:
    """The filename under `resources/` of a design's off-design defining sites.

    Args:
        design_key: The `[references]` key naming the design.
        genome: Reference build of the sites subtracted from.

    Returns:
        `bg_off_design_sites.<design segment>.<genome>.bed`, the design rendered as
        `stage_support.design_segment` renders it into output paths.
    """
    return f'bg_off_design_sites.{stage_support.design_segment(design_key)}.{genome}.bed'


def file_md5(path: Path) -> str:
    """Hex MD5 of a file's bytes."""
    return hashlib.md5(path.read_bytes(), usedforsecurity=False).hexdigest()


def read_manifest(path: Path) -> dict[tuple[str, str], ManifestRow]:
    """Read the manifest, keyed by (design_key, genome).

    Raises:
        ValueError: The header is not exactly `MANIFEST_COLUMNS`, or two rows name one design
            and build.
    """
    with path.open(newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        if tuple(reader.fieldnames or ()) != MANIFEST_COLUMNS:
            raise ValueError(f'{path} has columns {reader.fieldnames}, expected {list(MANIFEST_COLUMNS)}')
        rows = {}
        for raw in reader:
            row = ManifestRow(**{**raw, 'n_sites': int(raw['n_sites']), 'n_off_design': int(raw['n_off_design'])})
            key = (row.design_key, row.genome)
            if key in rows:
                raise ValueError(f'{path} lists {row.design_key} on {row.genome} twice')
            rows[key] = row
    return rows


def write_manifest(path: Path, rows: dict[tuple[str, str], ManifestRow]) -> None:
    """Write the manifest, one row per design and build, in a stable order."""
    with path.open('w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=MANIFEST_COLUMNS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for key in sorted(rows):
            writer.writerow(dataclasses.asdict(rows[key]))


def resource_path(design_key: str) -> str:
    """The committed off-design defining sites for this design and the configured build.

    Resolved at graph build, so a run whose design has no committed subtraction, or whose
    subtraction predates the shipped sites, fails before any job starts and says what to run.

    Args:
        design_key: The `[references]` key naming the design, from `workflow.exome_design_bed`.

    Returns:
        The absolute path of the shipped BED.

    Raises:
        FileNotFoundError: No resource, or no manifest row, is committed for the design.
        ValueError: The manifest row was built from a different sites BED than the one shipped,
            or records no off-design site, so there is nothing for the recall to fill.
    """
    genome = cpg_utils.config.genome_build()
    regenerate = f'Generate it with {GENERATOR} against that design and commit it under resources/.'
    path = stage_support.blood_group_resource(
        resource_name(design_key, genome),
        missing_hint=(
            f'No off-design defining sites are committed for {stage_support.DESIGN_CONFIG_PATH} = {design_key} '
            f'on {genome}. {regenerate}'
        ),
    )
    manifest = read_manifest(Path(stage_support.blood_group_resource(MANIFEST_NAME)))
    row = manifest.get((design_key, genome))
    if row is None:
        raise FileNotFoundError(
            f'{MANIFEST_NAME} has no row for {design_key} on {genome}, so {Path(path).name} cannot be checked '
            f'against the shipped defining sites. {regenerate}'
        )
    sites_md5 = file_md5(Path(stage_support.blood_group_resource(f'bg_defining_sites.{genome}.bed')))
    if row.sites_md5 != sites_md5:
        raise ValueError(
            f'{Path(path).name} was subtracted from a bg_defining_sites.{genome}.bed with MD5 {row.sites_md5}, '
            f'but the shipped one is {sites_md5}: the defining sites were regenerated without it. {regenerate}'
        )
    if row.n_off_design == 0:
        # The generator refuses to write this, so reaching it means a hand-edited manifest or
        # resource. Refused here too because the merge hands the file to `bcftools -T`, which
        # aborts on an empty targets file with a message naming a temp path and nothing about
        # the design, once per exome sequencing group.
        raise ValueError(
            f'{Path(path).name} is empty: {design_key} targets every defining site on {genome}, so post-hoc '
            'calling has nothing to fill. Run this cohort without the recall: leave '
            f'{stage_support.DESIGN_CONFIG_PATH} unset and keep the post-hoc stages out of only_stages.'
        )
    return path


def subtraction_commands(sites_bed: str, design_bed: str, out_bed: str, design_key: str, design_path: str) -> str:
    """Shell that writes the defining sites outside the design, or fails on a design that cannot be right.

    A pure function of its arguments so the test suite can run it under real bedtools; the
    subtraction is bedtools' interval semantics and nothing of ours, and a test that spelled
    the command out again would stay green if this one changed.

    An empty *result* is a valid subtraction, and the generator, not this shell, refuses to
    commit it: a design that targets every defining site leaves nothing to fill, and the
    right run for such a cohort is one without the recall configured, not one that localises
    an empty targets file to every exome (`bcftools -T` aborts on one). Two inputs are not
    valid, and both would
    otherwise surface the same way, as every defining site off-design, which the per-sample
    check on DRAGEN's records then fails on one sample at a time, naming the design file rather
    than what was wrong with it. Both are caught here instead, once, before any of that:

    - an empty design BED;
    - a design that shares no interval with any site, which for a real exome design means the
      two BEDs name their contigs differently (`1` against `chr1`). bedtools warns about that
      on stderr and exits 0 with every site off-design.

    A third input is refused for the opposite reason, that it would narrow the fill silently:
    a design row whose end is not greater than its start. BED is half-open, so such a row
    describes no bases, but `bedtools intersect` treats it as covering the base at its
    coordinate and the one before (probed on 2.31.1), so a defining site there would count as
    in-design and stay NOCOV instead of being filled. The awk loop this replaced ignored such
    rows. Neither real design has one (0 of 229,273 Twist rows, 0 of 275,017 CREv2 rows), so
    rather than pick a meaning for a malformed row the subtraction fails naming it: a design
    with one is not the file DRAGEN was given, or has been mangled since.

    Args:
        sites_bed: Path of `bg_defining_sites.<genome>.bed`.
        design_bed: Path of the capture design BED.
        out_bed: Where to write the off-design sites.
        design_key: The `[references]` key the design came from, for the messages.
        design_path: The path it was read from, for the messages.

    Returns:
        Shell lines, for `bash -c`.
    """
    return f"""
            if [ ! -s {design_bed} ]; then
                echo "ERROR: the capture design BED is empty." >&2
                echo "design {design_key}" >&2
                echo "read from {design_path}" >&2
                exit 1
            fi
            awk -F'\\t' '!/^(track|browser|#)/ && $3 <= $2' {design_bed} > zero_length_rows.bed
            if [ -s zero_length_rows.bed ]; then
                n_zero=$(wc -l < zero_length_rows.bed | tr -d ' ')
                echo "ERROR: $n_zero capture design row(s) have end <= start, so describe no bases." >&2
                echo "bedtools would count each as covering a base, narrowing the fill silently. First rows:" >&2
                head -n 5 zero_length_rows.bed >&2
                echo "design {design_key}" >&2
                echo "read from {design_path}" >&2
                exit 1
            fi
            bedtools intersect -v -a {sites_bed} -b {design_bed} > {out_bed}
            # awk's NR, not `wc -l`: wc counts newlines, so a file whose last line lacks one
            # comes up one short and the two counts below could disagree for a reason that is
            # nothing to do with the data.
            n_off=$(awk 'END {{print NR}}' {out_bed})
            n_sites=$(awk 'END {{print NR}}' {sites_bed})
            if [ "$n_off" -eq "$n_sites" ]; then
                echo "ERROR: every defining site is outside the capture design, which no exome design leaves." >&2
                echo "The two BEDs most likely name their contigs differently:" >&2
                echo "  defining sites: $(cut -f1 {sites_bed} | sort -u | tr '\\n' ' ')" >&2
                design_contigs=$(awk '!/^(track|browser|#)/ {{print $1}}' {design_bed} | sort -u | tr '\\n' ' ')
                echo "  design:         $design_contigs" >&2
                echo "design {design_key}" >&2
                echo "read from {design_path}" >&2
                exit 1
            fi
            echo "off-design: $n_off of $n_sites defining site(s) are outside {design_key}," >&2
            echo "and so are the only sites a post-hoc record may fill" >&2
    """
