#!/usr/bin/env python3
"""Generate the committed off-design defining sites for one exome capture design.

Run by hand when a cohort arrives on a design with no committed subtraction, or after
`gen_bg_resources.py` regenerates the defining sites, and commit the results:

    gen_off_design_sites.py <design_key> <design_bed> GRCh38 src/popgen_rbceq2/resources

`design_key` is the `[references]` key a run sets as `workflow.exome_design_bed`, and
`design_bed` is the BED it names in the references repo, as a `gs://` or local path. The README
section "An exome run must name its capture design" lists both for the designs seen so far.

Writes `bg_off_design_sites.<design>.<genome>.bed`, the defining sites bedtools finds outside
the design, and adds or replaces the design's row in `bg_off_design_sites.manifest.tsv`, which
records the design and sites files it was built from. `off_design.resource_path` reads that row
at graph build and refuses a resource whose sites have since been regenerated.

Needs `bedtools` on PATH and, for a `gs://` design, application-default credentials that can
read it.
"""

import argparse
import datetime
import hashlib
import logging
import subprocess
import sys
import tempfile
from pathlib import Path

from cpg_utils import to_path

from popgen_rbceq2 import off_design

logger = logging.getLogger(__name__)


def bedtools_version() -> str:
    """The `bedtools --version` string, e.g. `bedtools v2.31.1`."""
    return subprocess.run(['bedtools', '--version'], capture_output=True, text=True, check=True).stdout.strip()  # noqa: S607


def generate(design_key: str, design_path: str, genome: str, out_dir: Path) -> off_design.ManifestRow:
    """Subtract one design from the committed sites and write the resource and its manifest row.

    Args:
        design_key: The `[references]` key naming the design.
        design_path: The design BED, `gs://` or local.
        genome: Reference build; selects `bg_defining_sites.<genome>.bed` in `out_dir`.
        out_dir: The resources directory, holding the sites BED and receiving the output.

    Returns:
        The manifest row written.

    Raises:
        FileNotFoundError: `out_dir` has no sites BED for the build.
        RuntimeError: bedtools, or one of the design checks around it, failed; stderr is in
            the message.
    """
    # Absolute, because the shell below runs in a scratch directory.
    sites_bed = (out_dir / f'bg_defining_sites.{genome}.bed').resolve()
    if not sites_bed.is_file():
        raise FileNotFoundError(f'{sites_bed} is not there; generate it with gen_bg_resources.py first')
    design_bytes = to_path(design_path).read_bytes()
    resource = out_dir / off_design.resource_name(design_key, genome)

    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        local_design = work / 'design.bed'
        local_design.write_bytes(design_bytes)
        out_bed = work / 'off_design.bed'
        script = 'set -euo pipefail\n' + off_design.subtraction_commands(
            str(sites_bed),
            str(local_design),
            str(out_bed),
            design_key,
            design_path,
        )
        result = subprocess.run(['bash', '-c', script], cwd=work, capture_output=True, text=True, check=False)  # noqa: S603, S607
        if result.returncode != 0:
            raise RuntimeError(f'the subtraction failed for {design_key}:\n{result.stderr}')
        logger.info(result.stderr.strip())
        resource.write_bytes(out_bed.read_bytes())

    row = off_design.ManifestRow(
        design_key=design_key,
        genome=genome,
        resource=resource.name,
        design_path=design_path,
        design_md5=hashlib.md5(design_bytes, usedforsecurity=False).hexdigest(),
        sites_md5=off_design.file_md5(sites_bed),
        n_sites=sum(1 for _ in sites_bed.open()),
        n_off_design=sum(1 for _ in resource.open()),
        bedtools_version=bedtools_version(),
        generated=datetime.date.today().isoformat(),  # noqa: DTZ011
    )
    manifest_path = out_dir / off_design.MANIFEST_NAME
    rows = off_design.read_manifest(manifest_path) if manifest_path.is_file() else {}
    rows[(design_key, genome)] = row
    off_design.write_manifest(manifest_path, rows)
    return row


def main(argv: list[str] | None = None) -> None:
    """Generate the resource for one design.

    Args:
        argv: Command-line arguments; defaults to `sys.argv[1:]`.
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('design_key', help="the design's [references] key, as set in workflow.exome_design_bed")
    parser.add_argument('design_bed', help='the capture design BED that key names, gs:// or local')
    parser.add_argument('genome', choices=['GRCh37', 'GRCh38'], help='build of the committed sites to subtract from')
    parser.add_argument('out_dir', type=Path, help='resources directory holding the sites BED; receives the output')
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s', stream=sys.stderr)
    row = generate(args.design_key, args.design_bed, args.genome, args.out_dir)
    logger.info(
        f'Wrote {row.resource}: {row.n_off_design} of {row.n_sites} defining sites are outside '
        f'{row.design_key} (design MD5 {row.design_md5}); manifest updated'
    )


if __name__ == '__main__':
    main()
