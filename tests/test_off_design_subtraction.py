"""bedtools' interval semantics, where SelectOffDesignDefiningSites depends on them.

The stage is one `bedtools intersect -v` over a vendor capture BED, so what has to hold is not
our arithmetic but bedtools' agreement with the BED convention the committed defining sites are
written in. Three properties carry the stage:

- **half-open intervals.** A design row `chr1 999 1100` targets 1-based bases 1000 to 1100. A
  site one base outside either end is off-design and eligible to be filled; the two ends are
  not. An off-by-one here either fills a site the capture did target, which answers a question
  about that sample's DRAGEN run with a second caller, or leaves a real hole unfilled.
- **`track` and `browser` lines are skipped.** Vendor BEDs carry them, and a header read as an
  interval would land on a contig no site is on, or worse, look like a target.
- **columns past the third are ignored.** Vendor BEDs carry a name, score and strand.

These were the awk's properties before the subtraction moved, and they are re-asserted here
against the tool that replaced it rather than assumed to carry over.

Skipped where bedtools is not installed, which includes CI: the pinned image carries it and the
driver image does not. Every case here was also run against the pinned
`cpg-common/images/bedtools:2.30.0-1` and passes there, and its subtraction of the real Twist
design was diffed identical to the awk's over all 1,625 committed sites.
"""

import shutil
import subprocess
from pathlib import Path

import pytest

pytestmark = [
    pytest.mark.fast,
    pytest.mark.skipif(not shutil.which('bedtools'), reason='bedtools is not on PATH'),
]


def _sites_bed(sites: list[tuple[str, int]]) -> str:
    """Render 1-based defining sites as the committed BED does, 0-based and half-open."""
    return ''.join(f'{chrom}\t{pos - 1}\t{pos}\n' for chrom, pos in sites)


def _off_design(tmp_path: Path, design: str, sites: list[tuple[str, int]]) -> list[tuple[str, int]]:
    """Run the stage's subtraction and return the sites it reports as off-design."""
    design_bed = tmp_path / 'design.bed'
    design_bed.write_text(design)
    sites_bed = tmp_path / 'sites.bed'
    sites_bed.write_text(_sites_bed(sites))
    out = subprocess.run(  # noqa: S603
        ['bedtools', 'intersect', '-v', '-a', str(sites_bed), '-b', str(design_bed)],  # noqa: S607
        capture_output=True,
        text=True,
        check=True,
    )
    return [(line.split('\t')[0], int(line.split('\t')[2])) for line in out.stdout.splitlines()]


@pytest.mark.parametrize('pos', [1000, 1100])
def test_both_ends_of_a_design_interval_are_targeted(tmp_path, pos):
    assert _off_design(tmp_path, 'chr1\t999\t1100\n', [('chr1', pos)]) == []


@pytest.mark.parametrize('pos', [999, 1101])
def test_one_base_outside_a_design_interval_is_off_design(tmp_path, pos):
    assert _off_design(tmp_path, 'chr1\t999\t1100\n', [('chr1', pos)]) == [('chr1', pos)]


def test_a_site_on_a_contig_the_design_never_mentions_is_off_design(tmp_path):
    assert _off_design(tmp_path, 'chr1\t999\t1100\n', [('chr2', 1050)]) == [('chr2', 1050)]


def test_a_design_bed_with_extra_columns_and_a_track_line_is_read_as_intervals(tmp_path):
    design = 'track name="Covered" description="probe footprint"\nchr1\t999\t1100\tTARGET_1\t0\t+\n'

    assert _off_design(tmp_path, design, [('chr1', 1050), ('chr1', 2000)]) == [('chr1', 2000)]


def test_the_reported_sites_keep_the_committed_beds_own_rows(tmp_path):
    # `-v` emits rows from -a untouched, so the output is a subset of the committed sites BED
    # and stays readable by everything downstream that reads that file's format.
    sites = [('chr1', 500), ('chr1', 1050), ('chr2', 700)]

    assert _off_design(tmp_path, 'chr1\t999\t1100\n', sites) == [('chr1', 500), ('chr2', 700)]
