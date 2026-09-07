"""bedtools' interval semantics, where the committed off-design sites depend on them.

`scripts/gen_off_design_sites.py` is one `bedtools intersect -v` over a vendor capture BED, so
what has to hold is not our arithmetic but bedtools' agreement with the BED convention the
committed defining sites are written in. Three properties carry the resources it writes:

- **half-open intervals.** A design row `chr1 999 1100` targets 1-based bases 1000 to 1100. A
  site one base outside either end is off-design and eligible to be filled; the two ends are
  not. An off-by-one here either fills a site the capture did target, which answers a question
  about that sample's DRAGEN run with a second caller, or leaves a real hole unfilled.
- **`track` and `browser` lines are skipped.** Vendor BEDs carry them, and a header read as an
  interval would land on a contig no site is on, or worse, look like a target.
- **columns past the third are ignored.** Vendor BEDs carry a name, score and strand.

These were the awk's properties before the subtraction moved, and they are re-asserted here
against the tool that replaced it rather than assumed to carry over.

Skipped where bedtools is not installed; CI installs it. Every case here was also run against
`cpg-common/images/bedtools:2.30.0-1` and passes there, and the subtraction of the real Twist
design was diffed identical to the awk's over all 1,625 committed sites. The committed
resources themselves, and the manifest that ties them to the sites BED, are covered in
test_off_design_resources.
"""

import shutil
import subprocess
from pathlib import Path

import pytest

from popgen_rbceq2.off_design import subtraction_commands

pytestmark = [
    pytest.mark.fast,
    pytest.mark.skipif(not shutil.which('bedtools'), reason='bedtools is not on PATH'),
]


def _sites_bed(sites: list[tuple[str, int]]) -> str:
    """Render 1-based defining sites as the committed BED does, 0-based and half-open."""
    return ''.join(f'{chrom}\t{pos - 1}\t{pos}\n' for chrom, pos in sites)


def _subtract(tmp_path: Path, design: str, sites: list[tuple[str, int]]) -> subprocess.CompletedProcess[str]:
    """Run the generator's own subtraction command, verbatim, under real bedtools."""
    design_bed = tmp_path / 'design.bed'
    design_bed.write_text(design)
    sites_bed = tmp_path / 'sites.bed'
    sites_bed.write_text(_sites_bed(sites))
    script = 'set -euo pipefail\n' + subtraction_commands(
        str(sites_bed),
        str(design_bed),
        str(tmp_path / 'off_design.bed'),
        'exome_probesets_hg38/test_design_bed',
        'gs://bucket/test_design.bed',
    )
    return subprocess.run(['bash', '-c', script], cwd=tmp_path, capture_output=True, text=True, check=False)  # noqa: S603, S607


def _off_design(tmp_path: Path, design: str, sites: list[tuple[str, int]]) -> list[tuple[str, int]]:
    """Run the subtraction and return the sites it reports as off-design."""
    result = _subtract(tmp_path, design, sites)
    assert result.returncode == 0, result.stderr
    lines = (tmp_path / 'off_design.bed').read_text().splitlines()
    return [(line.split('\t')[0], int(line.split('\t')[2])) for line in lines]


@pytest.mark.parametrize('pos', [1000, 1100])
def test_both_ends_of_a_design_interval_are_targeted(tmp_path, pos):
    assert _off_design(tmp_path, 'chr1\t999\t1100\n', [('chr1', pos)]) == []


@pytest.mark.parametrize('pos', [999, 1101])
def test_one_base_outside_a_design_interval_is_off_design(tmp_path, pos):
    # A second, targeted site alongside, because a set that is entirely off-design is the
    # contig-mismatch failure below, and the command refuses it.
    assert _off_design(tmp_path, 'chr1\t999\t1100\n', [('chr1', 1050), ('chr1', pos)]) == [('chr1', pos)]


def test_a_site_on_a_contig_the_design_never_mentions_is_off_design(tmp_path):
    assert _off_design(tmp_path, 'chr1\t999\t1100\n', [('chr1', 1050), ('chr2', 1050)]) == [('chr2', 1050)]


def test_a_design_that_targets_every_site_leaves_nothing_to_fill(tmp_path):
    # A valid subtraction, and the opposite of the two failures below. The generator refuses
    # to commit it (test_off_design_resources), because a run on such a design has nothing to
    # fill and should not configure the recall.
    assert _off_design(tmp_path, 'chr1\t0\t5000\n', [('chr1', 1050), ('chr1', 2000)]) == []


def test_an_empty_design_bed_fails_the_subtraction(tmp_path):
    result = _subtract(tmp_path, '', [('chr1', 1050)])

    assert result.returncode == 1
    assert 'the capture design BED is empty' in result.stderr


def test_a_design_naming_its_contigs_differently_fails_the_subtraction(tmp_path):
    # `1` against `chr1`: bedtools warns on stderr and exits 0 with every site off-design, the
    # same observable as an empty design. Left alone it fails one sample at a time in the
    # conversion job, with a message blaming the design file rather than its contig names.
    result = _subtract(tmp_path, '1\t0\t5000\n2\t0\t5000\n', [('chr1', 1050), ('chr2', 2000)])

    assert result.returncode == 1
    assert 'every defining site is outside the capture design' in result.stderr
    assert 'defining sites: chr1 chr2' in result.stderr
    assert 'design:         1 2' in result.stderr


def test_a_design_bed_with_extra_columns_and_a_track_line_is_read_as_intervals(tmp_path):
    design = 'track name="Covered" description="probe footprint"\nchr1\t999\t1100\tTARGET_1\t0\t+\n'

    assert _off_design(tmp_path, design, [('chr1', 1050), ('chr1', 2000)]) == [('chr1', 2000)]


def test_the_reported_sites_keep_the_committed_beds_own_rows(tmp_path):
    # `-v` emits rows from -a untouched, so the output is a subset of the committed sites BED
    # and stays readable by everything downstream that reads that file's format.
    sites = [('chr1', 500), ('chr1', 1050), ('chr2', 700)]

    assert _off_design(tmp_path, 'chr1\t999\t1100\n', sites) == [('chr1', 500), ('chr2', 700)]


def test_a_design_row_with_end_not_greater_than_start_fails_the_subtraction(tmp_path):
    # BED is half-open, so such a row describes no bases, but bedtools treats it as covering
    # the base at its coordinate and the one before: probed here, the zero-length row at 1005
    # would otherwise make the site at 1005 in-design and leave its hole NOCOV. The awk this
    # replaced ignored such rows. Neither real design has one, so the generator refuses the
    # file rather than choosing a meaning for a malformed row.
    design = 'chr1\t999\t1100\nchr1\t1004\t1004\nchr1\t3000\t2990\n'

    result = _subtract(tmp_path, design, [('chr1', 1050), ('chr1', 1005), ('chr1', 2000)])

    assert result.returncode == 1
    assert '2 capture design row(s) have end <= start' in result.stderr
    assert 'chr1\t1004\t1004' in result.stderr
    assert 'chr1\t3000\t2990' in result.stderr
    assert not (tmp_path / 'off_design.bed').exists()


def test_a_track_line_is_not_mistaken_for_a_zero_length_row(tmp_path):
    # The header line has no numeric columns, so `$3 <= $2` compares strings there; it must be
    # excluded rather than reported.
    design = 'track name="Covered"\nchr1\t999\t1100\n'

    assert _off_design(tmp_path, design, [('chr1', 1050), ('chr1', 2000)]) == [('chr1', 2000)]
