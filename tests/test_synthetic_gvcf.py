"""The fixture generator writes the ploidy it claims to, at sites the pipeline ships.

A fixture is wrong in two ways that look fine. It can put the ploidy boundary in the wrong
place, in which case the single-copy case is never exercised and the run looks clean. Or it
can call a coordinate rbceq2's database has nothing at, in which case every system comes back
reference and the diff is empty for the wrong reason.

So these tests check the boundary against PAR1 and PAR2, and check every called coordinate
against the same committed defining-sites BED the pipeline hands bcftools. No external
tools: this is about what the generator writes, not what bcftools then does with it.
"""

import csv
from pathlib import Path

import pytest

from popgen_rbceq2 import constants, stage_support
from popgen_rbceq2.scripts import bg_db, gen_synthetic_gvcf

GENOME = 'GRCh38'

# PAR1 ends at 2,781,479 and PAR2 begins at 155,701,383 on GRCh38. A male sample is diploid
# in both and haploid between them, which is the distinction the whole fixture rests on.
PAR1_LAST = 2_781_479
PAR2_FIRST = 155_701_383

pytestmark = pytest.mark.fast


def _records(text: str) -> list[list[str]]:
    return [line.split('\t') for line in text.splitlines() if not line.startswith('#')]


def _variants(text: str) -> list[list[str]]:
    """Variant records only, dropping the <NON_REF> reference blocks."""
    return [r for r in _records(text) if r[4] != '<NON_REF>']


def _gt(record: list[str]) -> str:
    return record[9].split(':')[0]


def _is_one_token(gt: str) -> bool:
    return '/' not in gt and '|' not in gt


def test_the_production_par_bounds_are_where_grch38_puts_them():
    """The constant the merge's gate reads, pinned against the literals.

    This is the copy that matters. `is_non_par_x` reads it rather than holding its own, so
    without this assertion the production bounds could be wrong by 100 Mb and every test below
    would still pass, because they would all be wrong together.
    """
    assert constants.NON_PAR_X[GENOME] == (PAR1_LAST + 1, PAR2_FIRST - 1)


@pytest.mark.parametrize(
    ('chrom', 'pos', 'expected'),
    [
        ('chrX', PAR1_LAST, False),
        ('chrX', PAR1_LAST + 1, True),
        ('chrX', PAR2_FIRST - 1, True),
        ('chrX', PAR2_FIRST, False),
        # A contig test dropped from the window would leave these single-copy, and the
        # autosomal half of every fixture would be quietly wrong.
        ('chr1', PAR1_LAST + 1, False),
        ('chr9', PAR1_LAST + 1, False),
    ],
)
def test_the_window_is_chrx_between_par1_and_par2(chrom, pos, expected):
    """The window is a contig and a range, starting one base past PAR1."""
    assert gen_synthetic_gvcf.is_non_par_x(chrom, pos, GENOME) is expected


def test_a_build_with_no_recorded_bounds_raises_rather_than_guessing():
    """GRCh37 moves both boundaries, so falling back to GRCh38's would be silently wrong."""
    # The merge's gate reads the same table, where a wrong window means filling single-copy
    # chrX or gating PAR. Refusing is the only safe answer for a build we ship no bounds for.
    with pytest.raises(KeyError, match='GRCh37'):
        constants.non_par_x('GRCh37')
    with pytest.raises(KeyError, match='GRCh37'):
        gen_synthetic_gvcf.is_non_par_x('chrX', PAR1_LAST + 1, 'GRCh37')


def test_the_native_fixture_writes_one_copy_on_chrx_and_two_in_par():
    """The default output is what DRAGEN writes for a single-copy sample.

    Named coordinates, not `is_non_par_x`. Classifying the fixture with the same function that
    built it cannot fail: moving the window to 40 Mb would leave that version green.
    """
    gts = {(r[0], int(r[1])): _gt(r) for r in _variants(gen_synthetic_gvcf.build(GENOME))}
    assert gts[('chrX', 37686082)] == '1'  # XK, between PAR1 and PAR2
    assert gts[('chrX', 2748343)] == '0/1'  # XG, inside PAR1
    assert gts[('chr1', 159205564)] == '1/1'  # FY, autosomal


def test_the_diploidised_fixture_has_no_haploid_genotype_left():
    """`--diploidise` leaves no single-copy genotype anywhere in the file."""
    text = gen_synthetic_gvcf.build(GENOME, diploidise=True)
    assert not any(_is_one_token(_gt(r)) for r in _records(text))
    # Allele duplication with a phased separator, which is the shape bcftools 1.24 writes.
    assert any(_gt(r) == '1|1' for r in _variants(text))


def test_the_mixed_fixture_disagrees_with_itself_on_one_contig():
    """`--mixed` is the half-rewritten state, and has to be mixed on non-PAR chrX itself."""
    # Mixing across contigs would not reproduce it: rbceq2 derives its copy count per
    # chromosome, so a haploid chrX beside a diploid chr1 is an ordinary male sample.
    native = _variants(gen_synthetic_gvcf.build(GENOME))
    mixed = _variants(gen_synthetic_gvcf.build(GENOME, mixed=True))
    outside = [_gt(r) for r in mixed if gen_synthetic_gvcf.is_non_par_x(r[0], int(r[1]), GENOME)]
    assert any(_is_one_token(gt) for gt in outside)
    assert any(not _is_one_token(gt) for gt in outside)
    # Adds to the native fixture rather than replacing it, so a diff of the two isolates
    # the contradiction rather than mixing it with a changed call somewhere else.
    assert {(r[0], r[1]) for r in native} < {(r[0], r[1]) for r in mixed}


@pytest.mark.parametrize('kwargs', [{}, {'diploidise': True}, {'mixed': True}])
def test_every_called_coordinate_is_a_committed_defining_site(kwargs):
    """A fixture that calls a site rbceq2 does not define would produce an empty diff."""
    bed = Path(stage_support.blood_group_resource(f'bg_defining_sites.{GENOME}.bed'))
    sites = {(row[0], int(row[2])) for row in csv.reader(bed.read_text().splitlines(), delimiter='\t') if row}

    for record in _variants(gen_synthetic_gvcf.build(GENOME, **kwargs)):
        assert (record[0], int(record[1])) in sites, f'{record[0]}:{record[1]} is not a defining site'


@pytest.mark.parametrize('kwargs', [{}, {'mixed': True}])
def test_the_records_are_sorted(kwargs):
    """An unsorted file is refused by bcftools index, and the stage indexes what it is given."""
    records = _records(gen_synthetic_gvcf.build(GENOME, **kwargs))
    keys = [(bg_db.chrom_key(r[0]), int(r[1])) for r in records]
    assert keys == sorted(keys)


def test_a_reference_block_precedes_each_call_and_does_not_cover_it():
    """The blocks exercise the <NON_REF> drop without masking the variant they sit before."""
    # A block whose END reached the call would make the site covered twice, and the QC
    # extract would read the block's MIN_DP where it should read the variant's DP.
    #
    # The comparison has to be against the variant's own POS. Comparing the block's END to
    # anything derived from BLOCK_LEN checks the generator against its own arithmetic, which
    # holds whatever the block and the variant do to each other.
    text = gen_synthetic_gvcf.build(GENOME)
    records = _records(text)
    blocks = [r for r in records if r[4] == '<NON_REF>']
    variants = _variants(text)
    assert len(blocks) == len(variants)

    # Each block is emitted immediately before the variant it belongs to, and the whole file
    # is then sorted by coordinate, so pairing by position is what the reader sees.
    for block in blocks:
        following = [v for v in variants if v[0] == block[0] and int(v[1]) > int(block[1])]
        assert following, f'block at {block[0]}:{block[1]} precedes no variant on its contig'
        variant_pos = min(int(v[1]) for v in following)
        end = int(block[7].removeprefix('END='))
        assert end < variant_pos, f'block {block[0]}:{block[1]}-{end} covers the call at {variant_pos}'


@pytest.mark.parametrize(
    ('calls', 'match'),
    [
        # A system that left the resources. A silent skip would shrink the fixture and
        # weaken every comparison built on it.
        ((('NOT_A_BLOOD_GROUP', 1, 1),), 'NOT_A_BLOOD_GROUP'),
        # A coordinate that left them, which is the case naming sites by coordinate exists to
        # catch. An index into a system's site list stays valid when the list is reshuffled
        # and simply points at a different allele, and every other test here would still pass.
        ((('XK', 1, 1),), 'XK has no variant site at GRCh38 position 1'),
    ],
)
def test_a_call_that_left_the_resources_fails_by_name(monkeypatch, calls, match):
    """`CALLS` is a hard-coded table, and a resource bump must not silently shrink it."""
    monkeypatch.setattr(gen_synthetic_gvcf, 'CALLS', calls)
    with pytest.raises(LookupError, match=match):
        gen_synthetic_gvcf.build(GENOME)
