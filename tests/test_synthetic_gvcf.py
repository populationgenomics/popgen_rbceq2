"""The fixture generator writes the ploidy it claims to, at sites the pipeline ships.

`gen_synthetic_gvcf.py` exists to answer one question locally: what do the rbceq2 TSVs say
when non-PAR chrX is encoded the way DRAGEN encodes it. That answer is only worth anything
if the fixture is right, and a fixture is wrong in two ways that look fine. It can put the
ploidy boundary in the wrong place, in which case the haploid case is never exercised and
the run looks clean. Or it can call a coordinate rbceq2's database has nothing at, in which
case every system comes back reference and the diff is empty for the wrong reason.

So these tests check the boundary against PAR1 and PAR2, and check every called coordinate
against the same committed defining-sites BED the pipeline hands bcftools. No external
tools: this is about what the generator writes, not what bcftools then does with it.
"""

import csv
from pathlib import Path

import pytest

from popgen_rbceq2 import stage_support
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


def test_the_par_boundary_is_where_grch38_puts_it():
    """The haploid window starts one base past PAR1 and ends one base before PAR2."""
    assert not gen_synthetic_gvcf.is_non_par_x('chrX', PAR1_LAST)
    assert gen_synthetic_gvcf.is_non_par_x('chrX', PAR1_LAST + 1)
    assert gen_synthetic_gvcf.is_non_par_x('chrX', PAR2_FIRST - 1)
    assert not gen_synthetic_gvcf.is_non_par_x('chrX', PAR2_FIRST)


def test_no_autosome_is_ever_treated_as_haploid():
    """Only chrX carries the one-token encoding; chr1 at the same offsets does not."""
    # The window is a coordinate range, so a contig check that was dropped would leave
    # chr1:37686082 haploid and the autosomal half of every fixture quietly wrong.
    assert not gen_synthetic_gvcf.is_non_par_x('chr1', PAR1_LAST + 1)
    assert not gen_synthetic_gvcf.is_non_par_x('chr9', PAR1_LAST + 1)


def test_the_native_fixture_is_haploid_outside_par_and_diploid_inside_it():
    """The default output is what DRAGEN writes for a male sample."""
    text = gen_synthetic_gvcf.build(GENOME)
    outside = [r for r in _variants(text) if gen_synthetic_gvcf.is_non_par_x(r[0], int(r[1]))]
    inside = [r for r in _variants(text) if not gen_synthetic_gvcf.is_non_par_x(r[0], int(r[1]))]

    assert outside, 'the fixture must exercise the non-PAR chrX path or it tests nothing'
    assert inside, 'the fixture must keep a diploid comparison in the same file'
    assert all(_is_one_token(_gt(r)) for r in outside)
    assert not any(_is_one_token(_gt(r)) for r in inside)


def test_the_diploidised_fixture_has_no_haploid_genotype_left():
    """`--diploidise` reproduces what the removed +fixploidy step left behind."""
    text = gen_synthetic_gvcf.build(GENOME, diploidise=True)
    assert not any(_is_one_token(_gt(r)) for r in _records(text))
    # Allele duplication with a phased separator is what bcftools 1.24 actually writes.
    assert any(_gt(r) == '1|1' for r in _variants(text))


def test_the_mixed_fixture_disagrees_with_itself_on_one_contig():
    """`--mixed` is the half-rewritten state, and has to be mixed on non-PAR chrX itself."""
    # Mixing across contigs would not reproduce it: rbceq2 derives its copy count per
    # chromosome, so a haploid chrX beside a diploid chr1 is an ordinary male sample.
    text = gen_synthetic_gvcf.build(GENOME, mixed=True)
    outside = [_gt(r) for r in _variants(text) if gen_synthetic_gvcf.is_non_par_x(r[0], int(r[1]))]
    assert any(_is_one_token(gt) for gt in outside)
    assert any(not _is_one_token(gt) for gt in outside)


def test_the_mixed_fixture_adds_to_the_native_one_rather_than_replacing_it():
    """Everything the native fixture calls is still called, so the diff isolates the mix."""
    native = {(r[0], r[1]) for r in _variants(gen_synthetic_gvcf.build(GENOME))}
    mixed = {(r[0], r[1]) for r in _variants(gen_synthetic_gvcf.build(GENOME, mixed=True))}
    assert native < mixed


@pytest.mark.parametrize('kwargs', [{}, {'diploidise': True}, {'mixed': True}])
def test_every_called_coordinate_is_a_committed_defining_site(kwargs):
    """A fixture that calls a site rbceq2 does not define would produce an empty diff."""
    bed = Path(stage_support.blood_group_resource(f'bg_defining_sites.{GENOME}.bed'))
    sites = {(row[0], int(row[2])) for row in csv.reader(bed.read_text().splitlines(), delimiter='\t') if row}

    for record in _variants(gen_synthetic_gvcf.build(GENOME, **kwargs)):
        assert (record[0], int(record[1])) in sites, f'{record[0]}:{record[1]} is not a defining site'


@pytest.mark.parametrize('kwargs', [{}, {'diploidise': True}, {'mixed': True}])
def test_the_records_are_sorted(kwargs):
    """An unsorted file is refused by bcftools index, and the stage indexes what it is given."""
    records = _records(gen_synthetic_gvcf.build(GENOME, **kwargs))
    keys = [(bg_db.chrom_key(r[0]), int(r[1])) for r in records]
    assert keys == sorted(keys)


def test_a_reference_block_precedes_each_call_and_does_not_cover_it():
    """The blocks exercise the <NON_REF> drop without masking the variant they sit before."""
    # A block whose END reached the call would make the site covered twice, and the QC
    # extract would read the block's MIN_DP where it should read the variant's DP.
    text = gen_synthetic_gvcf.build(GENOME)
    blocks = [r for r in _records(text) if r[4] == '<NON_REF>']
    assert len(blocks) == len(_variants(text))
    for block in blocks:
        end = int(block[7].removeprefix('END='))
        assert end < int(block[1]) + gen_synthetic_gvcf.BLOCK_LEN


def test_an_unknown_system_fails_loudly():
    """A resource bump that removes a system this script calls must not pass silently."""
    # The indices in CALLS point into a system's site list. If the system disappeared, a
    # silent skip would shrink the fixture and weaken every comparison built on it.
    original = gen_synthetic_gvcf.CALLS
    gen_synthetic_gvcf.CALLS = (('NOT_A_BLOOD_GROUP', 0, 1),)
    try:
        with pytest.raises(KeyError, match='NOT_A_BLOOD_GROUP'):
            gen_synthetic_gvcf.build(GENOME)
    finally:
        gen_synthetic_gvcf.CALLS = original
