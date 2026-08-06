"""Tests for the blood-group call QC annotation.

Covers the db parse (`scripts/bg_db`), the committed resources it generates, and
the coverage resolution and flagging in `jobs/rbceq2_call_qc_job`.
"""

import re
from importlib.resources import files

import pytest

from popgen_rbceq2.jobs.rbceq2_call_qc_job import (
    DELETED,
    LOWQ,
    NOCOV,
    NOT_ASSESSED,
    PASS,
    Coverage,
    GvcfRecord,
    build_qc_tsv,
    flag_site,
    flags_by_system,
    load_site_systems,
    parse_extract,
    resolve_coverage,
)
from popgen_rbceq2.scripts import bg_db
from popgen_rbceq2.scripts.bg_db import DefiningSite
from popgen_rbceq2.stage_support import blood_group_resource

pytestmark = pytest.mark.fast

MIN_DEPTH = 10
MIN_GQ = 20

# GRCh38 primary-assembly contig lengths, to assert every shipped coordinate exists. The
# v2.4.1 db has ABCC1*01N.01 on chr16 carrying ABCC4's chr13 coordinate verbatim, which is
# past the end of chr16 in both builds.
GRCH38_CONTIG_LENGTHS = {
    'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555,
    'chr5': 181538259, 'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636,
    'chr9': 138394717, 'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
    'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 'chr16': 90338345,
    'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
    'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415,
}  # fmt: skip


def _site(chrom='chr1', pos=3774964, ref='A', alt='G', kind: bg_db.SiteKind = 'var', system='VEL') -> DefiningSite:
    return DefiningSite(chrom=chrom, pos=pos, ref=ref, alt=alt, kind=kind, system=system)


def _record(chrom, pos, ref, alt, end=None, gt='0/1', dp=None, gq=None, min_dp=None) -> GvcfRecord:
    """Build one extracted GVCF record, defaulting the fields a test is not exercising."""
    return GvcfRecord(chrom=chrom, pos=pos, ref=ref, alt=alt, end=end, gt=gt, dp=dp, gq=gq, min_dp=min_dp)


def _db_row(chrom='chr1', genotype='VEL*01.01', genome='3774964_A_G', genotype_alt='') -> dict[str, str]:
    return {'Chrom': chrom, 'Genotype': genotype, 'Genotype_alt': genotype_alt, 'GRCh38': genome}


# --- db parse -------------------------------------------------------------------------


@pytest.mark.parametrize(
    ('token', 'expected'),
    [
        ('3774964_A_G', (3774964, 'var', 'A', 'G')),
        # A lane site: defined by the reference base being present, no ref/alt of its own.
        ('207331122_ref', (207331122, 'ref', '.', '.')),
        # Trailing db notes must not be mistaken for alleles.
        ('25284544_G_C_no_phenotype', (25284544, 'var', 'G', 'C')),
        ('25321858_ref_no_phenotype', (25321858, 'ref', '.', '.')),
        # A large indel, which the db writes with `>` rather than `_`.
        ('159205730_TGTCCTGGCACAGCTG>T', (159205730, 'var', 'TGTCCTGGCACAGCTG', 'T')),
        # Structural variants, in the db's two spellings: `del`/`ins` with a size, and
        # RHD/RHCE's `DEL`/`INS` with a base count. Neither word is a sequence.
        ('95018451_del_21kb', (95018451, 'sv', 'del', '21kb')),
        ('25272546_DEL_148', (25272546, 'sv', 'DEL', '148')),
        ('144001284_ins_83bp', (144001284, 'sv', 'ins', '83bp')),
        # RBCeq2's own unsupported notations; parse_positions skips these too, so they are
        # absent from the regions BED and must be absent from the map.
        ('unknown_del', None),
        ('unknown_del_27kb', None),
        ('na', None),
    ],
)
def test_parse_defining_token(token, expected):
    assert bg_db.parse_defining_token(token) == expected


def test_parse_defining_token_raises_on_unknown_allele_form():
    """A numeric position with an unrecognised allele form raises."""
    # Staying silent would drop the site from the QC map while keeping it in the regions
    # BED, so the site would be converted but never assessed.
    with pytest.raises(ValueError, match='Unparseable allele'):
        bg_db.parse_defining_token('3774964_somethingnew')


def test_parse_defining_token_raises_on_an_allele_that_is_not_a_sequence():
    """A two-word token whose halves are not sequences raises rather than shipping."""
    # This is how `95018451_del_21kb` used to reach the map as ref='del', alt='21kb': a
    # notation read as an allele renders as `del>21kb` in a flag and gets a DP/GQ reported
    # for something the check never assessed.
    with pytest.raises(ValueError, match='Unrecognised allele'):
        bg_db.parse_defining_token('3774964_A_Z')


def test_structural_variant_sites_are_excluded_from_the_shipped_resources():
    """An SV-defined allele contributes no site to the map or the extraction BED."""
    # One base of DP/GQ cannot say whether a sample carries a 21kb deletion, so assessing
    # the SV's start coordinate would report a quality for something never looked at.
    rows = [
        _db_row(chrom='chr16', genotype='ABCC1*01N.01', genome='95018451_del_21kb'),
        _db_row(chrom='chr13', genotype='ABCC4*01N.02', genome='95163161_C_T'),
    ]
    assert [(s.chrom, s.pos, s.system) for s in bg_db.site_system_map(rows, 'GRCh38')] == [
        ('chr13', 95163161, 'ABCC4'),
    ]
    assert bg_db.defining_site_positions(rows, 'GRCh38') == [('chr13', 95163161)]


def test_structural_variant_positions_stay_in_the_regions_bed():
    """Excluding an SV site from the map does not shrink the regions BED."""
    # build_intervals mirrors rbceq2's own parse_positions, which keeps these coordinates.
    # The BED has to stay a superset of everything rbceq2 reads or calls go silently wrong.
    rows = [_db_row(chrom='chr16', genotype='ABCC1*01N.01', genome='95018451_del_21kb')]
    assert bg_db.build_intervals(rows, 'GRCh38', flank=100) == {'chr16': [[95018351, 95018551]]}


@pytest.mark.parametrize(
    ('genotype', 'genotype_alt', 'expected'),
    [
        ('VEL*01.01', '', 'VEL'),
        ('JK*02N.17', '', 'JK'),
        # RBCeq2 keys these alleles as KLF1 but emits the system as KLF.
        ('KLF1*BGM12', '', 'KLF'),
        ('KLF1*01', '', 'KLF'),
        # Falls back to Genotype_alt when Genotype is absent, as Allele.blood_group does.
        ('.', 'CR*A', 'CR'),
        ('', 'CR*A', 'CR'),
    ],
)
def test_blood_group_system(genotype, genotype_alt, expected):
    assert bg_db.blood_group_system(genotype, genotype_alt) == expected


def test_multi_site_allele_yields_one_site_per_defining_snp():
    """An allele with several defining coordinates yields one site each."""
    # JK*02N.17 needs both its SNPs; losing one collapses the call to JK*02.
    rows = [_db_row(chrom='chr18', genotype='JK*02N.17', genome='45739309_G_A,45739554_G_A')]
    sites = bg_db.site_system_map(rows, 'GRCh38')
    assert [(s.chrom, s.pos, s.system) for s in sites] == [
        ('chr18', 45739309, 'JK'),
        ('chr18', 45739554, 'JK'),
    ]


def test_site_shared_by_two_systems_is_kept_once_per_system():
    """A coordinate used by two systems yields one row per system."""
    # The RHD and RHCE coordinates overlap, so a site is not uniquely keyed by position.
    rows = [
        _db_row(chrom='chr1', genotype='RHD*01', genome='25370466_C_T'),
        _db_row(chrom='chr1', genotype='RHCE*01', genome='25370466_C_T'),
    ]
    assert [s.system for s in bg_db.site_system_map(rows, 'GRCh38')] == ['RHCE', 'RHD']


def test_duplicate_tokens_across_rows_are_deduplicated():
    rows = [_db_row(genome='3774964_A_G'), _db_row(genotype='VEL*01.02', genome='3774964_A_G')]
    assert len(bg_db.site_system_map(rows, 'GRCh38')) == 1


def test_lowercase_chrom_is_normalised_to_gvcf_contig_form():
    """A lower-case db contig is normalised to GVCF form."""
    # One v2.4.1 db row carries `chrx`; it has to match the GVCF's chrX.
    rows = [_db_row(chrom='chrx', genotype='XK*01', genome='37685689_G_A')]
    assert bg_db.site_system_map(rows, 'GRCh38')[0].chrom == 'chrX'


def test_rows_without_coordinates_are_skipped():
    """A db row with no usable coordinate contributes no sites."""
    rows = [_db_row(genome=''), _db_row(genome='.'), _db_row(genome='na')]
    assert bg_db.site_system_map(rows, 'GRCh38') == []


# --- committed resources --------------------------------------------------------------


def _resource(name: str) -> str:
    return files('popgen_rbceq2').joinpath('resources').joinpath(name).read_text()


@pytest.mark.parametrize(
    'name',
    ['bg_regions.GRCh38.bed', 'bg_defining_sites.GRCh38.bed', 'bg_site_systems.GRCh38.tsv'],
)
def test_committed_resources_are_shipped(name):
    assert blood_group_resource(name).endswith(name)


def test_blood_group_resource_raises_for_a_build_we_have_not_generated():
    """Requesting a resource for an unshipped reference build raises, naming how to make it."""
    # GRCh37 is a valid rbceq2 reference_genome, but we ship no resources for it.
    with pytest.raises(FileNotFoundError, match=re.escape('gen_bg_resources.py')):
        blood_group_resource('bg_site_systems.GRCh37.tsv')


def test_committed_bed_and_map_describe_the_same_sites():
    """The committed defining-sites BED and site-system map cover one site set."""
    # A site in one but not the other is either extracted and never assessed, or assessed
    # and never extracted (so flagged NOCOV for no reason). Guards a hand-edit.
    bed_sites = {
        (chrom, int(end))
        for chrom, _start, end in (line.split('\t') for line in _resource('bg_defining_sites.GRCh38.bed').splitlines())
    }
    map_sites = {(s.chrom, s.pos) for s in load_site_systems(_resource('bg_site_systems.GRCh38.tsv'))}
    assert bed_sites == map_sites


def test_committed_map_carries_only_sequence_alleles():
    """Every committed site is defined on sequences, or is a lane site with none."""
    # A structural-variant notation read as an allele renders as `del>21kb` in a flag, and
    # gets a DP/GQ reported for a multi-kb event that one base cannot speak to.
    for site in load_site_systems(_resource('bg_site_systems.GRCh38.tsv')):
        if site.kind == 'ref':
            assert (site.ref, site.alt) == ('.', '.')
            continue
        assert site.kind == 'var'
        where = f'{site.chrom}:{site.pos} ({site.system}) is defined on {site.ref}/{site.alt}, not sequences'
        assert set(site.ref) <= set('ACGTN'), where
        assert set(site.alt) <= set('ACGTN'), where


def test_committed_sites_fall_inside_their_contigs():
    """No committed defining site is past the end of its contig."""
    # An off-contig coordinate is unreachable, so it is NOCOV in every sample forever and
    # reads as a mappability problem rather than the db error it is. v2.4.1 puts ABCC4's
    # chr13 coordinate on ABCC1's chr16 row, which is past the end of chr16 in both builds.
    for site in load_site_systems(_resource('bg_site_systems.GRCh38.tsv')):
        length = GRCH38_CONTIG_LENGTHS[site.chrom]
        assert site.pos <= length, f'{site.chrom}:{site.pos} ({site.system}) is past the end of {site.chrom} ({length})'


def test_committed_bed_sites_fall_inside_the_committed_regions():
    """Every committed defining site falls inside the committed regions BED."""
    # Restriction is unconditional, so a defining site outside the regions BED would be
    # dropped from the conversion and never reach rbceq2 or the QC pass.
    regions: dict[str, list[tuple[int, int]]] = {}
    for line in _resource('bg_regions.GRCh38.bed').splitlines():
        chrom, start, end = line.split('\t')
        regions.setdefault(chrom, []).append((int(start), int(end)))
    for site in load_site_systems(_resource('bg_site_systems.GRCh38.tsv')):
        assert any(start < site.pos <= end for start, end in regions[site.chrom]), (
            f'{site.chrom}:{site.pos} ({site.system}) is outside bg_regions.GRCh38.bed'
        )


# --- extract parsing ------------------------------------------------------------------


def test_parse_extract_reads_blocks_and_variants():
    """Reference blocks and variant records both parse, and blank lines are skipped."""
    text = (
        'chr20\t19999580\tT\t<NON_REF>\t20001310\t0/0\t31\t38\t25\n'
        'chr20\t20003793\tC\tG\t.\t0/1\t25\t38\t.\n'
        '\n'  # bcftools output ends with a newline; a blank line is not a record
    )
    records = parse_extract(text)
    assert records == [
        _record('chr20', 19999580, 'T', '<NON_REF>', end=20001310, gt='0/0', dp=31, gq=38, min_dp=25),
        _record('chr20', 20003793, 'C', 'G', gt='0/1', dp=25, gq=38),
    ]
    assert records[0].is_block
    assert not records[1].is_block
    assert (records[0].span, records[1].span) == (1731, 1)


def test_parse_extract_raises_on_unexpected_column_count():
    with pytest.raises(ValueError, match='columns'):
        parse_extract('chr20\t19999580\tT\n')


def test_load_site_systems_raises_on_wrong_columns():
    with pytest.raises(ValueError, match='must have columns'):
        load_site_systems('chrom\tpos\tsystem\nchr1\t3774964\tVEL\n')


def test_load_site_systems_rejects_a_structural_variant_row():
    """A hand-added SV row in the map raises rather than being assessed."""
    # gen_bg_resources.py never writes one; this catches an edited resource.
    text = 'chrom\tpos\tref\talt\tkind\tsystem\nchr16\t95018451\tdel\t21kb\tsv\tABCC1\n'
    with pytest.raises(ValueError, match='structural-variant site'):
        load_site_systems(text)


# --- coverage resolution --------------------------------------------------------------


def _resolved(records, chrom, pos) -> Coverage:
    """resolve_coverage for a site the test expects to resolve, narrowed to the Coverage."""
    coverage = resolve_coverage(records, chrom, pos)
    assert coverage is not None, f'no record resolved for {chrom}:{pos}'
    return coverage


def test_site_on_its_own_variant_record_uses_that_records_dp_gq():
    """A site with its own record resolves to that record's DP and GQ."""
    records = [_record('chr1', 3774964, 'A', 'G', dp=19, gq=45)]
    coverage = _resolved(records, 'chr1', 3774964)
    assert (coverage.dp, coverage.gq, coverage.source) == (19, 45, 'site')


def test_site_inside_a_reference_block_is_judged_on_the_blocks_dp():
    """A site with no record of its own is judged on the covering block's DP."""
    # The case that makes the flag worth having: with no variant record the call looks like
    # confident reference, and only the covering block says anything about the data under it.
    # DP is the band median, so it estimates this site's depth; MIN_DP is a span floor one
    # shallow base anywhere in the block can set, and is reported rather than tested.
    records = [_record('chr20', 19999580, 'T', '<NON_REF>', end=20001310, dp=31, gq=38, min_dp=8)]
    coverage = _resolved(records, 'chr20', 20000000)
    assert (coverage.dp, coverage.gq, coverage.source) == (31, 38, 'block')
    assert coverage.record.min_dp == 8


def test_a_thin_block_min_dp_alone_does_not_flag_a_site():
    """A block whose MIN_DP is below threshold but whose DP is not is not flagged."""
    # One shallow base in the span cannot condemn every defining site inside the block.
    records = [_record('chr20', 19999580, 'T', '<NON_REF>', end=20001310, dp=31, gq=38, min_dp=8)]
    site = _site(chrom='chr20', pos=20000000, ref='T', alt='C', system='JK')
    coverage = _resolved(records, 'chr20', 20000000)
    assert flag_site(site, coverage, MIN_DEPTH, MIN_GQ) is None


def test_one_base_reference_block_is_assessed_on_its_dp():
    """A block whose END equals its POS covers its site and reads DP, not MIN_DP."""
    # END == POS is not treated as a block, so DP is used. Over a single base the block
    # average (DP) and the block minimum (MIN_DP) are the same number, so a one-base block
    # is assessed identically whichever field is read.
    records = [_record('chr20', 20000000, 'T', '<NON_REF>', end=20000000, dp=31, gq=38, min_dp=31)]
    coverage = _resolved(records, 'chr20', 20000000)
    assert (coverage.dp, coverage.gq, coverage.source) == (31, 38, 'site')


def test_multi_base_reference_block_starting_on_the_site_is_still_a_block():
    """A block that happens to begin on the defining site is a span, not a call at it."""
    # Blocks are banded on GQ, so a block boundary can land on any coordinate, including a
    # defining one. Reading the source off the branch that found the record rather than off
    # the record rendered these as a call at the site, dropping the span and MIN_DP.
    records = [_record('chr20', 20000000, 'T', '<NON_REF>', end=20001310, dp=31, gq=38, min_dp=8)]
    coverage = _resolved(records, 'chr20', 20000000)
    assert (coverage.dp, coverage.gq, coverage.source) == (31, 38, 'block')
    assert coverage.record.min_dp == 8


def test_a_block_starting_on_the_site_reports_its_span_and_min_dp():
    """The flag for such a block names it as a block, as it does for one spanning in."""
    records = [_record('chr1', 3774964, 'A', '<NON_REF>', end=3775064, dp=26, gq=15, min_dp=19)]
    site = _site(chrom='chr1', pos=3774964, ref='A', alt='G', system='VEL')
    coverage = _resolved(records, 'chr1', 3774964)
    assert flag_site(site, coverage, MIN_DEPTH, MIN_GQ) == 'LOWQ:1:3774964(A>G,block=101bp,DP=26,MIN_DP=19,GQ=15)'


def test_a_variant_at_the_site_is_still_a_call_even_beside_a_block_starting_there():
    """The real ALT is preferred first, and a variant record is a call at the site."""
    records = [
        _record('chr1', 3774964, 'A', '<NON_REF>', end=3775064, dp=26, gq=15, min_dp=19),
        _record('chr1', 3774964, 'A', 'G', dp=19, gq=45),
    ]
    coverage = _resolved(records, 'chr1', 3774964)
    assert (coverage.dp, coverage.gq, coverage.source) == (19, 45, 'site')


def test_real_alt_record_wins_over_the_non_ref_twin_at_the_same_position():
    """Where a split leaves two records at one position, the real ALT is chosen."""
    # `bcftools norm -m -any` splits every GVCF variant into the real ALT plus a
    # <NON_REF> twin at the same POS, both carrying identical DP/GQ.
    records = [
        _record('chr20', 20003793, 'C', '<NON_REF>', dp=25, gq=38),
        _record('chr20', 20003793, 'C', 'G', dp=25, gq=38),
    ]
    coverage = _resolved(records, 'chr20', 20003793)
    assert (coverage.dp, coverage.gq) == (25, 38)
    assert coverage.record.alt == 'G'


def test_record_at_the_site_wins_over_a_block_spanning_it():
    records = [
        _record('chr1', 3774900, 'T', '<NON_REF>', end=3775000, dp=40, gq=50, min_dp=9),
        _record('chr1', 3774964, 'A', 'G', dp=19, gq=45),
    ]
    assert _resolved(records, 'chr1', 3774964).dp == 19


def test_innermost_block_wins_when_several_span_the_site():
    records = [
        _record('chr1', 3774900, 'T', '<NON_REF>', end=3775000, dp=40, gq=50, min_dp=40),
        _record('chr1', 3774950, 'T', '<NON_REF>', end=3774980, dp=12, gq=22, min_dp=11),
    ]
    assert _resolved(records, 'chr1', 3774964).dp == 12


def test_deletion_covers_sites_inside_its_ref_allele():
    records = [_record('chr1', 3774960, 'ACGTT', 'A', dp=14, gq=30)]
    coverage = _resolved(records, 'chr1', 3774964)
    assert (coverage.dp, coverage.source) == (14, 'deletion')
    assert resolve_coverage(records, 'chr1', 3774965) is None


def test_deletion_starting_at_the_site_is_not_a_spanning_deletion():
    """A deletion whose POS is the site itself leaves the defining base in place."""
    # The first REF base of a deletion is the retained anchor, so the site is called, not
    # removed. Only a deletion reaching the site from an earlier POS deletes it.
    records = [_record('chr1', 3774964, 'CATGA', 'C', gt='0/1', dp=30, gq=50)]
    assert _resolved(records, 'chr1', 3774964).source == 'site'


def test_records_on_another_contig_do_not_cover_the_site():
    records = [_record('chr2', 3774964, 'A', 'G', dp=30, gq=40)]
    assert resolve_coverage(records, 'chr1', 3774964) is None


def test_site_with_no_covering_record_resolves_to_none():
    assert resolve_coverage([], 'chr1', 3774964) is None


# --- flagging -------------------------------------------------------------------------


def test_site_clearing_both_thresholds_is_not_flagged():
    coverage = _resolved([_record('chr1', 3774964, 'A', 'G', dp=19, gq=45)], 'chr1', 3774964)
    assert flag_site(_site(), coverage, MIN_DEPTH, MIN_GQ) is None


def test_thresholds_are_quiet_at_the_documented_borderline_sites():
    """None of the borderline sites from the design spec are flagged at DP 10 / GQ 20."""
    # Those calls are fixed by no longer dropping the genotypes; flagging them all would
    # make the annotation useless. Values from SPEC-blood-group-call-qc.md section 2.3.
    for dp, gq in [(19, 45), (48, 33), (14, 38), (28, 28), (24, 43), (18, 39), (21, 26), (52, 21)]:
        records = [_record('chr1', 3774964, 'A', 'G', dp=dp, gq=gq)]
        coverage = _resolved(records, 'chr1', 3774964)
        assert flag_site(_site(), coverage, MIN_DEPTH, MIN_GQ) is None, f'DP={dp} GQ={gq}'


@pytest.mark.parametrize(
    ('dp', 'gq', 'expected'),
    [
        (8, 45, f'{LOWQ}:1:3774964(A>G,DP=8,GQ=45)'),
        (30, 15, f'{LOWQ}:1:3774964(A>G,DP=30,GQ=15)'),
        # A missing DP or GQ is not evidence of a good call, so it flags like a low value.
        (None, 45, f'{LOWQ}:1:3774964(A>G,DP=.,GQ=45)'),
        (30, None, f'{LOWQ}:1:3774964(A>G,DP=30,GQ=.)'),
    ],
)
def test_low_quality_site_is_flagged_with_its_coordinate_and_metrics(dp, gq, expected):
    coverage = _resolved([_record('chr1', 3774964, 'A', 'G', dp=dp, gq=gq)], 'chr1', 3774964)
    assert flag_site(_site(), coverage, MIN_DEPTH, MIN_GQ) == expected


def test_uncovered_site_is_flagged_nocov():
    assert flag_site(_site(), None, MIN_DEPTH, MIN_GQ) == f'{NOCOV}:1:3774964(A>G)'


@pytest.mark.parametrize('span', [(20001310, '1.7kb'), (19999600, '21bp')])
def test_a_flag_from_a_spanning_block_names_the_span_and_both_depths(span):
    """A flag whose depth came from a block gives its span, its DP and its MIN_DP."""
    # `LOWQ:...(T>C,DP=31,GQ=15)` alone reads as a low-quality T>C call at DP 31. The 31 is
    # the block's median depth, so the flag also gives the span and the span's floor: a reader
    # can see the depth is not this site's own, and how low any base in the block went.
    end, rendered = span
    site = _site(chrom='chr20', pos=19999590, ref='T', alt='C', system='JK')
    records = [_record('chr20', 19999580, 'T', '<NON_REF>', end=end, dp=31, gq=15, min_dp=19)]
    coverage = _resolved(records, 'chr20', 19999590)
    expected = f'{LOWQ}:20:19999590(T>C,block={rendered},DP=31,MIN_DP=19,GQ=15)'
    assert flag_site(site, coverage, MIN_DEPTH, MIN_GQ) == expected


def test_a_block_carrying_no_min_dp_renders_it_absent():
    """A block with no MIN_DP field renders it as `.`, not as a fabricated number."""
    records = [_record('chr20', 19999580, 'T', '<NON_REF>', end=19999600, dp=8, gq=45)]
    site = _site(chrom='chr20', pos=19999590, ref='T', alt='C', system='JK')
    coverage = _resolved(records, 'chr20', 19999590)
    assert flag_site(site, coverage, MIN_DEPTH, MIN_GQ) == f'{LOWQ}:20:19999590(T>C,block=21bp,DP=8,MIN_DP=.,GQ=45)'


@pytest.mark.parametrize('gt', ['0/1', '1/1'])
def test_a_site_a_deletion_swallows_is_flagged_deleted(gt):
    """A confidently-called deletion over the defining base is DEL, not PASS."""
    # The sample carries a 5bp deletion at 3774961 spanning 3774961-3774965, so the A that
    # VEL's allele is defined on is not in this sample. rbceq2 sees no variant at 3774964
    # and calls VEL reference; the deletion's own DP and GQ are excellent, so on quality
    # alone this site reads as a well-supported call.
    records = [_record('chr1', 3774961, 'CATGA', 'C', gt=gt, dp=30, gq=50)]
    coverage = _resolved(records, 'chr1', 3774964)
    expected = f'{DELETED}:1:3774964(A>G,del=CATGA>C,GT={gt},DP=30,GQ=50)'
    assert flag_site(_site(), coverage, MIN_DEPTH, MIN_GQ) == expected


def test_deleted_outranks_low_quality():
    """A deletion over the defining base is reported as DEL even when its own DP is thin."""
    # The base being absent is the finding; its depth cannot make the call better supported.
    records = [_record('chr1', 3774961, 'CATGA', 'C', gt='1/1', dp=4, gq=8)]
    coverage = _resolved(records, 'chr1', 3774964)
    assert flag_site(_site(), coverage, MIN_DEPTH, MIN_GQ) == f'{DELETED}:1:3774964(A>G,del=CATGA>C,GT=1/1,DP=4,GQ=8)'


def test_lane_site_renders_as_ref_having_no_allele_of_its_own():
    site = _site(pos=207331122, ref='.', alt='.', kind='ref', system='CROM')
    coverage = _resolved([_record('chr1', 207331122, 'G', '<NON_REF>', dp=5, gq=12)], 'chr1', 207331122)
    assert flag_site(site, coverage, MIN_DEPTH, MIN_GQ) == f'{LOWQ}:1:207331122(ref,DP=5,GQ=12)'


def test_long_alleles_are_rendered_as_lengths():
    """An allele too long to show is rendered as its length."""
    # The FY and XK large indels carry ~200-base REFs that would otherwise swamp the cell.
    site = _site(chrom='chr1', pos=159205730, ref='T' * 198, alt='T', system='FY')
    coverage = _resolved([_record('chr1', 159205730, 'T' * 198, 'T', dp=4, gq=30)], 'chr1', 159205730)
    assert flag_site(site, coverage, MIN_DEPTH, MIN_GQ) == f'{LOWQ}:1:159205730(198bp>T,DP=4,GQ=30)'


def test_a_long_deletion_in_a_del_flag_is_rendered_as_a_length():
    """The deletion named in a DEL flag is length-capped like any other allele."""
    records = [_record('chr1', 3774900, 'A' + 'C' * 200, 'A', gt='1/1', dp=30, gq=50)]
    coverage = _resolved(records, 'chr1', 3774964)
    assert flag_site(_site(), coverage, MIN_DEPTH, MIN_GQ) == f'{DELETED}:1:3774964(A>G,del=201bp>A,GT=1/1,DP=30,GQ=50)'


# --- aggregation ----------------------------------------------------------------------


def test_clean_system_is_pass_and_flagged_sites_are_joined_in_coordinate_order():
    sites = [
        _site(chrom='chr18', pos=45739554, ref='G', alt='A', system='JK'),
        _site(chrom='chr18', pos=45739309, ref='G', alt='A', system='JK'),
        _site(system='VEL'),
    ]
    records = [
        _record('chr18', 45739309, 'G', 'A', dp=4, gq=37),
        _record('chr18', 45739554, 'G', 'A', dp=30, gq=8),
        _record('chr1', 3774964, 'A', 'G', dp=40, gq=45),
    ]
    flags, uncovered = flags_by_system(sites, records, MIN_DEPTH, MIN_GQ)
    assert flags == {
        'JK': f'{LOWQ}:18:45739309(G>A,DP=4,GQ=37);{LOWQ}:18:45739554(G>A,DP=30,GQ=8)',
        'VEL': PASS,
    }
    assert uncovered == []


def test_uncovered_sites_are_reported_for_logging():
    """Sites with no covering record are returned alongside the flags."""
    sites = [_site(), _site(chrom='chr18', pos=45739309, system='JK')]
    flags, uncovered = flags_by_system(sites, [], MIN_DEPTH, MIN_GQ)
    assert [(s.chrom, s.pos) for s in uncovered] == [('chr1', 3774964), ('chr18', 45739309)]
    assert flags['VEL'].startswith(NOCOV)


# --- TSV rendering --------------------------------------------------------------------


def test_qc_tsv_mirrors_the_geno_tsv_header_and_sample_row():
    """The UUID header cell and the sample id are copied from the geno TSV verbatim."""
    # That is what lets the QC TSV join to the calls 1:1 by column and concatenate across
    # a cohort the same way.
    geno = 'UUID: abc123\tJK\tVEL\nSAMPLE1\tJK*01/JK*02\tVEL*01/VEL*01\n'
    qc = build_qc_tsv(geno, {'JK': PASS, 'VEL': f'{LOWQ}:1:3774964(A>G,DP=8,GQ=45)'})
    assert qc == ('UUID: abc123\tJK\tVEL\nSAMPLE1\tPASS\tLOWQ:1:3774964(A>G,DP=8,GQ=45)\n')


def test_system_called_by_rbceq2_but_absent_from_the_map_is_not_assessed():
    """A geno-TSV column with no defining site in the map is NOT_ASSESSED."""
    # rbceq2 emits 48 systems by default while the map covers all 87 in the db, so the
    # reverse can happen after a db bump. Such a system is neither PASS nor LOWQ.
    geno = 'UUID: abc123\tJK\tSOMETHINGNEW\nSAMPLE1\tJK*01/JK*02\tX/Y\n'
    assert build_qc_tsv(geno, {'JK': PASS}).splitlines()[1] == f'SAMPLE1\t{PASS}\t{NOT_ASSESSED}'


def test_build_qc_tsv_raises_without_a_sample_row():
    with pytest.raises(ValueError, match='single sample row'):
        build_qc_tsv('UUID: abc123\tJK\tVEL\n', {})


def test_build_qc_tsv_raises_on_ragged_geno_tsv():
    with pytest.raises(ValueError, match='header has'):
        build_qc_tsv('UUID: abc123\tJK\tVEL\nSAMPLE1\tJK*01/JK*02\n', {})
