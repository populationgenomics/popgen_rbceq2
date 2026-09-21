"""What the merge does with a post-hoc record that reaches past the hole it was kept for.

A record is selected for covering a fillable hole and is selected whole, so one anchored in a
hole can extend over a neighbouring defining site it may not fill: one the primary caller did
call, or a hole inside the capture design. Blood-group defining sites are dense enough that
this is the ordinary case near a capture edge. What the merge does next differs by record type,
and both halves rest on bcftools semantics rather than on any Python here:

- a post-hoc **variant** reaching such a site is dropped, because keeping it hands rbceq2 two
  callers' alleles at one base, or lets the second caller decide a site the design targeted;
- a post-hoc **reference block** reaching one is kept, because it asserts nothing rbceq2 sees
  and dropping it would throw away the hole it was kept for. The QC then reads the fillable
  sites and disregards the block anywhere else, which is what keeps an in-design hole NOCOV.

So these tests run the real merge shell under real bcftools, then read the result with the real
extract format and the real QC functions, handed the same off-design BED the merge read. A
Python stand-in would test a restatement of `annotate -m`'s overlap rule rather than the rule,
and that rule is the whole mechanism: `-m` marks on a record's span, so it catches a deletion
anchored on a hole that reaches a called site one base away, and it runs before `norm -m -any`
splits a variant from its `<NON_REF>` allele, so a variant is recognised by its alleles and no
split-off twin is around to pass as a block. The mark is the unfillable defining sites rather
than the DRAGEN records' spans, which is the last test of the called-site group here.

Skipped where bcftools is not installed. Local bcftools is expected to be the pinned image's
1.24, so the semantics checked here are the ones the job will meet.
"""

import shutil
import subprocess
from pathlib import Path

import pytest

from popgen_rbceq2.constants import POSTHOC_CALLER
from popgen_rbceq2.jobs.rbceq2_call_qc_job import flags_by_system, load_fillable_sites, parse_extract
from popgen_rbceq2.scripts.bg_db import DefiningSite
from popgen_rbceq2.stages.blood_group_genotyping.filter_and_convert import (
    _EXTRACT_FORMAT,
    _POSTHOC_HEADER_LINE,
    _merge_posthoc_commands,
    _primary_records_guard,
    _sample_check_commands,
)

pytestmark = [
    pytest.mark.fast,
    pytest.mark.skipif(
        not all(shutil.which(tool) for tool in ('bcftools', 'bgzip', 'tabix')),
        reason='bcftools, bgzip and tabix are not all on PATH',
    ),
]

# A defining site the capture design targeted, and one just outside it. The design below puts
# the boundary between them, which is what makes the second a hole and the first not.
IN_DESIGN = 2005
OFF_DESIGN = 2000
SITES_BED = f'chr1\t{OFF_DESIGN - 1}\t{OFF_DESIGN}\nchr1\t{IN_DESIGN - 1}\t{IN_DESIGN}\n'
# What the committed off-design resource hands this job: the defining sites the design never
# targeted, which here is the second one and not the first. The subtraction itself is tested in
# test_off_design_subtraction; this file starts from its result.
OFF_DESIGN_BED = f'chr1\t{OFF_DESIGN - 1}\t{OFF_DESIGN}\n'
# A third site the design did target and DRAGEN still has no record at: the in-design hole the
# fill must leave alone. Tests that need it pass this in place of SITES_BED; the off-design set
# is unchanged, because being a hole does not make a site fillable.
IN_DESIGN_HOLE = 2010
SITES_WITH_IN_DESIGN_HOLE = SITES_BED + f'chr1\t{IN_DESIGN_HOLE - 1}\t{IN_DESIGN_HOLE}\n'

# DRAGEN calls the in-design site and has no record at all at the off-design one, which is the
# capture edge this feature exists for.
DRAGEN_RECORDS = f'chr1\t{IN_DESIGN}\t.\tA\tG,<NON_REF>\t200\tPASS\tSPARE=1\tGT:DP:GQ\t0/1:50:50\n'

# A gVCF header declares far more tags than any one record carries, which matters: the merge's
# `annotate -x '^INFO/END,...'` reads the header, and errors outright if the header has nothing
# outside its keep list. SPARE stands in for the BaseQRankSum, MLEAC and friends a real
# HaplotypeCaller header declares.
HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=250000000>
##contig=<ID=chrX,length=156040895>
##ALT=<ID=NON_REF,Description="Represents any possible alternative allele">
##INFO=<ID=END,Number=1,Type=Integer,Description="Block end position">
##INFO=<ID=SPARE,Number=1,Type=Integer,Description="A tag nothing downstream reads">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum depth over the block">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
"""


def _bgzip(
    tmp_path: Path, name: str, body: str, *, index: bool, extra_header: str = '', sample: str = 'SAMPLE1'
) -> Path:
    """Write a VCF body under HEADER and bgzip it, indexing only when asked.

    `extra_header` is one more `##` line, placed before the column header. `sample` renames
    the one sample column.
    """
    plain = tmp_path / name
    header = HEADER if not extra_header else HEADER.replace('#CHROM', f'{extra_header}\n#CHROM')
    header = header.replace('\tSAMPLE1\n', f'\t{sample}\n')
    plain.write_text(header + body)
    packed = tmp_path / f'{name}.gz'
    with packed.open('wb') as out:
        subprocess.run(['bgzip', '-c', str(plain)], stdout=out, check=True)  # noqa: S603, S607
    if index:
        subprocess.run(['bcftools', 'index', '-t', str(packed)], check=True)  # noqa: S603, S607
    return packed


def _merge(
    tmp_path: Path,
    posthoc_records: str,
    *,
    dragen_records: str = DRAGEN_RECORDS,
    sites_bed: str = SITES_BED,
    off_design_bed: str = OFF_DESIGN_BED,
    posthoc_sample: str = 'SAMPLE1',
    check: bool = True,
) -> subprocess.CompletedProcess[str]:
    """Run the real merge shell over one post-hoc record set, in tmp_path.

    `posthoc_sample` names the post-hoc file's sample; the DRAGEN file is always SAMPLE1.
    `check=False` for a merge expected to fail, so the test can read the exit code and stderr.
    """
    posthoc = _bgzip(tmp_path, 'posthoc.g.vcf', posthoc_records, index=True, sample=posthoc_sample)
    sites = tmp_path / 'sites.bed'
    sites.write_text(sites_bed)
    off_design = tmp_path / 'off_design_defining_sites.bed'
    off_design.write_text(off_design_bed)
    # What the caller leaves behind for the fragment: the POSTHOC declaration both sides of
    # the concat have to carry, and the bg-regions intermediate, which the stage has already
    # stamped with that declaration so the extract can read the tag whether or not anything
    # was merged. Not indexed, because the fragment's own first line indexes it.
    (tmp_path / 'posthoc_hdr.txt').write_text(f'{_POSTHOC_HEADER_LINE}\n')
    dragen = _bgzip(tmp_path, 'dragen.vcf', dragen_records, index=False, extra_header=_POSTHOC_HEADER_LINE)

    # The sample check runs first in the stage's job and reads the raw gVCF; the intermediate
    # carries the same header, so it stands in for the raw file here.
    script = (
        'set -euo pipefail\n'
        + _sample_check_commands(str(posthoc), str(dragen))
        + _merge_posthoc_commands(
            str(posthoc),
            str(sites),
            str(off_design),
            'exome_probesets_hg38/test_design_bed',
            cpu=1,
            genome='GRCh38',
        )
    )
    return subprocess.run(  # noqa: S603
        ['bash', '-c', script],  # noqa: S607
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=check,
    )


def _guard_then_merge(
    tmp_path: Path,
    posthoc_records: str,
    *,
    dragen_records: str,
) -> subprocess.CompletedProcess[str]:
    """Run the wrong-build guard and then the merge, the order the stage runs them in."""
    posthoc = _bgzip(tmp_path, 'posthoc.g.vcf', posthoc_records, index=True)
    sites = tmp_path / 'sites.bed'
    sites.write_text(SITES_BED)
    off_design = tmp_path / 'off_design_defining_sites.bed'
    off_design.write_text(OFF_DESIGN_BED)
    (tmp_path / 'posthoc_hdr.txt').write_text(f'{_POSTHOC_HEADER_LINE}\n')
    _bgzip(tmp_path, 'dragen.vcf', dragen_records, index=False, extra_header=_POSTHOC_HEADER_LINE)

    script = (
        'set -euo pipefail\n'
        + _primary_records_guard(str(sites), 'GRCh38')
        + _merge_posthoc_commands(
            str(posthoc),
            str(sites),
            str(off_design),
            'exome_probesets_hg38/test_design_bed',
            cpu=1,
            genome='GRCh38',
        )
    )
    return subprocess.run(  # noqa: S603
        ['bash', '-c', script],  # noqa: S607
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )


def _flags(tmp_path: Path, *, sites_bed: str = SITES_BED, off_design_bed: str = OFF_DESIGN_BED) -> dict[int, str]:
    """Extract merged.vcf.gz the way the stage does and flag every defining site.

    Runs the QC's own aggregation, handed the off-design BED as the fillable set the way the
    QC stage hands it to the job, so the result is what the QC TSV would say at each site.
    """
    out = subprocess.run(  # noqa: S603
        [  # noqa: S607
            'bcftools',
            'query',
            '-T',
            str(tmp_path / 'sites.bed'),
            '--targets-overlap',
            '2',
            '-f',
            _EXTRACT_FORMAT.replace(r'\t', '\t').replace(r'\n', '\n'),
            str(tmp_path / 'merged.vcf.gz'),
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    records = parse_extract(out.stdout)
    # The ref/alt here only shape how a flag renders the site, not which record resolves it.
    # One system per site, so a system's cell is that site's flag.
    sites = [
        DefiningSite(
            chrom='chr1', pos=int(line.split('\t')[2]), ref='C', alt='T', kind='var', system=f'S{line.split()[2]}'
        )
        for line in sites_bed.splitlines()
    ]
    cells, _ = flags_by_system(sites, records, 10, 20, load_fillable_sites(off_design_bed))
    return {site.pos: cells[site.system] for site in sites}


def test_a_posthoc_deletion_reaching_a_called_base_is_dropped_and_the_site_stays_a_hole(tmp_path):
    # The finding this file exists for. Kept, the deletion would remove the base DRAGEN called
    # a SNP on, in the file rbceq2 reads, and the QC would still say PASS there because
    # resolve_coverage prefers the primary record. Dropped, the hole it filled goes back to
    # NOCOV, which is the honest answer: DRAGEN spoke at the base this record reaches.
    deletion = f'chr1\t{OFF_DESIGN}\t.\tCATGAAAAA\tC,<NON_REF>\t60\t.\tSPARE=1\tGT:DP:GQ\t0/1:44:80\n'

    _merge(tmp_path, deletion)

    assert _flags(tmp_path) == {OFF_DESIGN: 'NOCOV:1:2000(C>T)', IN_DESIGN: 'PASS'}


def test_a_posthoc_deletion_carrying_an_end_tag_is_still_dropped(tmp_path):
    # Harper's second-review repro. HaplotypeCaller never writes INFO/END on a variant, so a
    # drop keyed on "no END" worked on its output alone; a variant that carries END must go
    # too, because what makes it a variant is its alleles, not the absence of a tag one caller
    # happens not to write. Its <NON_REF> twin would inherit the END and pass as a block, which
    # is why the drop runs before norm splits one off.
    deletion = f'chr1\t{OFF_DESIGN}\t.\tACGTAC\tA,<NON_REF>\t60\t.\tEND={IN_DESIGN};SPARE=1\tGT:DP:GQ\t0/1:44:80\n'

    stderr = _merge(tmp_path, deletion).stderr

    assert '1 variant(s) dropped for' in stderr
    assert _flags(tmp_path) == {OFF_DESIGN: 'NOCOV:1:2000(C>T)', IN_DESIGN: 'PASS'}


def test_a_posthoc_reference_block_reaching_a_called_base_is_kept_and_fills_its_hole(tmp_path):
    # The case Harper's one-line `-T ^covered.bed` would have broken. A block straddling the
    # capture edge is how many holes get filled at all, and it asserts nothing rbceq2 sees:
    # the conversion drops every <NON_REF>-only record before rbceq2 reads the file.
    block = f'chr1\t1990\t.\tG\t<NON_REF>\t.\t.\tEND={IN_DESIGN + 5};SPARE=1\tGT:DP:GQ:MIN_DP\t0/0:40:60:35\n'

    _merge(tmp_path, block)

    flags = _flags(tmp_path)
    assert flags[OFF_DESIGN].startswith('POSTHOC:1:2000(')
    assert 'src=gatk-hc-' in flags[OFF_DESIGN]
    # The primary record still wins at the site DRAGEN called, which is the existing rule the
    # kept block relies on.
    assert flags[IN_DESIGN] == 'PASS'


def test_a_multiallelic_trespasser_is_counted_once(tmp_path):
    # The log promises variants dropped, not alleles. HaplotypeCaller does emit multiallelic
    # records, and after `norm -m -any` one would count as two.
    deletion = f'chr1\t{OFF_DESIGN}\t.\tCATGAAAAA\tC,CA,<NON_REF>\t60\t.\tSPARE=1\tGT:DP:GQ\t1/2:44:80\n'

    stderr = _merge(tmp_path, deletion).stderr

    assert '1 variant(s) dropped for' in stderr
    assert _flags(tmp_path) == {OFF_DESIGN: 'NOCOV:1:2000(C>T)', IN_DESIGN: 'PASS'}


def test_a_posthoc_variant_confined_to_its_hole_is_kept(tmp_path):
    # The exclusion has to be narrow. A variant that reaches no called base is exactly what
    # this feature is for, and must survive untouched.
    snp = f'chr1\t{OFF_DESIGN}\t.\tC\tT,<NON_REF>\t60\t.\tSPARE=1\tGT:DP:GQ\t0/1:44:80\n'

    _merge(tmp_path, snp)

    flags = _flags(tmp_path)
    assert flags[OFF_DESIGN].startswith('POSTHOC:1:2000(C>T,src=gatk-hc-')
    assert flags[IN_DESIGN] == 'PASS'


def test_the_merge_reports_what_it_kept_and_what_it_dropped(tmp_path):
    # The old log line reported the hole count and nothing else, so a run that merged nothing
    # at all — every hole zero-depth, or every record dropped for trespass — was silent.
    deletion = f'chr1\t{OFF_DESIGN}\t.\tCATGAAAAA\tC,<NON_REF>\t60\t.\tSPARE=1\tGT:DP:GQ\t0/1:44:80\n'

    stderr = _merge(tmp_path, deletion).stderr

    assert '0 record(s) kept over 1 hole(s)' in stderr
    assert '1 variant(s) dropped for' in stderr


def test_a_posthoc_variant_overlapping_a_dragen_block_but_no_called_site_is_kept(tmp_path):
    # The mark is the DRAGEN-covered defining *sites*, not the DRAGEN records' spans. A long
    # reference block reaches far past the site it covers, and a post-hoc variant that only
    # clips its tail contradicts nothing rbceq2 reads: rbceq2 looks at defining sites, and
    # this record touches none that DRAGEN called. Marking on spans would drop it and lose the
    # hole it was kept for.
    dragen = (
        f'chr1\t{IN_DESIGN}\t.\tA\tG,<NON_REF>\t200\tPASS\tSPARE=1\tGT:DP:GQ\t0/1:50:50\n'
        f'chr1\t{IN_DESIGN + 1}\t.\tT\t<NON_REF>\t.\t.\tEND={OFF_DESIGN + 900};SPARE=1'
        '\tGT:DP:GQ:MIN_DP\t0/0:45:70:40\n'
    )
    # A hole beyond that block, and a deletion anchored in it whose REF clips the block's tail.
    far_site = OFF_DESIGN + 902
    sites = SITES_BED + f'chr1\t{far_site - 1}\t{far_site}\n'
    off_design = OFF_DESIGN_BED + f'chr1\t{far_site - 1}\t{far_site}\n'
    deletion = f'chr1\t{far_site - 3}\t.\tCATGA\tC,<NON_REF>\t60\t.\tSPARE=1\tGT:DP:GQ\t0/1:44:80\n'

    _merge(tmp_path, deletion, dragen_records=dragen, sites_bed=sites, off_design_bed=off_design)

    # Kept, so the site is no longer a hole. It reads DEL+POSTHOC rather than plain POSTHOC
    # because the hole sits inside the deletion rather than on its anchor, which it must: the
    # anchor has to be inside DRAGEN's block for this record to clip the block at all.
    flags = _flags(tmp_path, sites_bed=sites, off_design_bed=off_design)
    assert flags[far_site] == f'DEL+POSTHOC:1:{far_site}(C>T,src=gatk-hc-4.6.2.0,del=CATGA>C,GT=0/1,DP=44,GQ=80)'


# --- a hole inside the design, which the fill must leave alone ---


def test_a_posthoc_block_reaching_an_in_design_hole_fills_only_the_off_design_one(tmp_path):
    # The headline of the design bound: a site the design targeted and DRAGEN still said nothing
    # about stays NOCOV. The block is kept whole for the off-design hole, so it is the only
    # record at the in-design hole beside it, and the merge cannot clip it. The QC is what
    # keeps the promise, by disregarding a post-hoc record at a site the merge could not fill.
    block = f'chr1\t1990\t.\tG\t<NON_REF>\t.\t.\tEND={IN_DESIGN_HOLE + 5};SPARE=1\tGT:DP:GQ:MIN_DP\t0/0:40:60:35\n'

    stderr = _merge(tmp_path, block, sites_bed=SITES_WITH_IN_DESIGN_HOLE).stderr

    assert 'capture design and filled, 1 inside it and left for the QC to flag NOCOV' in stderr
    flags = _flags(tmp_path, sites_bed=SITES_WITH_IN_DESIGN_HOLE)
    assert flags[OFF_DESIGN].startswith('POSTHOC:1:2000(')
    assert flags[IN_DESIGN] == 'PASS'
    assert flags[IN_DESIGN_HOLE] == f'NOCOV:1:{IN_DESIGN_HOLE}(C>T)'


def test_a_posthoc_deletion_reaching_an_in_design_hole_is_dropped(tmp_path):
    # Kept, this deletion would reach rbceq2 and decide the genotype at a site the design
    # targeted, from a caller the design bound exists to keep out of it. It reaches no site
    # DRAGEN called, so the mark has to be the unfillable sites and not the called ones.
    # DRAGEN's own record sits far away so that neither hole has a primary record.
    dragen = 'chr1\t2500\t.\tA\tG,<NON_REF>\t200\tPASS\tSPARE=1\tGT:DP:GQ\t0/1:50:50\n'
    sites = SITES_WITH_IN_DESIGN_HOLE.replace(f'chr1\t{IN_DESIGN - 1}\t{IN_DESIGN}\n', 'chr1\t2499\t2500\n')
    # REF runs from the off-design hole over the in-design one.
    deletion = f'chr1\t{OFF_DESIGN}\t.\tCATGAAAAAAA\tC,<NON_REF>\t60\t.\tSPARE=1\tGT:DP:GQ\t0/1:44:80\n'

    stderr = _merge(tmp_path, deletion, dragen_records=dragen, sites_bed=sites).stderr

    assert '1 variant(s) dropped for' in stderr
    assert _flags(tmp_path, sites_bed=sites) == {
        OFF_DESIGN: 'NOCOV:1:2000(C>T)',
        2500: 'PASS',
        IN_DESIGN_HOLE: f'NOCOV:1:{IN_DESIGN_HOLE}(C>T)',
    }


# --- FILTER on the supplement: what rbceq2 will and will not use ---


def _filters(tmp_path: Path) -> dict[tuple[int, str], str]:
    """FILTER of every record in merged.vcf.gz, keyed by (POS, ALT)."""
    out = subprocess.run(  # noqa: S603
        ['bcftools', 'query', '-f', '%POS\\t%ALT\\t%FILTER\\n', str(tmp_path / 'merged.vcf.gz')],  # noqa: S607
        capture_output=True,
        text=True,
        check=True,
    )
    rows = (line.split('\t') for line in out.stdout.splitlines())
    return {(int(pos), alt): flt for pos, alt, flt in rows}


def test_a_recovered_variant_reaches_rbceq2_as_pass(tmp_path):
    # HaplotypeCaller leaves FILTER `.` and rbceq2 keeps an allele only when its defining
    # variant is literally PASS, so without this every recovered alternate allele was silently
    # discarded and the site typed as reference by absence, under a flag saying otherwise.
    snp = f'chr1\t{OFF_DESIGN}\t.\tC\tT,<NON_REF>\t60\t.\tSPARE=1\tGT:DP:GQ\t0/1:44:80\n'

    _merge(tmp_path, snp)

    filters = _filters(tmp_path)
    assert filters[(OFF_DESIGN, 'T')] == 'PASS'
    assert filters[(OFF_DESIGN, '<NON_REF>')] == 'PASS'
    assert filters[(IN_DESIGN, 'G,<NON_REF>')] == 'PASS'


def test_a_single_read_recovered_variant_is_marked_low_depth_as_dragen_would(tmp_path):
    # DRAGEN's LowDepth rule, DP<=1, is the one filter on these cohorts that means the same on
    # both callers, so it is the one copied. rbceq2 then excludes the allele exactly as it
    # excludes a DRAGEN LowDepth call, and the QC still reports the site's DP and GQ.
    snp = f'chr1\t{OFF_DESIGN}\t.\tC\tT,<NON_REF>\t12\t.\tSPARE=1\tGT:DP:GQ\t0/1:1:3\n'

    _merge(tmp_path, snp)

    assert _filters(tmp_path)[(OFF_DESIGN, 'T')] == 'LowDepth'
    assert _flags(tmp_path)[OFF_DESIGN] == 'LOWQ+POSTHOC:1:2000(C>T,src=gatk-hc-4.6.2.0,DP=1,GQ=3)'


# --- the capture-design gate: what proves the wrong design, and what does not ---


def test_a_dragen_deletion_running_off_the_capture_edge_is_a_carrier_not_a_wrong_design(tmp_path):
    # DRAGEN anchors a variant inside the target but its REF can run past the edge, so a
    # deletion reaching an off-design site is biology. Gating on every record's span failed
    # this sample's whole conversion, blaming a config key that was right. The site is covered,
    # so not filled, and reads DEL from DRAGEN's own record with no post-hoc provenance.
    # As the stage's `norm -m -any` leaves the deletion: the real ALT and its <NON_REF> twin.
    dragen = (
        f'chr1\t{OFF_DESIGN - 3}\t.\tCATGA\tC\t200\tPASS\tSPARE=1\tGT:DP:GQ\t0/1:50:50\n'
        f'chr1\t{OFF_DESIGN - 3}\t.\tCATGA\t<NON_REF>\t200\tPASS\tSPARE=1\tGT:DP:GQ\t0/0:50:50\n'
        f'chr1\t{IN_DESIGN}\t.\tA\tG,<NON_REF>\t200\tPASS\tSPARE=1\tGT:DP:GQ\t0/1:50:50\n'
    )
    snp = f'chr1\t{OFF_DESIGN}\t.\tC\tT,<NON_REF>\t60\t.\tSPARE=1\tGT:DP:GQ\t0/1:44:80\n'

    stderr = _merge(tmp_path, snp, dragen_records=dragen).stderr

    assert '0 defining site(s) with no DRAGEN record' in stderr
    assert _flags(tmp_path) == {OFF_DESIGN: 'DEL:1:2000(C>T,del=CATGA>C,GT=0/1,DP=50,GQ=50)', IN_DESIGN: 'PASS'}


def test_a_dragen_reference_block_over_an_off_design_site_fails_the_job_naming_the_key(tmp_path):
    # A reference block asserts hom-ref over bases DRAGEN evaluated, which stop at the target
    # edge. One reaching a site the configured design says is untargeted means the design is not
    # the gVCF's, and the job says which key to fix rather than quietly filling too little.
    dragen = (
        f'chr1\t{OFF_DESIGN - 10}\t.\tG\t<NON_REF>\t.\t.\tEND={IN_DESIGN + 10};SPARE=1\tGT:DP:GQ:MIN_DP\t0/0:40:60:35\n'
    )
    snp = f'chr1\t{OFF_DESIGN}\t.\tC\tT,<NON_REF>\t60\t.\tSPARE=1\tGT:DP:GQ\t0/1:44:80\n'

    result = _merge(tmp_path, snp, dragen_records=dragen, check=False)

    assert result.returncode == 1
    assert 'DRAGEN reference blocks cover 1 defining site(s) outside the capture design' in result.stderr
    assert 'exome_design_bed = exome_probesets_hg38/test_design_bed is not the BED' in result.stderr
    assert f'chr1\t{OFF_DESIGN - 1}\t{OFF_DESIGN}' in result.stderr
    assert not (tmp_path / 'merged.vcf.gz').exists()


def test_a_posthoc_file_naming_another_sample_fails_even_when_there_is_nothing_to_fill(tmp_path):
    # The check is about the two inputs, not about the data in them, so it must not depend on
    # the data: here DRAGEN covers every defining site, no hole is fillable, and the mismatch
    # still fails the job before anything is computed.
    dragen_covering_both = 'chr1\t1990\t.\tA\t<NON_REF>\t.\t.\tEND=2010\tGT:DP:GQ:MIN_DP\t0/0:50:50:50\n'
    snp = f'chr1\t{OFF_DESIGN}\t.\tC\tT,<NON_REF>\t80\t.\t.\tGT:DP:GQ\t0/1:44:80\n'

    result = _merge(tmp_path, snp, dragen_records=dragen_covering_both, posthoc_sample='SAMPLE2', check=False)

    assert result.returncode == 1
    assert 'name different samples' in result.stderr
    assert not (tmp_path / 'covered.bed').exists()


def test_a_posthoc_file_naming_another_sample_fails_the_job_rather_than_relabelling(tmp_path):
    # The CRAM's read group names the post-hoc sample and DRAGEN named the gVCF from the same
    # run, so a disagreement means the two inputs do not describe one individual. Relabelling
    # would satisfy `concat` and splice another person's genotypes in at exactly the sites
    # nothing else covers, looking like an ordinary recovery. The check runs after the supplement
    # is built, so a mismatch has to fail even when there is something to merge.
    snp = f'chr1\t{OFF_DESIGN}\t.\tC\tT,<NON_REF>\t60\t.\tSPARE=1\tGT:DP:GQ\t0/1:44:80\n'

    result = _merge(tmp_path, snp, posthoc_sample='SAMPLE2', check=False)

    assert result.returncode == 1
    assert 'the post-hoc calls and the gVCF name different samples' in result.stderr
    assert 'CRAM/post-hoc: SAMPLE2' in result.stderr
    assert 'primary gVCF:  SAMPLE1' in result.stderr
    assert not (tmp_path / 'merged.vcf.gz').exists()


# --- the wrong-build guard, and why it runs before the merge ---


def test_a_gvcf_with_no_dragen_record_at_any_defining_site_fails_before_the_merge(tmp_path):
    # A gVCF called against another build, or naming its contigs without `chr`, matches no
    # blood-group region and leaves an empty DRAGEN set with exit 0. The post-hoc caller takes
    # its contigs from the reference fasta, so its records are well-formed regardless, and
    # would fill every off-design hole and populate the extract. Asked after the merge, "does
    # any record overlap a defining site" is therefore answered yes by the supplement alone,
    # and the sample is typed from the second caller with everything else NOCOV. The guard
    # reads the DRAGEN intermediate, before the merge, so this input fails as it did before
    # post-hoc calling existed.
    posthoc = f'chr1\t{OFF_DESIGN}\t.\tC\tT,<NON_REF>\t80\t.\t.\tGT:DP:GQ\t0/1:44:80\n'

    result = _guard_then_merge(tmp_path, posthoc, dragen_records='')

    assert result.returncode == 1
    assert 'no DRAGEN gVCF record overlaps any blood-group defining site' in result.stderr
    assert not (tmp_path / 'merged.vcf.gz').exists()


def test_the_guard_passes_a_gvcf_dragen_called_a_defining_site_in(tmp_path):
    # The ordinary exome: DRAGEN called the in-design site, so the guard is satisfied by a
    # primary record and the merge goes on to fill the off-design hole from the supplement.
    posthoc = f'chr1\t{OFF_DESIGN}\t.\tC\tT,<NON_REF>\t80\t.\t.\tGT:DP:GQ\t0/1:44:80\n'

    result = _guard_then_merge(tmp_path, posthoc, dragen_records=DRAGEN_RECORDS)

    assert result.returncode == 0, result.stderr
    assert (tmp_path / 'merged.vcf.gz').exists()
    flags = _flags(tmp_path)
    assert flags[OFF_DESIGN].startswith('POSTHOC:')


# Single-copy chrX, where the two callers disagree about ploidy. XK:37694549 is a real
# GRCh38 defining site and is off-design on the Twist exome; XG:2748343 is a real one in
# PAR1, where a male sample is genuinely diploid and the fill must still work.
XK_SINGLE_COPY = 37694549
XG_IN_PAR1 = 2748343


def _chrx_beds(pos: int) -> tuple[str, str]:
    """Sites and off-design BEDs with one chrX hole at `pos`, beside the chr1 pair."""
    row = f'chrX\t{pos - 1}\t{pos}\n'
    return SITES_BED + row, OFF_DESIGN_BED + row


def test_a_posthoc_call_in_single_copy_chrx_is_not_filled(tmp_path):
    """HaplotypeCaller's diploid GT must not reach the merged VCF beside DRAGEN's haploid one.

    What this guards is not in this stage's own shell. HaplotypeCaller runs at its default
    ploidy of 2, so its record here says two chromosome copies while DRAGEN's chrX records say
    one, and rbceq2 answers a file claiming both by reporting the whole blood group
    Undetermined. A het post-hoc call is worse: it claims one copy, passes unchallenged, and
    flips the phenotype.
    """
    sites_bed, off_design_bed = _chrx_beds(XK_SINGLE_COPY)
    posthoc = f'chrX\t{XK_SINGLE_COPY}\t.\tG\tA,<NON_REF>\t410\t.\tSPARE=1\tGT:DP:GQ\t1/1:38:99\n'
    result = _merge(tmp_path, posthoc, sites_bed=sites_bed, off_design_bed=off_design_bed)

    assert result.returncode == 0, result.stderr
    merged = subprocess.run(  # noqa: S603
        ['bcftools', 'query', '-f', '%CHROM\t%POS\t[%GT]\n', str(tmp_path / 'merged.vcf.gz')],  # noqa: S607
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    assert f'chrX\t{XK_SINGLE_COPY}' not in merged
    # And the job says so, rather than dropping the site without a word.
    assert 'single-copy chrX left unfilled' in result.stderr


def test_a_posthoc_call_in_par1_is_still_filled(tmp_path):
    """The gate is bounded by PAR, so XG and CD99 keep their off-design recovery.

    A male sample has two copies of PAR1, so HaplotypeCaller's diploid GT is correct there
    and agrees with DRAGEN's. Gating the whole contig instead of the single-copy stretch
    would silently cost these sites for no ploidy reason at all.
    """
    sites_bed, off_design_bed = _chrx_beds(XG_IN_PAR1)
    posthoc = f'chrX\t{XG_IN_PAR1}\t.\tG\tC,<NON_REF>\t410\t.\tSPARE=1\tGT:DP:GQ\t0/1:38:99\n'
    result = _merge(tmp_path, posthoc, sites_bed=sites_bed, off_design_bed=off_design_bed)

    assert result.returncode == 0, result.stderr
    merged = subprocess.run(  # noqa: S603
        ['bcftools', 'query', '-f', '%CHROM\t%POS\t[%GT]\t%INFO/POSTHOC\n', str(tmp_path / 'merged.vcf.gz')],  # noqa: S607
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    assert f'chrX\t{XG_IN_PAR1}\t0/1\t{POSTHOC_CALLER}' in merged
    assert 'single-copy chrX left unfilled' not in result.stderr
