"""What the merge does with a post-hoc record that reaches past the hole it was kept for.

A record is selected for covering a hole and is selected whole, so one anchored in a hole can
extend over a neighbouring defining site the primary caller did call. Blood-group defining
sites are dense enough that this is the ordinary case near a capture edge. What the merge does
next differs by record type, and both halves rest on bcftools semantics rather than on any
Python here:

- a post-hoc **variant** reaching a defining site DRAGEN called is dropped, because keeping it
  hands rbceq2 two callers' alleles at one base and the QC reads that base as an ordinary PASS;
- a post-hoc **reference block** reaching one is kept, because it asserts nothing rbceq2 sees
  and dropping it would throw away the hole it was kept for.

So these tests run the real merge shell under real bcftools, then read the result with the real
extract format and the real QC functions. A Python stand-in would test a restatement of
`annotate -m`'s overlap rule rather than the rule, and that rule is the whole mechanism: `-m`
marks on a record's span, so it catches a deletion anchored on a hole that reaches a called
site one base away, and `INFO/END` is what separates a real reference block from the
`<NON_REF>` twin `norm -m -any` splits off a variant. The mark is the covered defining sites
rather than the DRAGEN records' spans, which is the last test here.

Skipped where bcftools is not installed. Local bcftools is expected to be the pinned image's
1.24, so the semantics checked here are the ones the job will meet.
"""

import shutil
import subprocess
from pathlib import Path

import pytest

from popgen_rbceq2.jobs.rbceq2_call_qc_job import flag_site, parse_extract, resolve_coverage
from popgen_rbceq2.scripts.bg_db import DefiningSite
from popgen_rbceq2.stages.blood_group_genotyping.filter_and_convert import (
    _EXTRACT_FORMAT,
    _POSTHOC_HEADER_LINE,
    _merge_posthoc_commands,
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
# What SelectOffDesignDefiningSites hands this job: the defining sites the design never
# targeted, which here is the second one and not the first. The subtraction itself is tested in
# test_exome_design_gate; this file starts from its result.
OFF_DESIGN_BED = f'chr1\t{OFF_DESIGN - 1}\t{OFF_DESIGN}\n'

# DRAGEN calls the in-design site and has no record at all at the off-design one, which is the
# capture edge this feature exists for.
DRAGEN_RECORDS = f'chr1\t{IN_DESIGN}\t.\tA\tG,<NON_REF>\t200\tPASS\tSPARE=1\tGT:DP:GQ\t0/1:50:50\n'

# A gVCF header declares far more tags than any one record carries, which matters: the merge's
# `annotate -x '^INFO/END,...'` reads the header, and errors outright if the header has nothing
# outside its keep list. SPARE stands in for the BaseQRankSum, MLEAC and friends a real
# HaplotypeCaller header declares.
HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=250000000>
##ALT=<ID=NON_REF,Description="Represents any possible alternative allele">
##INFO=<ID=END,Number=1,Type=Integer,Description="Block end position">
##INFO=<ID=SPARE,Number=1,Type=Integer,Description="A tag nothing downstream reads">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum depth over the block">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
"""


def _bgzip(tmp_path: Path, name: str, body: str, *, index: bool) -> Path:
    """Write a VCF body under HEADER and bgzip it, indexing only when asked."""
    plain = tmp_path / name
    plain.write_text(HEADER + body)
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
) -> subprocess.CompletedProcess[str]:
    """Run the real merge shell over one post-hoc record set, in tmp_path."""
    posthoc = _bgzip(tmp_path, 'posthoc.g.vcf', posthoc_records, index=True)
    sites = tmp_path / 'sites.bed'
    sites.write_text(sites_bed)
    off_design = tmp_path / 'off_design_defining_sites.bed'
    off_design.write_text(off_design_bed)
    # What the caller leaves behind for the fragment: the POSTHOC declaration both sides of
    # the concat have to carry, and the bg-regions intermediate. Not indexed, because the
    # fragment's own first line indexes it.
    (tmp_path / 'posthoc_hdr.txt').write_text(f'{_POSTHOC_HEADER_LINE}\n')
    _bgzip(tmp_path, 'dragen.vcf', dragen_records, index=False)

    script = 'set -euo pipefail\n' + _merge_posthoc_commands(
        str(posthoc),
        str(sites),
        str(off_design),
        'exome_probesets_hg38/test_design_bed',
        cpu=1,
    )
    return subprocess.run(  # noqa: S603
        ['bash', '-c', script],  # noqa: S607
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=True,
    )


def _flags(tmp_path: Path, *, sites_bed: str = SITES_BED) -> dict[int, str]:
    """Extract merged.vcf.gz the way the stage does and flag every defining site."""
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
    sites = [
        DefiningSite(chrom='chr1', pos=int(line.split('\t')[2]), ref='C', alt='T', kind='var', system='TEST')
        for line in sites_bed.splitlines()
    ]
    return {
        site.pos: flag_site(site, resolve_coverage(records, site.chrom, site.pos), 10, 20) or 'PASS' for site in sites
    }


def test_a_posthoc_deletion_reaching_a_called_base_is_dropped_and_the_site_stays_a_hole(tmp_path):
    # The finding this file exists for. Kept, the deletion would remove the base DRAGEN called
    # a SNP on, in the file rbceq2 reads, and the QC would still say PASS there because
    # resolve_coverage prefers the primary record. Dropped, the hole it filled goes back to
    # NOCOV, which is the honest answer: DRAGEN spoke at the base this record reaches.
    deletion = f'chr1\t{OFF_DESIGN}\t.\tCATGAAAAA\tC,<NON_REF>\t60\t.\tSPARE=1\tGT:DP:GQ\t0/1:44:80\n'

    _merge(tmp_path, deletion)

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
    assert '2 dropped for' in stderr


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
    flags = _flags(tmp_path, sites_bed=sites)
    assert flags[far_site] == f'DEL+POSTHOC:1:{far_site}(C>T,src=gatk-hc-4.6.2.0,del=CATGA>C,GT=0/1,DP=44,GQ=80)'
