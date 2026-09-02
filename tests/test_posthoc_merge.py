"""Tests for the post-hoc merge's two awk programs, in `stages/.../filter_and_convert.py`.

These run the real awk, not a Python re-implementation of it. The programs are the only logic
in this pipeline written outside Python, so a Python stand-in would test the wrong thing: what
matters is that `awk` itself, given what `bcftools query` really emits, decides holes the same
way the QC job's `GvcfRecord.covers` does. The two disagreeing is the failure that would put a
post-hoc record over a site DRAGEN had called, or leave a real hole unfilled.

No bcftools here, so nothing in this file needs an image or a cloud: the covered-spans BED is
rendered from `GvcfRecord`s using the same `%CHROM/%POS0/%END` fields bcftools would write.
"""

import shutil
import subprocess

import pytest

from popgen_rbceq2.jobs.rbceq2_call_qc_job import EXTRACT_COLUMNS, GvcfRecord
from popgen_rbceq2.stages.blood_group_genotyping.filter_and_convert import (
    _EXTRACT_FORMAT,
    _POSTHOC_HEADER_LINE,
    _TAG_POSTHOC_AWK,
    _UNCOVERED_SITES_AWK,
)

pytestmark = [
    pytest.mark.fast,
    pytest.mark.skipif(not shutil.which('awk'), reason='awk is not on PATH'),
]


def test_the_extract_format_has_one_field_per_parsed_column():
    """The bcftools query format and the parser's columns stay the same length."""
    # A field added to one and not the other is caught only when a real job runs and
    # parse_extract rejects the file it was handed, which is a whole batch late.
    fields = [f for f in _EXTRACT_FORMAT.removesuffix(r'\n').split(r'\t') if f]
    assert len(fields) == len(EXTRACT_COLUMNS)


def test_the_extract_format_declares_the_posthoc_tag_it_reads():
    """Anything the extract reads as INFO/POSTHOC must be declared on the intermediate."""
    # bcftools query fails outright on a tag the header does not define, rather than rendering
    # `.`, so the declaration is what keeps genome runs (where no record carries the tag)
    # working at all. Losing it breaks every run that has nothing to merge.
    assert '%INFO/POSTHOC' in _EXTRACT_FORMAT
    assert 'ID=POSTHOC' in _POSTHOC_HEADER_LINE


def _record(chrom, pos, ref='A', alt='G', end=None) -> GvcfRecord:
    return GvcfRecord(chrom=chrom, pos=pos, ref=ref, alt=alt, end=end, gt='0/0', dp=30, gq=40, min_dp=25)


def _covered_bed(records: list[GvcfRecord]) -> str:
    r"""Render records the way `bcftools query -f '%CHROM\t%POS0\t%END\n'` would.

    %POS0 is the 0-based position and %END is POS + rlen - 1, where rlen comes from INFO/END
    on a reference block and from len(REF) otherwise. That is the record's whole span, which
    is what `GvcfRecord.covers` tests against.
    """
    return ''.join(f'{r.chrom}\t{r.pos - 1}\t{max(r.end or r.pos, r.pos + len(r.ref) - 1)}\n' for r in records)


def _sites_bed(sites: list[tuple[str, int]]) -> str:
    """Render 1-based defining sites as the committed BED does, 0-based and half-open."""
    return ''.join(f'{chrom}\t{pos - 1}\t{pos}\n' for chrom, pos in sites)


def _run_uncovered(tmp_path, records: list[GvcfRecord], sites: list[tuple[str, int]]) -> list[tuple[str, int]]:
    """Run the hole-finding awk and return the sites it reports as uncovered."""
    covered = tmp_path / 'covered.bed'
    covered.write_text(_covered_bed(records))
    site_file = tmp_path / 'sites.bed'
    site_file.write_text(_sites_bed(sites))
    out = subprocess.run(  # noqa: S603
        ['awk', _UNCOVERED_SITES_AWK, str(covered), str(site_file)],  # noqa: S607
        capture_output=True,
        text=True,
        check=True,
    )
    return [(line.split('\t')[0], int(line.split('\t')[2])) for line in out.stdout.splitlines()]


def test_a_site_inside_a_reference_block_is_not_a_hole(tmp_path):
    records = [_record('chr1', 900, ref='T', alt='<NON_REF>', end=1100)]
    assert _run_uncovered(tmp_path, records, [('chr1', 1000)]) == []


def test_a_site_with_no_record_at_all_is_a_hole(tmp_path):
    records = [_record('chr1', 900, ref='T', alt='<NON_REF>', end=1100)]
    assert _run_uncovered(tmp_path, records, [('chr1', 2000)]) == [('chr1', 2000)]


def test_a_site_on_another_contig_is_a_hole(tmp_path):
    # Records are keyed by contig in the awk. A contig with no covered spans at all must fall
    # through to "print", not raise or silently count as covered.
    records = [_record('chr1', 900, ref='T', alt='<NON_REF>', end=1100)]
    assert _run_uncovered(tmp_path, records, [('chr2', 1000)]) == [('chr2', 1000)]


@pytest.mark.parametrize('pos', [900, 1100])
def test_the_ends_of_a_block_are_covered(tmp_path, pos):
    # Off-by-one at either edge would either fill a site that has a DRAGEN record, or leave a
    # real hole unfilled. Both ends are inclusive, matching covers().
    records = [_record('chr1', 900, ref='T', alt='<NON_REF>', end=1100)]
    assert _run_uncovered(tmp_path, records, [('chr1', pos)]) == []


@pytest.mark.parametrize('pos', [899, 1101])
def test_one_base_outside_a_block_is_a_hole(tmp_path, pos):
    records = [_record('chr1', 900, ref='T', alt='<NON_REF>', end=1100)]
    assert _run_uncovered(tmp_path, records, [('chr1', pos)]) == [('chr1', pos)]


def test_a_deletions_ref_span_covers_the_bases_it_removes(tmp_path):
    # A deletion is not a hole: the QC has something to say about it (DEL), and calling the
    # site again from the CRAM would replace a real finding with a second opinion.
    records = [_record('chr1', 3774961, ref='CATGA', alt='C')]
    assert _run_uncovered(tmp_path, records, [('chr1', 3774964)]) == []
    assert _run_uncovered(tmp_path, records, [('chr1', 3774966)]) == [('chr1', 3774966)]


def test_the_awk_agrees_with_the_qc_jobs_coverage_rule(tmp_path):
    """Over a mixed record set, awk marks exactly the sites `covers` says are uncovered."""
    # The property that matters. A site this calls a hole is a site the QC would otherwise
    # flag NOCOV, so the two definitions have to be the same one. They are written in
    # different languages, in different files, and nothing but this test ties them together.
    records = [
        _record('chr1', 900, ref='T', alt='<NON_REF>', end=1100),
        _record('chr1', 3774961, ref='CATGA', alt='C'),
        _record('chr1', 5000, ref='C', alt='G'),
        _record('chr2', 100, ref='A', alt='<NON_REF>', end=100),
        _record('chr2', 7000, ref='GGGG', alt='G'),
    ]
    sites = (
        [('chr1', pos) for pos in (899, 900, 1000, 1100, 1101, 3774960, 3774961, 3774964, 3774965, 3774966, 5000, 5001)]
        + [('chr2', pos) for pos in (99, 100, 101, 6999, 7000, 7003, 7004)]
        + [('chr3', 1)]
    )

    expected = [(chrom, pos) for chrom, pos in sites if not any(r.chrom == chrom and r.covers(pos) for r in records)]
    assert _run_uncovered(tmp_path, records, sites) == expected
    # The set is a real mix, not accidentally all-covered or all-uncovered.
    assert 0 < len(expected) < len(sites)


def _run_tag(tmp_path, vcf_body: str, tag: str = 'POSTHOC=gatk-hc-4.6.2.0') -> list[str]:
    """Run the INFO-tagging awk over VCF text and return its output lines."""
    vcf = tmp_path / 'in.vcf'
    vcf.write_text(vcf_body)
    out = subprocess.run(  # noqa: S603
        ['awk', '-v', 'OFS=\t', '-v', f'tag={tag}', _TAG_POSTHOC_AWK, str(vcf)],  # noqa: S607
        capture_output=True,
        text=True,
        check=True,
    )
    return out.stdout.splitlines()


def test_tagging_replaces_an_empty_info_column(tmp_path):
    # `.` is the empty INFO column. Appending to it would produce `.;POSTHOC=...`, which is
    # not valid INFO and which bcftools rejects on the next read.
    body = 'chr1\t2100\t.\tA\tT\t.\t.\t.\tGT:DP\t0/1:33\n'
    assert _run_tag(tmp_path, body) == ['chr1\t2100\t.\tA\tT\t.\t.\tPOSTHOC=gatk-hc-4.6.2.0\tGT:DP\t0/1:33']


def test_tagging_appends_to_an_existing_info_column(tmp_path):
    body = 'chr1\t950\t.\tG\t<NON_REF>\t.\t.\tEND=2050\tGT:DP\t0/0:45\n'
    expected = 'chr1\t950\t.\tG\t<NON_REF>\t.\t.\tEND=2050;POSTHOC=gatk-hc-4.6.2.0\tGT:DP\t0/0:45'
    assert _run_tag(tmp_path, body) == [expected]


def test_tagging_leaves_header_lines_untouched(tmp_path):
    # The header carries the ## lines and the #CHROM row. Rewriting field 8 of #CHROM would
    # rename the INFO column and corrupt the file.
    body = '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\n'
    assert _run_tag(tmp_path, body) == [
        '##fileformat=VCFv4.2',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1',
    ]
