"""What the conversion step does to a haploid genotype on the sex chromosomes.

DRAGEN calls non-PAR chrX and chrY at their real ploidy in a male sample, writing a one-token
`GT=1` rather than the pseudo-diploid `1/1` some callers emit. XK, GATA1 and ATP11C are defined
there. Until rbceq2 2.4.3 a bare `1` crashed the run, so this pipe ended in a parameter-free
`bcftools +fixploidy` that rewrote it as `1|1`. rbceq2 2.4.4 reads the one-token form natively
and scores it as one chromosome copy, so the rewrite is gone.

These tests exist because putting it back would be silent. A diploidised hemizygous null is a
well-formed VCF and a well-formed call: rbceq2 reports `XK*N.16/XK*N.16`, indistinguishable in
the genotype TSV from a female homozygote, with the correct phenotype beside it. Nothing fails,
nothing logs, and the only evidence is a genotype string no consumer can challenge. So the
check has to be on the bytes the converter emits, not on anything downstream.

They run the real `bcftools view` from the real helper, because what is being asserted is
bcftools' own behaviour: that `view --trim-alt-alleles` leaves FORMAT/GT alone. A Python
stand-in would assert that this file believes that, which is not the same claim.

Skipped where bcftools is not installed. Local bcftools is expected to be the pinned image's
1.24, as for the other shell tests here.
"""

import shutil
import subprocess
from pathlib import Path

import pytest

from popgen_rbceq2.stages.blood_group_genotyping.filter_and_convert import _convert_commands

pytestmark = [
    pytest.mark.fast,
    pytest.mark.skipif(
        not all(shutil.which(tool) for tool in ('bcftools', 'bgzip')),
        reason='bcftools and bgzip are not both on PATH',
    ),
]

# chrX:37686068 is XK's defining coordinate in GRCh38, one of the three non-PAR blood-group
# loci that carry a haploid call in a male sample. chrX:2748343 is XG's, inside PAR1, where the
# same sample is genuinely diploid. Using the real coordinates keeps the fixture honest about
# which case is which, though nothing in the converter reads them.
XK_NON_PAR = 37686068
XG_PAR1 = 2748343

HEADER = """##fileformat=VCFv4.2
##contig=<ID=chrX,length=156040895>
##ALT=<ID=NON_REF,Description="Represents any possible alternative allele">
##INFO=<ID=END,Number=1,Type=Integer,Description="Block end position">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
"""


def _record(pos: int, gt: str, *, ref: str = 'G', alt: str = 'A') -> str:
    return f'chrX\t{pos}\t.\t{ref}\t{alt}\t200\tPASS\t.\tGT:DP:GQ\t{gt}:50:50\n'


def _convert(tmp_path: Path, body: str) -> list[str]:
    """Run the real conversion shell over one record set and return its genotype column.

    Writes the body as `merged.vcf.gz`, which is the file the converter reads, then runs the
    helper the stage runs and reads FORMAT/GT back out of what it wrote.
    """
    plain = tmp_path / 'merged.vcf'
    plain.write_text(HEADER + body)
    packed = tmp_path / 'merged.vcf.gz'
    with packed.open('wb') as out:
        subprocess.run(['bgzip', '-c', str(plain)], stdout=out, check=True)  # noqa: S603, S607

    script = 'set -euxo pipefail\n' + _convert_commands('converted.vcf.gz', cpu=1)
    subprocess.run(['bash', '-c', script], cwd=tmp_path, check=True, capture_output=True)  # noqa: S603, S607

    query = subprocess.run(  # noqa: S603
        ['bcftools', 'query', '-f', '[%GT]\n', str(tmp_path / 'converted.vcf.gz')],  # noqa: S607
        check=True,
        capture_output=True,
        text=True,
    )
    return query.stdout.split()


def test_a_haploid_genotype_survives_the_conversion_unchanged(tmp_path):
    """A one-token GT reaches rbceq2 as one token."""
    # The regression this guards: any step reinstated here that rewrites `1` to `1|1` turns a
    # hemizygous male null into something rbceq2 reports as a homozygote, with no error anywhere.
    assert _convert(tmp_path, _record(XK_NON_PAR, '1')) == ['1']


def test_a_haploid_reference_call_survives_the_conversion_unchanged(tmp_path):
    """The `0` half of the same case, which fixploidy also used to rewrite."""
    # rbceq2 2.4.4 reads a single-token `0` as no data rather than as a hemizygous reference,
    # so this is not merely the symmetric case: diploidising it to `0|0` would assert a
    # reference genotype the caller never made.
    assert _convert(tmp_path, _record(XK_NON_PAR, '0')) == ['0']


def test_a_diploid_genotype_in_par_survives_the_conversion_unchanged(tmp_path):
    """A PAR call is diploid in every sample and must stay that way."""
    assert _convert(tmp_path, _record(XG_PAR1, '0/1')) == ['0/1']


def test_the_converter_does_not_normalise_ploidy_across_records(tmp_path):
    """Mixed ploidy in, mixed ploidy out: the converter has no opinion about it.

    This is the state rbceq2 2.4.4 refuses, dropping the blood group to Undetermined rather
    than mis-rendering it, and the converter must not be the thing that hides it. Making the
    file self-consistent here would mean choosing a ploidy for a sample from its own calls,
    which is rbceq2's job and is done against the PAR table this stage does not carry.
    """
    body = _record(XK_NON_PAR, '1') + _record(XK_NON_PAR + 10, '1/1')
    assert _convert(tmp_path, body) == ['1', '1/1']


def test_the_non_ref_drop_still_applies_to_a_haploid_record(tmp_path):
    """Passing GT through does not mean passing the symbolic allele through."""
    # A haploid reference block is the ordinary shape on non-PAR chrX, so the <NON_REF> drop
    # has to keep working on exactly the records whose genotypes are now left alone.
    body = _record(XK_NON_PAR, '1') + _record(XK_NON_PAR + 10, '0', ref='C', alt='<NON_REF>')
    assert _convert(tmp_path, body) == ['1']


def test_the_conversion_shell_runs_no_bcftools_plugin():
    """The converter is two bcftools subcommands and no plugin."""
    # The behavioural tests above cover what the current shell does. This one is about what a
    # future edit puts back. A plugin is how ploidy gets rewritten -- `+fixploidy` was the one
    # here, `+setGT` would do it too -- so the shape is worth pinning, not just the old name.
    # Named separately from the behavioural tests because it needs no bcftools to fail.
    shell = _convert_commands('converted.vcf.gz', cpu=1)
    assert 'fixploidy' not in shell
    assert ' bcftools +' not in shell
