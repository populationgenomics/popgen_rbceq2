"""That the conversion step leaves FORMAT/GT alone, whatever ploidy it is written at.

Runs the real helper under real bcftools, and is skipped where bcftools is not installed.
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

# XK's defining coordinate in GRCh38, one of the three loci that carry a one-token call in a
# single-copy sample. Nothing in the converter reads the coordinate.
XK_NON_PAR = 37686068

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


@pytest.mark.parametrize('gts', [['1'], ['0'], ['0/1'], ['1', '0/1']])
def test_the_converter_passes_genotypes_through(tmp_path, gts):
    """Every GT the converter is given comes out of it byte-identical."""
    # `0` is not the inert case it looks like. rbceq2 drops a one-token `0` as hom-ref exactly
    # as it drops `0/0`, so the row contributes no allele either way, but ploidy inference runs
    # before that drop and reads one token as single-copy evidence where two read as two.
    #
    # The mixed pair is deliberate: rbceq2 refuses a file claiming both one chromosome copy and
    # two, and the converter must not be what hides that. Making the file self-consistent here
    # would mean choosing a ploidy from the sample's own calls, which is the merge's job.
    body = ''.join(_record(XK_NON_PAR + 10 * i, gt) for i, gt in enumerate(gts))
    assert _convert(tmp_path, body) == gts


def test_a_non_ref_record_is_dropped(tmp_path):
    """A <NON_REF>-only record is not in the converted VCF; the real call beside it is."""
    body = _record(XK_NON_PAR, '1') + _record(XK_NON_PAR + 10, '0', ref='C', alt='<NON_REF>')
    assert _convert(tmp_path, body) == ['1']
