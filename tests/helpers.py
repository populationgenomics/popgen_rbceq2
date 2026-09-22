"""Helper functions for the tests."""

import pathlib
import subprocess
from typing import Any, Protocol, runtime_checkable

import toml
from cpg_utils import Path
from cpg_utils.config import set_config_paths


@runtime_checkable
class IDictRepresentable(Protocol):
    def as_dict(self) -> dict[str, Any]: ...


class TomlAnyPathEncoder(toml.TomlEncoder):
    """Support for CPG path objects in TOML.

    A CPG path is either a regular pathlib path or a cloud path.
    """

    def dump_value(self, v):
        if isinstance(v, Path):
            v = str(v)
        return super().dump_value(v)


def set_config(
    config: str | dict[str, Any] | IDictRepresentable,
    path: Path,
    merge_with: list[Path] | None = None,
) -> None:
    """Write a config to `path` and point `CPG_CONFIG_PATH` at it.

    That environment variable is how `cpg_utils` finds a config. If `merge_with` is provided,
    the config is merged with the configs at the given paths. Merging happens right to left, so
    values in the right config override values in the left one.

    Args:
        config (str | dict[str, Any] | IDictRepresentable):
            A valid TOML string, a dictionary to be converted to TOML, or an object
            which implements the `IDictRepresentable` protocol.

        path (Path):
            Path to write the config to.

        merge_with (list[Path] | None, optional):
            A list of paths to merge with the config. Merging happens right to left,
            so that values in the right config will override values in the left config.
            Defaults to `None`.
    """
    with path.open('w') as f:
        if isinstance(config, dict):
            toml.dump(config, f, encoder=TomlAnyPathEncoder())
        elif isinstance(config, IDictRepresentable):
            toml.dump(config.as_dict(), f, encoder=TomlAnyPathEncoder())
        elif isinstance(config, str):
            f.write(config)
        else:
            raise TypeError(f'Expected config to be a string, dict, or IDictRepresentable, butgot {type(config)}')

        f.flush()

    return set_config_paths([*[str(s) for s in (merge_with or [])], str(path)])


# A gVCF header for the shell tests, declaring more tags than any one record carries. That
# matters for the merge: `annotate -x '^INFO/END,...'` reads the header and errors outright if
# it has nothing outside its keep list, which is what SPARE stands in for. The tests that do
# not run annotate are unaffected by the extra declarations.
VCF_HEADER = """##fileformat=VCFv4.2
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


def write_bgzipped_vcf(
    tmp_path: pathlib.Path,
    name: str,
    body: str,
    *,
    index: bool = False,
    extra_header: str = '',
    sample: str = 'SAMPLE1',
) -> pathlib.Path:
    """Write a VCF body under VCF_HEADER and bgzip it, indexing only when asked.

    Args:
        tmp_path: Directory to write into.
        name: Uncompressed file name; the bgzipped file gets `.gz` on the end.
        body: Record lines, no header.
        index: Whether to write a `.tbi` beside it.
        extra_header: One more `##` line, placed before the column header.
        sample: Name for the one sample column.

    Returns:
        Path to the bgzipped file.
    """
    plain = tmp_path / name
    header = VCF_HEADER if not extra_header else VCF_HEADER.replace('#CHROM', f'{extra_header}\n#CHROM')
    header = header.replace('\tSAMPLE1\n', f'\t{sample}\n')
    plain.write_text(header + body)
    packed = tmp_path / f'{name}.gz'
    with packed.open('wb') as out:
        subprocess.run(['bgzip', '-c', str(plain)], stdout=out, check=True)  # noqa: S603, S607
    if index:
        subprocess.run(['bcftools', 'index', '-t', str(packed)], check=True)  # noqa: S603, S607
    return packed
