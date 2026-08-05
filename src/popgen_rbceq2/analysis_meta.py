"""The ``update_analysis_meta`` callbacks that turn a written TSV into Metamist Analysis meta.

Each function takes the path of the output its Analysis is registered against and returns extra
``Analysis.meta``. They run in the Metamist status-reporter PythonJob, which depends on the job
that wrote the file, so the file is there to be read.

They are module-level functions, not methods: cpg_flow dill-pickles the callable into that job,
and a bound method would carry ``self`` and fail there, long after the compute has run.

None of them set ``meta.stage`` — ``stage_support.wire`` adds it from the stage class name, so a
renamed stage cannot leave a stale literal behind. See its docstring.
"""

from typing import Any

import cpg_utils
from cpg_utils.config import config_retrieve, genome_build

from popgen_rbceq2.constants import RBCEQ2_VERSION


def parse_single_row_rbceq2_tsv(text: str) -> dict[str, str]:
    """Parse a one-data-row rbceq2 TSV into ``{blood_group_system: value}``.

    rbceq2 per-sample TSVs are a header row of blood-group systems plus a single
    data row; the first column is the sample id and is dropped.
    """
    rows = [line.split('\t') for line in text.splitlines() if line.strip()]
    if len(rows) < 2:
        raise ValueError(
            'Missing rows in Rbceq2 output. Single sample Rbceq2 outputs have two rows; a header row and sample row'
        )
    header, values = rows[0], rows[1]
    return dict(zip(header[1:], values[1:], strict=True))


def blood_group_calls(output: str) -> dict[str, Any]:
    """Meta for the per-SG rbceq2 calls Analysis.

    Args:
        output: The ``<sg>.geno.tsv`` path; the alphanumeric phenotype TSV sits beside it.

    Returns:
        The inferred blood-group genotypes and phenotypes, plus the tool version and reference
        build they were called with.
    """
    genotypes: dict[str, str] = parse_single_row_rbceq2_tsv(cpg_utils.to_path(output).read_text())
    pheno_path: cpg_utils.Path = cpg_utils.to_path(output.replace('.geno.tsv', '.pheno_alphanumeric.tsv'))
    phenotypes: dict[str, str] = parse_single_row_rbceq2_tsv(pheno_path.read_text())
    return {
        'rbceq2_version': RBCEQ2_VERSION,
        'reference_genome': genome_build(),
        'blood_group_genotypes': genotypes,
        'blood_group_phenotypes': phenotypes,
    }


def call_qc(output: str) -> dict[str, Any]:
    """Meta for the per-SG blood-group call QC Analysis.

    Args:
        output: Path to the ``<sg>.qc.tsv``, in the same wide one-row layout as rbceq2's own
            TSVs.

    Returns:
        The Analysis meta. ``blood_group_qc`` maps each system to ``PASS``, a
        semicolon-joined ``LOWQ``/``NOCOV`` flag naming the failing site and its DP and GQ,
        or ``NA`` for a system rbceq2 called that has no defining site in the map, e.g.
        ``{'JK': 'PASS', 'VEL': 'LOWQ:1:3774964(A>G,DP=8,GQ=45)', 'FUT2': 'NA'}``. The
        thresholds are recorded alongside, since a flag means nothing without them.
    """
    cfg = 'flag_blood_group_call_qc'
    return {
        'reference_genome': genome_build(),
        'min_depth': config_retrieve(['workflow', cfg, 'min_depth'], 10),
        'min_gq': config_retrieve(['workflow', cfg, 'min_gq'], 20),
        'blood_group_qc': parse_single_row_rbceq2_tsv(cpg_utils.to_path(output).read_text()),
    }


def cohort_calls(output: str) -> dict[str, str]:
    """Meta for the combined-cohort rbceq2 Analysis.

    Args:
        output: The combined geno TSV.

    Returns:
        The sibling phenotype and QC TSV paths, which are otherwise undiscoverable: the
        Analysis is registered against the geno TSV alone.
    """
    return {
        'pheno_numeric_path': output.replace('.geno.tsv', '.pheno_numeric.tsv'),
        'pheno_alphanumeric_path': output.replace('.geno.tsv', '.pheno_alphanumeric.tsv'),
        'qc_path': output.replace('.geno.tsv', '.qc.tsv'),
    }
