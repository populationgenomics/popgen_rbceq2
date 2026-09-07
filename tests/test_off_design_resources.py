"""The committed off-design defining sites, their manifest, and the generator that writes them.

An exome run fills post-hoc calls only at defining sites outside its capture design, and reads
that set from `resources/bg_off_design_sites.<design>.<genome>.bed`, committed once per design
by `scripts/gen_off_design_sites.py`. Committed derived data can go stale: regenerating the
defining sites without re-running the generator would leave a subtraction of the old sites in
place, and the fill would quietly widen or narrow. The manifest's `sites_md5` is the tie, and
`off_design.resource_path` refuses a row that does not match the shipped sites. This file holds
every committed row to the same checks, so the staleness is a red suite, not a wrong fill.

The generator tests run real bedtools over tiny inputs and are skipped without it; the
subtraction's interval semantics are test_off_design_subtraction's.
"""

import shutil
from pathlib import Path

import pytest

from popgen_rbceq2 import off_design, stage_support
from popgen_rbceq2.scripts import gen_off_design_sites
from tests.helpers import set_config

pytestmark = pytest.mark.fast

RESOURCES = Path(stage_support.blood_group_resource(off_design.MANIFEST_NAME)).parent
GENOME = 'GRCh38'
TWIST_KEY = 'exome_probesets_hg38/twist_vcgs_custom_exome_covered_targets_bed'
CREV2_KEY = 'exome_probesets_hg38/agilent_sureselect_clinical_research_exome_v2_covered_by_probes_bed'


def _manifest() -> dict[tuple[str, str], off_design.ManifestRow]:
    return off_design.read_manifest(RESOURCES / off_design.MANIFEST_NAME)


def _bed_rows(path: Path) -> list[str]:
    return path.read_text().splitlines()


# --- the committed resources ---


def test_both_validated_designs_are_committed():
    # The two designs the recall has been run against end to end (README). A run on either
    # must resolve without anyone running the generator first.
    assert {TWIST_KEY, CREV2_KEY} <= {key for key, _ in _manifest()}


def test_every_committed_row_was_subtracted_from_the_shipped_sites():
    # The staleness check itself, over every row: a bg_defining_sites regeneration that did not
    # re-run the generator fails here rather than filling a stale set.
    for (_, genome), row in _manifest().items():
        sites_md5 = off_design.file_md5(RESOURCES / f'bg_defining_sites.{genome}.bed')
        assert row.sites_md5 == sites_md5, f'{row.resource} was built from a different bg_defining_sites.{genome}.bed'


def test_every_committed_row_describes_the_file_beside_it():
    for (design_key, genome), row in _manifest().items():
        assert row.resource == off_design.resource_name(design_key, genome)
        resource = RESOURCES / row.resource
        assert resource.is_file(), f'{row.resource} is in the manifest but not committed'
        assert len(_bed_rows(resource)) == row.n_off_design
        assert len(_bed_rows(RESOURCES / f'bg_defining_sites.{genome}.bed')) == row.n_sites


def test_every_committed_resource_has_a_manifest_row():
    # The other direction: a BED dropped in by hand has no record of what it came from, and
    # resource_path would refuse it. Better to hear that here.
    committed = {path.name for path in RESOURCES.glob('bg_off_design_sites.*.bed')}
    assert committed == {row.resource for row in _manifest().values()}


def test_every_committed_resource_is_a_subset_of_the_shipped_sites():
    # `bedtools intersect -v` emits -a rows untouched, so each resource is rows of the sites
    # BED and nothing else; downstream reads it in that file's format.
    for (_, genome), row in _manifest().items():
        sites = set(_bed_rows(RESOURCES / f'bg_defining_sites.{genome}.bed'))
        off = _bed_rows(RESOURCES / row.resource)
        assert off, f'{row.resource} is empty: the design targets every defining site, which no exome design does'
        assert set(off) <= sites, f'{row.resource} has rows that are not defining sites'
        assert len(off) < len(sites), f'{row.resource} holds every defining site, so the design covered none'


def test_the_manifest_is_in_generator_order():
    # The generator writes rows sorted, so a hand edit that reorders them shows as a diff here
    # rather than as churn on the next generator run.
    path = RESOURCES / off_design.MANIFEST_NAME
    rows = _manifest()
    rewritten = path.parent / '.manifest_check.tsv'
    try:
        off_design.write_manifest(rewritten, rows)
        assert rewritten.read_text() == path.read_text()
    finally:
        rewritten.unlink(missing_ok=True)


# --- resolving one at graph build ---


def test_a_run_resolves_the_committed_file_for_its_design(shm_tmp_path):
    set_config({'references': {'genome_build': GENOME}}, shm_tmp_path / 'build.toml')

    path = Path(off_design.resource_path(TWIST_KEY))

    assert path == RESOURCES / off_design.resource_name(TWIST_KEY, GENOME)


def test_a_stale_resource_is_refused_at_graph_build(shm_tmp_path, monkeypatch):
    # The same check the suite makes above, applied where a run would read the file, so a
    # stale resource that somehow reached an image fails the submission and names the fix.
    set_config({'references': {'genome_build': GENOME}}, shm_tmp_path / 'build.toml')
    rows = _manifest()
    stale = {key: (row if key[0] != TWIST_KEY else _with_sites_md5(row, '0' * 32)) for key, row in rows.items()}
    monkeypatch.setattr(off_design, 'read_manifest', lambda _path: stale)

    with pytest.raises(ValueError, match='regenerated without it') as raised:
        off_design.resource_path(TWIST_KEY)

    assert off_design.GENERATOR in str(raised.value)


def test_a_resource_without_a_manifest_row_is_refused_at_graph_build(shm_tmp_path, monkeypatch):
    set_config({'references': {'genome_build': GENOME}}, shm_tmp_path / 'build.toml')
    rows = {key: row for key, row in _manifest().items() if key[0] != TWIST_KEY}
    monkeypatch.setattr(off_design, 'read_manifest', lambda _path: rows)

    with pytest.raises(FileNotFoundError, match='has no row'):
        off_design.resource_path(TWIST_KEY)


def _with_sites_md5(row: off_design.ManifestRow, sites_md5: str) -> off_design.ManifestRow:
    return off_design.ManifestRow(**{**row.__dict__, 'sites_md5': sites_md5})


# --- the generator ---

needs_bedtools = pytest.mark.skipif(not shutil.which('bedtools'), reason='bedtools is not on PATH')


def _resources_dir(tmp_path: Path) -> Path:
    out = tmp_path / 'resources'
    out.mkdir()
    (out / f'bg_defining_sites.{GENOME}.bed').write_text('chr1\t999\t1000\nchr1\t1499\t1500\nchr2\t9\t10\n')
    return out


@needs_bedtools
def test_the_generator_writes_the_resource_and_its_manifest_row(tmp_path):
    out = _resources_dir(tmp_path)
    design = tmp_path / 'design.bed'
    design.write_text('chr1\t990\t1010\n')

    row = gen_off_design_sites.generate('exome_probesets_hg38/test_bed', str(design), GENOME, out)

    assert (out / row.resource).read_text() == 'chr1\t1499\t1500\nchr2\t9\t10\n'
    assert row.resource == off_design.resource_name('exome_probesets_hg38/test_bed', GENOME)
    assert (row.n_sites, row.n_off_design) == (3, 2)
    assert row.sites_md5 == off_design.file_md5(out / f'bg_defining_sites.{GENOME}.bed')
    assert row.design_md5 == off_design.file_md5(design)
    assert row.bedtools_version.startswith('bedtools v')
    assert off_design.read_manifest(out / off_design.MANIFEST_NAME) == {('exome_probesets_hg38/test_bed', GENOME): row}


@needs_bedtools
def test_the_generator_adds_a_second_design_and_replaces_a_rerun_one(tmp_path):
    # One manifest for every design: a new design adds a row, re-running a design replaces its
    # row rather than duplicating it, and rows nobody touched are kept.
    out = _resources_dir(tmp_path)
    design_a = tmp_path / 'a.bed'
    design_a.write_text('chr1\t990\t1010\n')
    design_b = tmp_path / 'b.bed'
    design_b.write_text('chr1\t990\t1010\nchr1\t1490\t1510\n')

    gen_off_design_sites.generate('exome_probesets_hg38/a_bed', str(design_a), GENOME, out)
    gen_off_design_sites.generate('exome_probesets_hg38/b_bed', str(design_b), GENOME, out)
    design_a.write_text('chr1\t990\t1010\nchr2\t0\t100\n')
    row_a = gen_off_design_sites.generate('exome_probesets_hg38/a_bed', str(design_a), GENOME, out)

    manifest = off_design.read_manifest(out / off_design.MANIFEST_NAME)
    assert sorted(key for key, _ in manifest) == ['exome_probesets_hg38/a_bed', 'exome_probesets_hg38/b_bed']
    assert manifest[('exome_probesets_hg38/a_bed', GENOME)] == row_a
    assert row_a.n_off_design == 1
    assert manifest[('exome_probesets_hg38/b_bed', GENOME)].n_off_design == 1


@needs_bedtools
def test_the_generator_fails_on_a_design_the_subtraction_refuses(tmp_path):
    # The design checks are the generator's now, and a refused design must leave no resource
    # and no manifest row behind to be committed by mistake.
    out = _resources_dir(tmp_path)
    design = tmp_path / 'design.bed'
    design.write_text('1\t990\t1010\n')

    with pytest.raises(RuntimeError, match='every defining site is outside the capture design'):
        gen_off_design_sites.generate('exome_probesets_hg38/wrong_contigs_bed', str(design), GENOME, out)

    assert not list(out.glob('bg_off_design_sites.*'))


def test_the_generator_demands_the_sites_bed_it_subtracts_from(tmp_path):
    with pytest.raises(FileNotFoundError, match=r'gen_bg_resources\.py'):
        gen_off_design_sites.generate('exome_probesets_hg38/x_bed', str(tmp_path / 'd.bed'), GENOME, tmp_path)
