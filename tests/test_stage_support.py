"""Locks the conventions stage_support enforces, and the DAG declared in stages/pipeline.py.

A stage's `[workflow.<section>]` key is derived from its class name; a rename or a new stage
that desyncs the two would silently fall through to in-code defaults, so these tests fail
loudly instead. Same for a stage that is wired but never reached, and for a stage that reads an
upstream it forgot to declare.
"""

import argparse
import ast
import importlib.resources
import inspect
import pathlib
import re
import shlex
import tomllib
from collections.abc import Callable, Mapping
from types import SimpleNamespace
from typing import Any, cast

import cpg_flow.stage
import cpg_utils
import cpg_utils.config
import pytest

from popgen_rbceq2 import run_workflow, stage_support
from popgen_rbceq2.stages import pipeline
from tests import helpers

pytestmark = pytest.mark.fast


# --- wire: Metamist registration is declared here, so a typo is an error ---


class _Impl(cpg_flow.stage.CohortStage):
    """Minimal undecorated stage; wire() is what turns one of these into a DAG node."""

    def expected_outputs(self, cohort) -> dict:  # noqa: ARG002
        return {'vcf': cpg_utils.to_path('gs://bucket/out.vcf.gz')}

    def queue_jobs(self, cohort, inputs) -> None:  # noqa: ARG002
        return None


def _meta(output: str) -> dict:
    return {'n_systems': 42, 'path': output}


def _hook(wired) -> Callable[[str], dict]:
    """The update_analysis_meta cpg_flow will call, asserted present."""
    hook = wired().update_analysis_meta
    assert hook is not None
    return hook


def _as_stage(name: str) -> Any:
    """A stand-in carrying only the attribute stage_support reads off a stage: its name."""
    return cast('Any', SimpleNamespace(name=name))


@pytest.fixture
def _workflow(mocker, tmp_path) -> Any:
    """Constructing a wired stage reads the active workflow and config; neither exists here."""
    helpers.set_config({'workflow': {'name': 'popgen_rbceq2', 'dataset': 'ourdna'}}, tmp_path / 'config.toml')
    mock_wf = mocker.MagicMock()
    mock_wf.name = 'popgen_rbceq2'
    mock_wf.status_reporter = None
    mocker.patch('cpg_flow.stage.get_workflow', return_value=mock_wf)
    return mock_wf


@pytest.mark.usefixtures('_workflow')
def test_wire_records_the_stage_name_without_anyone_typing_it():
    # The whole point: cpg_flow hands update_analysis_meta only the output path, so a stage
    # that wants its own name in Analysis.meta would otherwise hardcode a string literal that
    # drifts the moment the class is renamed.
    wired = stage_support.wire(
        type('CallSomething', (_Impl,), {}),
        analysis_type='qc',
        analysis_keys=['vcf'],
        update_analysis_meta=_meta,
    )
    assert _hook(wired)('gs://bucket/out.tsv') == {
        'stage': 'CallSomething',
        'n_systems': 42,
        'path': 'gs://bucket/out.tsv',
    }


@pytest.mark.usefixtures('_workflow')
def test_wire_records_the_stage_name_when_the_stage_adds_no_meta_of_its_own():
    wired = stage_support.wire(type('CallPlain', (_Impl,), {}), analysis_type='qc', analysis_keys=['vcf'])
    assert _hook(wired)('gs://bucket/out.tsv') == {'stage': 'CallPlain'}


@pytest.mark.usefixtures('_workflow')
def test_wire_keeps_the_class_name_even_when_the_stage_sets_stage_itself():
    # The pre-wire house style hardcoded {'stage': 'rbceq2'}. Copy such a function into a new
    # stage, forget to delete that key, and a merge in the other order would file every
    # Analysis under someone else's stage name with no error anywhere.
    def stale(output: str) -> dict:  # noqa: ARG001
        return {'stage': 'rbceq2', 'n': 1}

    wired = stage_support.wire(
        type('CallCopied', (_Impl,), {}),
        analysis_type='qc',
        analysis_keys=['vcf'],
        update_analysis_meta=stale,
    )
    assert _hook(wired)('gs://b/o.tsv') == {'stage': 'CallCopied', 'n': 1}


@pytest.mark.usefixtures('_workflow')
def test_wire_meta_hook_rejects_a_non_dict_return():
    # A function that falls off the end returns None, and `dict | None` is a TypeError with no
    # indication of which stage produced it.
    wired = stage_support.wire(
        type('CallNoReturn', (_Impl,), {}),
        analysis_type='qc',
        analysis_keys=['vcf'],
        # Cast past the type checker: a meta function that falls off the end is exactly
        # what the runtime guard exists for, and only untyped callers can produce it.
        update_analysis_meta=cast('Callable[[str], dict]', lambda output: None),  # noqa: ARG005
    )
    with pytest.raises(TypeError, match='CallNoReturn returned NoneType'):
        _hook(wired)('gs://b/o.tsv')


@pytest.mark.usefixtures('_workflow')
def test_wire_adds_no_meta_hook_when_nothing_is_registered():
    # analysis_type unset means no Analysis at all, so cpg_flow must not be handed a hook.
    assert stage_support.wire(type('CallNone', (_Impl,), {}))().update_analysis_meta is None


def test_wire_rejects_analysis_keys_without_a_type():
    # cpg_flow gates registration on analysis_type, so the keys would be silently ignored.
    with pytest.raises(ValueError, match='without analysis_type'):
        stage_support.wire(type('CallOrphanKeys', (_Impl,), {}), analysis_keys=['vcf'])


def test_wire_rejects_a_meta_function_that_takes_self():
    # A method binds self to the output path and blows up inside the Metamist status job, long
    # after the compute has run and only where the status reporter is enabled.
    with pytest.raises(TypeError, match='callable with one argument'):
        stage_support.wire(
            type('CallMethodMeta', (_Impl,), {}),
            analysis_type='qc',
            analysis_keys=['vcf'],
            update_analysis_meta=cast('Callable[[str], dict]', lambda self, output: {}),  # noqa: ARG005
        )


def test_wire_rejects_an_unknown_registration_argument():
    # The reason for moving these off the class: Python checks the keyword for us.
    with pytest.raises(TypeError):
        stage_support.wire(type('CallTypo', (_Impl,), {}), analysis_typ='qc')


@pytest.mark.usefixtures('_workflow')
def test_wire_passes_other_stage_options_through():
    wired = stage_support.wire(type('CallTolerant', (_Impl,), {}), analysis_type='qc', tolerate_missing_output=True)
    assert wired().tolerate_missing_output


@pytest.mark.usefixtures('_workflow')
def test_wire_meta_hook_survives_the_pickling_cpg_flow_does_to_it():
    # cpg_flow ships this callable into a Hail PythonJob, so it is dill-pickled and run in
    # another container. A closure over the class would be the obvious implementation and a
    # riskier one to serialise.
    dill = pytest.importorskip('dill')
    cls = type('CallPickle', (_Impl,), {})
    wired = stage_support.wire(cls, analysis_type='qc', analysis_keys=['vcf'], update_analysis_meta=_meta)
    assert dill.loads(dill.dumps(_hook(wired)))('gs://b/o.tsv')['stage'] == 'CallPickle'


# --- REQUESTED_STAGES: the one list that decides what actually runs ---


def _wired_stages() -> dict[str, Any]:
    """Every stage wired into the DAG, by name.

    `@stage` returns a function wrapping the class under functools.wraps, so `__wrapped__`
    pointing at a Stage subclass is what distinguishes a wired stage from anything else in the
    pipeline module's namespace. pipeline is the only place wire() is called, so it holds all
    of them.
    """
    wired = {}
    for name, obj in vars(pipeline).items():
        wrapped = getattr(obj, '__wrapped__', None)
        if callable(obj) and isinstance(wrapped, type) and issubclass(wrapped, cpg_flow.stage.Stage):
            wired[name] = obj
    return wired


@pytest.mark.usefixtures('_workflow')
def test_every_wired_stage_is_reachable_from_requested_stages():
    # A stage that is neither requested nor required by something requested never enters the
    # DAG: imports succeed, tests pass, dry-run passes, and it silently produces nothing. This
    # is the one step in adding a stage that has no loud failure of its own.
    wired = _wired_stages()
    unknown = [s.__name__ for s in pipeline.REQUESTED_STAGES if s.__name__ not in wired]
    assert not unknown, f'{unknown} are in REQUESTED_STAGES but are not wired stages'

    reachable: set[str] = set()
    pending = [s.__name__ for s in pipeline.REQUESTED_STAGES]
    while pending:
        name = pending.pop()
        if name in reachable:
            continue
        reachable.add(name)
        pending.extend(dep.__name__ for dep in wired[name]().required_stages_classes)

    assert not set(wired) - reachable, (
        f'{sorted(set(wired) - reachable)} are wired but unreachable from REQUESTED_STAGES. '
        f'Add each to REQUESTED_STAGES in stages/pipeline.py, or to the requires= of a stage '
        f'already there.'
    )


# Which positional argument of each StageInput accessor names the upstream stage. The
# *_by_target variants drop the target, so the stage moves to the front.
_STAGE_ARG_INDEX = {
    'as_path': 1,
    'as_str': 1,
    'as_dict': 1,
    'as_path_by_target': 0,
    'as_dict_by_target': 0,
    'as_path_dict_by_target': 0,
}

_STAGE_SOURCES = sorted(
    # One level down from pipeline.py, so a new stage subpackage is picked up on its own.
    pathlib.Path(pipeline.__file__).parent.glob('*/*.py'),
)
assert len(_STAGE_SOURCES) > 0


def test_stage_sources_cover_every_wired_stage():
    # The glob pins where implementations live. A stage module it misses — deeper nesting, or
    # a module beside pipeline.py — would drop out of
    # test_stages_only_read_from_stages_they_require without failing anything.
    sources = {source.resolve() for source in _STAGE_SOURCES}
    for name, stage_decorator in _wired_stages().items():
        impl_file = pathlib.Path(inspect.getfile(stage_decorator.__wrapped__)).resolve()
        assert impl_file in sources, f'{name} is implemented in {impl_file}, which _STAGE_SOURCES misses'


def _stages_read_by(class_node: ast.ClassDef) -> set[str]:
    """The stages this class body passes to inputs.as_*, read from the source.

    A stage is named `<module>.<Class>` at the call site, so the attribute is what identifies
    it — the module alias in front carries no information the DAG cares about.
    """
    names: set[str] = set()
    for node in ast.walk(class_node):
        func = getattr(node, 'func', None)
        if not (isinstance(node, ast.Call) and isinstance(func, ast.Attribute)):
            continue
        if func.attr not in _STAGE_ARG_INDEX or not (isinstance(func.value, ast.Name) and func.value.id == 'inputs'):
            continue
        index = _STAGE_ARG_INDEX[func.attr]
        stage_arg = node.args[index] if len(node.args) > index else None
        if isinstance(stage_arg, ast.Attribute):
            names.add(stage_arg.attr)
        elif isinstance(stage_arg, ast.Name):
            names.add(stage_arg.id)
    return names


@pytest.mark.usefixtures('_workflow')
@pytest.mark.parametrize('source', _STAGE_SOURCES, ids=lambda p: p.name)
def test_stages_only_read_from_stages_they_require(source):
    # Splitting the DAG out of the implementations means two files have to agree, and nothing
    # makes them. cpg_flow does catch the drift, but only at graph-build time on a real run,
    # and only on a branch that executes — a read inside a conditional would go unnoticed on
    # any run that skipped it. Reading the calls out of the source catches it at collection
    # time instead.
    wired = _wired_stages()
    for node in ast.walk(ast.parse(source.read_text())):
        if not isinstance(node, ast.ClassDef) or node.name not in wired:
            continue
        required = {dep.__name__ for dep in wired[node.name]().required_stages_classes}
        undeclared = _stages_read_by(node) - required
        assert not undeclared, (
            f'{node.name} reads from {sorted(undeclared)} but does not require them. '
            f'Add them to the requires= for this stage in stages/pipeline.py.'
        )


def test_requested_stages_has_no_duplicates():
    names = [s.__name__ for s in pipeline.REQUESTED_STAGES]
    assert len(names) == len(set(names))


def test_the_entry_point_runs_exactly_the_requested_stages(mocker):
    # An entry point holding its own copy of this list would mean adding a stage to pipeline.py
    # was not enough to make it run — and run_workflow.py is the file a newcomer reading
    # pipeline.py for the pattern never opens. Drive cli_main rather than checking the import,
    # so that filtering the list on the way past would fail too.
    run = mocker.patch('cpg_flow.workflow.run_workflow')
    mocker.patch('sys.argv', ['run_workflow'])
    mocker.patch.object(run_workflow.logging_setup, 'setup_logging')

    run_workflow.cli_main()

    assert run.call_args.kwargs['stages'] is pipeline.REQUESTED_STAGES
    assert run.call_args.kwargs['name'] == 'popgen_rbceq2'


# --- the stage -> config-section convention ---

# Keyed on the real objects, so renaming a class breaks this import instead of leaving the
# section stranded.
_EXPECTED = {
    pipeline.FilterAndConvertGvcfsForRbceq2: 'filter_and_convert_gvcfs_for_rbceq2',
    pipeline.GenotypeBloodGroupsWithRbceq2: 'genotype_blood_groups_with_rbceq2',
    pipeline.FlagBloodGroupCallQc: 'flag_blood_group_call_qc',
    pipeline.CombineRbceq2OutputsPerCohort: 'combine_rbceq2_outputs_per_cohort',
}


def _default_config_stage_sections() -> set[str]:
    """The [workflow.<section>] tables in the default config that a stage should be reading."""
    toml_path = importlib.resources.files('popgen_rbceq2') / 'config' / 'popgen_rbceq2_default_config.toml'
    with toml_path.open('rb') as fh:
        workflow = tomllib.load(fh).get('workflow', {})
    return {key for key, value in workflow.items() if isinstance(value, dict)}


def _fake_config_retrieve(values: Mapping[tuple, object]) -> Callable[..., object]:
    """Stand-in for config_retrieve: returns values[key] if set, else the caller's default."""

    def _cr(keys: list[str], default: object = None) -> object:
        return values.get(tuple(keys), default)

    return _cr


class _FakeJob:
    """Records the image/resource calls stage_support.configure_job makes."""

    def __init__(self) -> None:
        self.set: dict[str, object] = {}

    def image(self, value) -> '_FakeJob':
        self.set['image'] = value
        return self

    def cpu(self, value) -> '_FakeJob':
        self.set['cpu'] = value
        return self

    def memory(self, value) -> '_FakeJob':
        self.set['memory'] = value
        return self

    def storage(self, value) -> '_FakeJob':
        self.set['storage'] = value
        return self


def _configure(monkeypatch, values: Mapping[tuple, object], stage_name: str, **kwargs: Any) -> dict[str, object]:
    """Run configure_job for a stage named stage_name against a fake config; return what it set."""
    monkeypatch.setattr(cpg_utils.config, 'config_retrieve', _fake_config_retrieve(values))
    job = _FakeJob()
    stage_support.configure_job(cast('Any', job), _as_stage(stage_name), **kwargs)
    return job.set


@pytest.mark.parametrize(('stage_obj', 'section'), _EXPECTED.items(), ids=lambda v: getattr(v, '__name__', v))
def test_config_section_matches_expected(stage_obj, section):
    # stage_support.config_section reads stage.name (== class name via the @stage decorator).
    assert stage_support.config_section(_as_stage(stage_obj.__name__)) == section


@pytest.mark.parametrize('section', sorted(_default_config_stage_sections()))
def test_default_config_section_is_read_by_a_stage(section):
    # A renamed stage leaves its tuned cpu/memory/storage stranded in a section nothing reads.
    assert section in set(_EXPECTED.values())


def test_every_stage_has_a_default_config_section():
    # The other direction: a stage with no section runs on the in-code defaults, which are not
    # where anyone looks to tune it.
    assert set(_EXPECTED.values()) <= _default_config_stage_sections()


# --- configure_job: the per-stage [workflow.<section>] override seam ---


def test_configure_job_uses_in_code_defaults_when_section_absent(monkeypatch):
    set_on_job = _configure(
        monkeypatch,
        {},
        'FlagBloodGroupCallQc',
        cpu=1,
        memory='standard',
        storage='10Gi',
        image='driver',
    )
    assert set_on_job == {'image': 'driver', 'cpu': 1, 'memory': 'standard', 'storage': '10Gi'}


def test_configure_job_section_overrides_defaults(monkeypatch):
    # A partial section must override only the keys it sets.
    values = {
        ('workflow', 'filter_and_convert_gvcfs_for_rbceq2', 'cpu'): 16,
        ('workflow', 'filter_and_convert_gvcfs_for_rbceq2', 'storage'): '500Gi',
    }
    set_on_job = _configure(
        monkeypatch,
        values,
        'FilterAndConvertGvcfsForRbceq2',
        cpu=4,
        memory='highmem',
        storage='40Gi',
        image='bcftools',
    )
    assert set_on_job == {'image': 'bcftools', 'cpu': 16, 'memory': 'highmem', 'storage': '500Gi'}


def test_configure_job_reads_its_own_section_not_a_sibling(monkeypatch):
    # A stage must read only its own section, never a neighbouring stage's.
    values = {('workflow', 'flag_blood_group_call_qc', 'cpu'): 16}
    set_on_job = _configure(
        monkeypatch,
        values,
        'FilterAndConvertGvcfsForRbceq2',
        cpu=4,
        memory='highmem',
        storage='40Gi',
        image='bcftools',
    )
    assert set_on_job['cpu'] == 4


def test_configure_job_falls_back_to_driver_image(monkeypatch):
    # image=None is the in-repo case; tool-image stages pass image_path(...) instead.
    values = {('workflow', 'driver_image'): 'driver:latest'}
    set_on_job = _configure(
        monkeypatch,
        values,
        'CombineRbceq2OutputsPerCohort',
        cpu=2,
        memory='standard',
        storage='10Gi',
    )
    assert set_on_job['image'] == 'driver:latest'


# --- resolving what the package ships ---


def test_job_script_resolves_a_committed_script():
    assert stage_support.job_script('rbceq2_gather_job.py').endswith('jobs/rbceq2_gather_job.py')


def test_job_script_rejects_a_script_that_is_not_there():
    # The name is a string, so nothing but this check stands between a typo and a batch
    # submission whose every job dies on a missing file.
    with pytest.raises(FileNotFoundError, match=re.escape('rbceq2_gather.py')):
        stage_support.job_script('rbceq2_gather.py')  # real file ends _job.py


def test_job_script_error_names_what_is_available():
    with pytest.raises(FileNotFoundError, match=re.escape('rbceq2_call_qc_job.py')):
        stage_support.job_script('nope.py')


def test_blood_group_resource_raises_for_a_build_we_have_not_generated():
    with pytest.raises(FileNotFoundError, match=re.escape('gen_bg_resources.py')):
        stage_support.blood_group_resource('bg_regions.GRCh37.bed')


# --- build_python_command: argument formatting for every value type we pass ---

_JOB = 'rbceq2_gather_job.py'


def test_build_python_command_prefixes_scalar_args():
    cmd = stage_support.build_python_command(_JOB, {'output-geno': '/data/x.tsv', 'n': 5})
    assert cmd.startswith(f'python3 {stage_support.job_script(_JOB)}')
    assert '--output-geno /data/x.tsv' in cmd
    assert '--n 5' in cmd  # non-string scalars are stringified


def test_build_python_command_rejects_a_script_that_is_not_there():
    with pytest.raises(FileNotFoundError):
        stage_support.build_python_command('not_a_job.py', {'x': 1})


@pytest.mark.parametrize('value', ['-Xms8g -Xmx8g', '-Xmx8g', '--already-a-flag'])
def test_build_python_command_attaches_dash_leading_values_with_equals(value):
    # argparse reads a leading-dash value as the next option and dies with 'expected one
    # argument'. It excuses strings containing a space; --flag=value holds for both.
    cmd = stage_support.build_python_command(_JOB, {'opt': value})
    assert f'--opt={shlex.quote(value)}' in cmd


def test_build_python_command_dash_leading_value_round_trips_through_argparse():
    parser = argparse.ArgumentParser()
    parser.add_argument('--opt', required=True)
    cmd = stage_support.build_python_command(_JOB, {'opt': '-Xmx8g'})
    argv = shlex.split(cmd.replace('\\\n', ''))[2:]  # drop 'python3 <script>'
    assert parser.parse_args(argv).opt == '-Xmx8g'


def test_build_python_command_quotes_values_with_spaces():
    cmd = stage_support.build_python_command(_JOB, {'path': '/my dir/x.tsv'})
    assert "--path '/my dir/x.tsv'" in cmd


def test_build_python_command_list_emits_repeated_flags():
    # click multiple=True consumes repeated flags, NOT space-joined values.
    cmd = stage_support.build_python_command(_JOB, {'geno': ['a.tsv', 'b.tsv', 'c.tsv']})
    assert cmd.count('--geno ') == 3
    for v in ('a.tsv', 'b.tsv', 'c.tsv'):
        assert f'--geno {v}' in cmd


def test_build_python_command_rejects_an_empty_list():
    # An empty list would drop the flag entirely, and the job would fail on a missing required
    # option far from the stage that built the command.
    with pytest.raises(ValueError, match='empty list'):
        stage_support.build_python_command(_JOB, {'geno': []})


def test_build_python_command_bool_true_is_bare_flag():
    cmd = stage_support.build_python_command(_JOB, {'output-pdfs': True})
    assert cmd.rstrip().endswith('--output-pdfs')  # is_flag: flag present, no value


def test_build_python_command_bool_false_and_none_are_omitted():
    cmd = stage_support.build_python_command(_JOB, {'output-pdfs': False, 'optional': None, 'keep': 'yes'})
    assert 'output-pdfs' not in cmd
    assert 'optional' not in cmd
    assert '--keep yes' in cmd


def test_build_python_command_preserves_arg_order():
    cmd = stage_support.build_python_command(_JOB, {'a': 1, 'b': 2, 'c': 3})
    assert cmd.index('--a') < cmd.index('--b') < cmd.index('--c')
