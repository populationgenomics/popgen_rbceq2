"""Shared helpers for stage orchestration in the stages package.

Convention: a stage's `[workflow.<section>]` config key is derived from its class name by
`config_section` (camelCase -> snake_case), so the class name and its config section must stay
in sync. Renaming a stage class silently repoints its config section — rename the TOML section
to match (see test_stage_support, which keeps the mapping under test).

The `jobs/<name>_job.py` a stage runs is *not* derived from the class name and does not have to
match it: FlagBloodGroupCallQc runs `rbceq2_call_qc_job.py`. The filename is passed explicitly,
and `build_python_command` resolves it through `job_script`, which fails on a name that is not
there rather than building a command that dies once the job starts.
"""

import functools
import importlib.resources
import inspect
import re
import shlex
from collections.abc import Callable, Mapping
from importlib.resources.abc import Traversable
from typing import Any, TypeAlias

import cpg_flow.stage
import cpg_flow.targets
import cpg_flow.workflow
import cpg_utils
import cpg_utils.config
import hailtop.batch.job
import hailtop.batch.resource

from popgen_rbceq2 import constants

# Values accepted by build_python_command: paths, localised batch resources, and scalars. None
# is a member because an optional argument is passed as None and omitted from the command, so a
# stage can hand over a config value it did not find without branching around the call.
JobArg: TypeAlias = str | cpg_utils.Path | hailtop.batch.resource.Resource | bool | int | list[str] | None

# What every stage here returns from expected_outputs: the dict member of cpg_flow's
# ExpectedResultT. Declaring the dict rather than the whole union is what lets `outputs[key]`
# type-check — the union also holds a bare Path and a str, which are not subscriptable by name.
# The value type is cpg_flow's own, left unnarrowed so the override stays compatible with the
# base signature (dict is invariant in its value type).
ExpectedOutputs: TypeAlias = dict[str, str | cpg_utils.Path | list[str | cpg_utils.Path]]


def _with_stage_name(output: str, *, fn: 'Callable[[str], dict] | None', stage_name: str) -> dict:
    """Apply the stage's own meta function, then re-assert the stage name over its result.

    cpg_flow already records ``stage=<class name>`` in every Analysis meta (via
    ``get_job_attrs``), but the status reporter merges the meta function's dict in last, so a
    function that set ``stage`` itself would override the framework's value with a hand-typed
    one. Re-applying the class name after the function closes that hole; deriving anything
    here is not the point.

    A partial over this is what wire hands to cpg_flow. It is a module-level function rather
    than a closure because cpg_flow ships the callable into a Hail PythonJob, where it is
    dill-pickled; a partial over a module-level function pickles by reference.

    Raises:
        TypeError: If the stage's function returns anything but a dict. Left to fail here
            rather than coerced, but note this runs in the status-reporter job, long after the
            compute — which is why wire checks what it can at import time instead.
    """
    extra = fn(output) if fn else {}
    if not isinstance(extra, dict):
        raise TypeError(
            f'update_analysis_meta for {stage_name} returned {type(extra).__name__}, expected dict',
        )
    return extra | {'stage': stage_name}


def wire(
    cls: type[cpg_flow.stage.Stage],
    requires: list[cpg_flow.stage.StageDecorator] | None = None,
    *,
    analysis_type: str | None = None,
    analysis_keys: list[cpg_utils.Path | str] | None = None,
    update_analysis_meta: 'Callable[[str], dict] | None' = None,
    **stage_kwargs: Any,
) -> cpg_flow.stage.StageDecorator:
    """Attach a stage implementation to the DAG and declare its Metamist registration.

    Args:
        cls: The stage implementation class, undecorated.
        requires: Stages this one consumes output from. Empty for an entry point.
        analysis_type: Metamist analysis type. Omit to record no Analysis.
        analysis_keys: Which expected_outputs keys to register. Required when analysis_type
            is set and expected_outputs returns a dict.
        update_analysis_meta: Module-level function (not a method) taking the output path and
            returning extra Analysis.meta. The stage name is added for you.
        **stage_kwargs: Passed through to cpg_flow's ``@stage`` — e.g. tolerate_missing_output.

    Returns:
        The class decorated with cpg_flow's ``@stage``.

    Raises:
        TypeError: If update_analysis_meta does not take exactly one argument.
        ValueError: If analysis_keys is given without analysis_type, which would silently
            record nothing.

    Example — registering the per-SG QC TSV:

        FlagBloodGroupCallQc = wire(
            call_qc.FlagBloodGroupCallQc,
            requires=[FilterAndConvertGvcfsForRbceq2, GenotypeBloodGroupsWithRbceq2],
            analysis_type='blood_group_qc',
            analysis_keys=['qc'],
            update_analysis_meta=analysis_meta.call_qc,
        )

    The recorded meta is ``analysis_meta.call_qc``'s dict plus
    ``{'stage': 'FlagBloodGroupCallQc'}`` — the name comes from the class, so renaming the
    stage moves it without anyone editing a string literal.

    Note: when cpg_flow reports "getting inputs from stage X, but X is not listed in
    required_stages. Consider adding it into the decorator: @stage(required_stages=[X])", the
    fix here is to add X to ``requires=`` in stages/pipeline.py. The stage classes carry no
    decorator, and adding one would double-decorate them.
    """
    if analysis_keys and not analysis_type:
        raise ValueError(
            f'{cls.__name__}: analysis_keys={analysis_keys} was given without analysis_type, '
            f'so no Analysis would be recorded and the keys would be ignored.',
        )
    if update_analysis_meta is not None:
        try:
            inspect.signature(update_analysis_meta).bind('<output path>')
        except TypeError as e:
            raise TypeError(
                f'{cls.__name__}: update_analysis_meta must be callable with one argument, the '
                f'output path. A method picks up self and fails inside the Metamist status job, '
                f'long after the compute has run — use a module-level function.',
            ) from e
        except ValueError:
            # No introspectable signature (a builtin, or a C function). Nothing to check; a bad
            # return value is still caught by _with_stage_name.
            pass

    meta_hook = None
    if analysis_type:
        meta_hook = functools.partial(_with_stage_name, fn=update_analysis_meta, stage_name=cls.__name__)

    return cpg_flow.stage.stage(
        required_stages=requires or [],
        analysis_type=analysis_type,
        analysis_keys=analysis_keys,
        update_analysis_meta=meta_hook,
        **stage_kwargs,
    )(cls)


def _package_file(subdirectory: str, name: str, missing_hint: str) -> Traversable:
    """Resolve a file shipped inside the installed package.

    Goes through importlib.resources rather than ``__file__``, so the path is the installed
    package's, whatever installed it. Anything reached this way has to be declared as a wheel
    artifact in pyproject.toml, or it resolves in a source checkout and vanishes in the image.

    Raises:
        FileNotFoundError: The file is not shipped. Callers run during graph construction, so
            this fails the submission rather than every job that needed the file.
    """
    directory = importlib.resources.files('popgen_rbceq2').joinpath(subdirectory)
    resource = directory.joinpath(name)
    if not resource.is_file():
        available = sorted(entry.name for entry in directory.iterdir() if entry.is_file())
        raise FileNotFoundError(
            f'{subdirectory}/{name} is not shipped. {missing_hint} {subdirectory}/ holds {available}',
        )
    return resource


def blood_group_resource(name: str) -> str:
    """Resolve a committed blood-group site resource to a path.

    Args:
        name: Resource filename, e.g. `bg_regions.GRCh38.bed`.

    Returns:
        The absolute path to the shipped resource.

    Raises:
        FileNotFoundError: The resource is not shipped, i.e. no resources have been generated
            for the configured reference build.
    """
    return str(
        _package_file(
            'resources',
            name,
            'Generate it with scripts/gen_bg_resources.py against the db.tsv from the pinned '
            'rbceq2 image, and commit it under resources/.',
        ),
    )


def job_script(name: str) -> str:
    """Resolve a jobs/ script to the path it has inside the driver image.

    Args:
        name: Job module filename, e.g. `rbceq2_call_qc_job.py`.

    Returns:
        The absolute path to the shipped script.

    Raises:
        FileNotFoundError: No such script. The name is a string a stage passes to
            build_python_command, so this catches a typo at graph-build time instead of
            letting every job in the stage start and die on a missing file.
    """
    return str(_package_file('jobs', name, 'Stages may only run a script committed to jobs/.'))


def _output_version(stage_name: str) -> str:
    """The version segment for a stage's outputs: `rbceq2_<tool version>_<release>`.

    The tool-version half is derived from constants.RBCEQ2_VERSION, so a tool bump always
    lands in a fresh tree and the segment can never drift from the version actually run. The
    release half is ours — the stage's own output_versions pin if set, else workflow.version —
    bumped only when a pipeline change alters the outputs (deliberately not the image tag,
    which moves on rebuilds that change nothing about the outputs).
    """
    pinned = cpg_utils.config.config_retrieve(['workflow', 'output_versions', stage_name], None)
    release = pinned or cpg_utils.config.config_retrieve(['workflow', 'version'], 'v1')
    return f'rbceq2_{constants.RBCEQ2_VERSION.replace(".", "_")}_{release}'


def get_output_prefix(cohort: cpg_flow.targets.Cohort, stage_name: str, category: str | None = None) -> cpg_utils.Path:
    """Standardised output prefix for CohortStage outputs.

    Format: cohort.dataset.prefix() / workflow.name / rbceq2_<tool version>_<release> / stage_name / cohort.id

    The version segment sits directly under the workflow name so one release is one browsable
    tree; see _output_version for what its two halves mean. A stage with its own
    output_versions pin writes under its pinned release's tree instead.

    cohort.id is a path segment, so a different set of sequencing groups is a different cohort
    and therefore a different tree — outputs from one cohort can never be mistaken for another's.
    """
    return (
        cohort.dataset.prefix(category=category)
        / cpg_flow.workflow.get_workflow().name
        / _output_version(stage_name)
        / stage_name
        / cohort.id
    )


def get_sg_output_prefix(
    sequencing_group: cpg_flow.targets.SequencingGroup,
    stage_name: str,
    category: str | None = None,
) -> cpg_utils.Path:
    """Standardised output prefix for SequencingGroupStage outputs.

    Format: sg.dataset.prefix() / workflow.name / rbceq2_<tool version>_<release> / stage_name / sg.id

    See get_output_prefix and _output_version for what the version segment means.
    """
    return (
        sequencing_group.dataset.prefix(category=category)
        / cpg_flow.workflow.get_workflow().name
        / _output_version(stage_name)
        / stage_name
        / sequencing_group.id
    )


def _camel_to_snake(name: str) -> str:
    s = re.sub(r'([A-Z]+)([A-Z][a-z])', r'\1_\2', name)
    s = re.sub(r'([a-z0-9])([A-Z])', r'\1_\2', s)
    return s.lower()


def config_section(stage: cpg_flow.stage.Stage) -> str:
    """The [workflow.<section>] this stage reads, derived from its class name."""
    return _camel_to_snake(stage.name)


def _resolved(stage: cpg_flow.stage.Stage, key: str, default: Any) -> Any:
    """A stage's compute value: its config section if set, else the in-code default."""
    return cpg_utils.config.config_retrieve(['workflow', config_section(stage), key], default)


def configure_job(
    job: hailtop.batch.job.BashJob,
    stage: cpg_flow.stage.Stage,
    *,
    cpu: int,
    memory: str,
    storage: str,
    image: str | None = None,
) -> hailtop.batch.job.BashJob:
    """Set the image and resources on a stage's job from its config section.

    image is a fully-resolved image string. When omitted (None), the driver image
    is used; tool-image stages pass e.g. image_path('bcftools', <version>).

    Typed to BashJob rather than Job because `image` is declared on BashJob and PythonJob, not
    on their base: `Job.__getattr__` would resolve it to a Resource, and the call would fail at
    graph-construction time instead of being a type error.
    """
    job.image(image if image is not None else cpg_utils.config.config_retrieve(['workflow', 'driver_image']))
    job.cpu(_resolved(stage, 'cpu', cpu))
    job.memory(_resolved(stage, 'memory', memory))
    job.storage(_resolved(stage, 'storage', storage))
    return job


def build_python_command(name: str, args: Mapping[str, JobArg]) -> str:
    """Build a `python3 <jobs/name> --flag value ...` shell command.

    Args:
        name: Job module filename, e.g. 'rbceq2_gather_job.py'. Resolved through job_script,
            so an unknown name fails here rather than inside the job.
        args: {hyphenated-flag: value}; values are stringified and shell-quoted.
            - list/tuple -> the flag repeated once per item (argparse/click multiple=True)
            - True       -> bare flag (a store_true/is_flag option); False/None -> omitted

    Returns:
        The shell command, one flag per line.
    """
    parts = [f'python3 {job_script(name)}']

    for flag, val in args.items():
        if val is None:
            continue
        if isinstance(val, bool):
            if val:
                parts.append(f'  --{flag}')
        elif isinstance(val, (list, tuple)):
            if not val:
                raise ValueError(f'--{flag} was given an empty list; the flag would be dropped silently')
            parts.extend(f'  {_arg(flag, v)}' for v in val)
        else:
            parts.append(f'  {_arg(flag, val)}')
    return ' \\\n'.join(parts)


def _arg(flag: str, value: Any) -> str:
    """One `--flag value` pair, quoted for the shell.

    A value starting with '-' is written as --flag=value. argparse otherwise reads it as the
    next option and fails with 'expected one argument' — it makes an exception only for
    strings containing a space.
    """
    text = str(value)
    separator = '=' if text.startswith('-') else ' '
    return f'--{flag}{separator}{shlex.quote(text)}'
