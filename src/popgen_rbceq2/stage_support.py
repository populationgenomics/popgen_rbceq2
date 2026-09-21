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
    """Apply the stage's own meta function, then add what a meta function cannot know.

    A meta function is handed the output path and nothing else, so five keys are added here.
    Two are the stage class:

    - ``stage``, the class name. cpg_flow's ``get_job_attrs`` puts ``stage=<class name>`` on
      the Hail *job*, not in ``Analysis.meta``, so this hook is the only thing that puts it in
      the meta — and adding it last also stops a meta function carrying a stale hand-typed
      one, which the pre-wire house style did.
    - ``stage_version``, the release the stage's outputs were written under, from
      `_release_version` — a config lookup keyed on the class name, so a stage with its own
      ``output_versions`` pin records the release it actually wrote into. It is the writing
      stage's own release, not its inputs'; see `_release_version`.

    Three are properties of the run, the same for every stage in it:

    - ``rbceq2_version``, the tool, from `constants`.
    - ``rbceq2_db_version``, the allele database, from `constants`. It moves independently of
      the tool: the committed ``resources/bg_*`` are built from the db, so a db bump changes
      the QC flags with no tool bump to show for it. Deliberately *not* an axis of the output
      tree — adding it to `_release_tree` would orphan every tree written so far — so a row's
      meta is the only place a db bump shows.
    - ``exome_design``, the capture design an exome run was called against, as configured, or
      None for a genome. Recorded on every row rather than only on exome rows so the meta
      shape does not depend on sequencing type (as call_qc's static meta already does).

    All but the db version are axes of the output tree (see `_release_tree`), and together the
    five are what make a row self-describing: Metamist inserts an Analysis per run and retires
    none, so a reader meeting two rows for one sequencing group has only the meta to tell which
    release, tool, database and design produced each.

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
    return extra | {
        'stage': stage_name,
        'stage_version': _release_version(stage_name),
        'rbceq2_version': constants.RBCEQ2_VERSION,
        'rbceq2_db_version': constants.RBCEQ2_DB_VERSION,
        'exome_design': exome_design_bed() if _is_exome_run() else None,
    }


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
            returning extra Analysis.meta. The stage name, its release version, the rbceq2
            tool and db versions and the exome design are added for you — see
            `_with_stage_name`.
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

    The recorded meta is ``analysis_meta.call_qc``'s dict plus the keys `_with_stage_name`
    adds, among them ``{'stage': 'FlagBloodGroupCallQc', 'stage_version': <release>}``:
    ``stage`` is the class name itself and ``stage_version`` a config lookup keyed on it, so
    renaming the stage moves both without anyone editing a string literal.

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


def blood_group_resource(name: str, missing_hint: str | None = None) -> str:
    """Resolve a committed blood-group site resource to a path.

    Args:
        name: Resource filename, e.g. `bg_regions.GRCh38.bed`.
        missing_hint: What the error should tell the reader to do if the file is not shipped.
            Defaults to regenerating the `bg_*.<genome>.*` set from the rbceq2 database.

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
            missing_hint
            or (
                'Generate it with scripts/gen_bg_resources.py against the db.tsv from the pinned '
                'rbceq2 image, and commit it under resources/.'
            ),
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


# The `[workflow]` key naming the capture design an exome cohort was called against, spelled as
# the design's key in the references repo, e.g.
# `exome_probesets_hg38/agilent_sureselect_clinical_research_exome_v2_covered_by_probes_bed`.
# Required for an exome run; a genome run never reads it. The pipeline does not resolve it to
# the vendor BED: it selects the committed subtraction of that BED from the defining sites
# (`off_design.resource_path`), and is a segment of every exome output path (`_release_tree`).
EXOME_DESIGN_KEY = 'exome_design_bed'
# The fully-qualified path, for the error messages that have to name it.
DESIGN_CONFIG_PATH = f'workflow.{EXOME_DESIGN_KEY}'


def exome_design_bed() -> str:
    """The key of the capture design an exome run fills holes outside of.

    Read at graph-build time by `_release_tree`, so a run missing it fails before a job starts
    rather than on the first exome sequencing group's merge. Also read by `_with_stage_name`,
    which runs in the status-reporter job after the compute — it reads the same config snapshot
    the driver did, so the graph-build check is what keeps the error below from surfacing
    there. Whether a subtraction is committed for the key is checked where it is read, in
    `off_design.resource_path`.

    Returns:
        The design key as configured.

    Raises:
        cpg_utils.config.ConfigError: The key is not set.
    """
    try:
        return cpg_utils.config.config_retrieve(['workflow', EXOME_DESIGN_KEY])
    except cpg_utils.config.ConfigError as e:
        raise cpg_utils.config.ConfigError(
            f'An exome run needs {DESIGN_CONFIG_PATH}: the references-repo key of the capture design '
            'BED the gVCFs were called against, e.g. '
            "'exome_probesets_hg38/twist_vcgs_custom_exome_covered_targets_bed'. Post-hoc calls "
            'fill defining sites only outside that design, and the sites outside each known design '
            'are committed under resources/.'
        ) from e


def _is_exome_run() -> bool:
    """Whether this run's configured sequencing type is exome.

    Read by `_release_tree` for the design segment of the path and by `_with_stage_name` for
    the design in the Analysis meta, so the two cannot disagree about which runs have one.
    """
    return cpg_utils.config.config_retrieve(['workflow', 'sequencing_type']) == constants.EXOME


def design_segment(design_key: str) -> str:
    """The output-path segment identifying an exome run's capture design.

    Args:
        design_key: The `[references]` key the design came from.

    Returns:
        The key with anything outside `[A-Za-z0-9._-]` replaced, so it is one path segment.
    """
    return re.sub(r'[^A-Za-z0-9._-]', '_', design_key)


def _release_version(stage_name: str) -> str:
    """The release a stage writes under: its own output_versions pin if set, else workflow.version.

    Read by `_release_tree` for the output path and by `_with_stage_name` for the Analysis
    meta, so the recorded release cannot disagree with the tree written into.

    It answers for the named stage only. A stage consuming another's outputs reads whatever
    tree that stage's own pin selected, so pinning one stage mid-release leaves every
    downstream stage writing and recording its own unpinned release over inputs from a pinned
    one. Pin the downstream stages too: merely forcing them rebuilds at the same path and
    records the same release, so the new Analysis cannot be told from the one it supersedes.

    A pin is honoured as set, not for being truthy: an empty or blank one is a mistake, and
    silently falling back to ``workflow.version`` would write the pinned stage into the shared
    tree and record the shared release, which is the opposite of what pinning was for. The
    value is stringified so a TOML integer cannot type ``meta.stage_version`` differently from
    one row to the next.

    Returns:
        The release, always a non-empty string.

    Raises:
        cpg_utils.config.ConfigError: Neither a pin nor ``workflow.version`` is set, or the
            one that is set is blank. Not defaulted: a release this run never declared would
            name an output tree and assert itself in Metamist on every row.
    """
    pinned = cpg_utils.config.config_retrieve(['workflow', 'output_versions', stage_name], None)
    source = 'workflow.version' if pinned is None else f'workflow.output_versions.{stage_name}'
    raw = cpg_utils.config.config_retrieve(['workflow', 'version'], None) if pinned is None else pinned
    version = str(raw).strip() if raw is not None else ''
    if not version:
        raise cpg_utils.config.ConfigError(
            f"{source} names the release, e.g. 'v4'; got {raw!r}. It is the segment of every "
            'output path after the rbceq2 version, and what every Analysis records as '
            'meta.stage_version. The shipped default config sets workflow.version, so a run '
            'reaching this either does not merge the defaults or pins a stage to a blank value.'
        )
    return version


def _release_tree(stage_name: str) -> str:
    """The release tree a stage's outputs land in, as the path below the workflow name.

    `rbceq2_<tool version>_<release>` for a genome run, and
    `rbceq2_<tool version>_<release>/<design>` for an exome run.

    The tool-version half is derived from constants.RBCEQ2_VERSION, so a tool bump always
    lands in a fresh tree and the segment can never drift from the version actually run. The
    release half is ours — the stage's own output_versions pin if set, else workflow.version —
    bumped only when a pipeline change alters the outputs (deliberately not the image tag,
    which moves on rebuilds that change nothing about the outputs).

    The design segment is there because cpg_flow reuses a stage whose expected outputs exist
    and asks nothing about how they were made. Every exome output from the conversion onward
    depends on the configured design: which holes the merge filled, so the converted VCF, the
    genotypes rbceq2 calls from it, the QC flags and the cohort tables. None of those stages
    reads the design itself, so none could put it in its own path, and the release segment
    cannot see a config change. Repointing EXOME_DESIGN_KEY inside one release would otherwise
    select the other design's committed off-design sites and then reuse every sample's outputs
    from the old design without a word in any log. With the design in the tree, repointing it
    starts a fresh tree for the whole run. A genome run never reads the key, so its tree is
    unchanged.
    """
    tree = f'rbceq2_{constants.RBCEQ2_VERSION.replace(".", "_")}_{_release_version(stage_name)}'
    if _is_exome_run():
        tree = f'{tree}/{design_segment(exome_design_bed())}'
    return tree


def get_output_prefix(cohort: cpg_flow.targets.Cohort, stage_name: str, category: str | None = None) -> cpg_utils.Path:
    """Standardised output prefix for CohortStage outputs.

    Format: cohort.dataset.prefix() / workflow.name / <release tree> / stage_name / cohort.id

    The release tree sits directly under the workflow name so one release is one browsable
    tree, and for an exome run one design within it; see _release_tree for its segments. A
    stage with its own output_versions pin writes under its pinned release's tree instead.

    cohort.id is a path segment, so a different set of sequencing groups is a different cohort
    and therefore a different tree — outputs from one cohort can never be mistaken for another's.
    """
    return (
        cohort.dataset.prefix(category=category)
        / cpg_flow.workflow.get_workflow().name
        / _release_tree(stage_name)
        / stage_name
        / cohort.id
    )


def get_sg_output_prefix(
    sequencing_group: cpg_flow.targets.SequencingGroup,
    stage_name: str,
    category: str | None = None,
) -> cpg_utils.Path:
    """Standardised output prefix for SequencingGroupStage outputs.

    Format: sg.dataset.prefix() / workflow.name / <release tree> / stage_name / sg.id

    See get_output_prefix and _release_tree for what the release tree's segments mean.
    """
    return (
        sequencing_group.dataset.prefix(category=category)
        / cpg_flow.workflow.get_workflow().name
        / _release_tree(stage_name)
        / stage_name
        / sequencing_group.id
    )


def camel_to_snake(name: str) -> str:
    s = re.sub(r'([A-Z]+)([A-Z][a-z])', r'\1_\2', name)
    s = re.sub(r'([a-z0-9])([A-Z])', r'\1_\2', s)
    return s.lower()


def config_section(stage: cpg_flow.stage.Stage) -> str:
    """The [workflow.<section>] this stage reads, derived from its class name."""
    return camel_to_snake(stage.name)


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
