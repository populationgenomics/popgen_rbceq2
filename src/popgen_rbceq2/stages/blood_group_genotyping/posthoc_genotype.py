"""Re-genotypes blood-group defining sites the exome capture left uncalled, from the CRAM."""

import typing

import cpg_flow.stage
import cpg_flow.targets
import cpg_utils.config
import cpg_utils.hail_batch
import hailtop.batch.resource

from popgen_rbceq2 import constants, stage_support

# The sequencing type this stage exists for. An exome gVCF is called against a capture-target
# BED, which is what puts defining sites outside it beyond reach; a genome gVCF has no such
# edge, so recalling from its CRAM would re-derive calls DRAGEN already made.
EXOME = 'exome'

# The memory tier this stage asks for when its config section does not set one.
MEMORY_TIER = 'standard'

# Container memory per cpu, in GB, for each Hail Batch memory tier this stage accepts, rounded
# down from what the tier really grants (standard is 3.75, highmem 6.5).
#
# The JVM heap is derived from this and not from `cpu` alone. `configure_job` takes the tier
# from this stage's own config section, so tuning memory there without also raising cpu would
# leave -Xmx asking for more than the container has. The job is then killed for exceeding its
# limit, with nothing in the log naming the heap as the cause.
#
# `lowmem` is not accepted. Hail grants it about 0.9 GiB per core, which rounds down to nothing
# a JVM can be sized from, and HaplotypeCaller is not a lowmem workload; a request for it fails
# at graph-build time with the message below rather than being killed mid-run.
_GB_PER_CPU = {'standard': 3, 'highmem': 6}

# Held back from the heap for the rest of the JVM: thread stacks, metaspace, GC structures and
# the direct buffers htsjdk uses for BGZF. One core on `standard` leaves 2GB after this, the
# smallest heap any accepted pairing can produce, and enough to start HaplotypeCaller on.
_JVM_OVERHEAD_GB = 1


def _heap_gb(cpu: int, memory: str, section: str) -> int:
    """The JVM heap, in GB, that fits the container this stage's cpu and memory tier give it.

    Args:
        cpu: Cores the job requests.
        memory: Hail Batch memory tier the job requests.
        section: This stage's config section, for the error messages.

    Returns:
        The heap size to pass as -Xmx.

    Raises:
        cpg_utils.config.ConfigError: `memory` is not an accepted tier.
    """
    if memory not in _GB_PER_CPU:
        raise cpg_utils.config.ConfigError(
            f'workflow.{section}.memory is {memory!r}; this stage sizes its JVM heap from the '
            f'tier it was given, so it needs one of {sorted(_GB_PER_CPU)}. An explicit size is '
            'not supported here, because nothing would then keep the heap inside it, and lowmem '
            'is not, because its ~0.9 GiB per core leaves no heap worth starting HaplotypeCaller on.'
        )
    return cpu * _GB_PER_CPU[memory] - _JVM_OVERHEAD_GB


def applies_to(sequencing_group: cpg_flow.targets.SequencingGroup) -> bool:
    """Whether post-hoc calling runs for this sequencing group.

    Both this stage's `expected_outputs` and the merge in `FilterAndConvertGvcfsForRbceq2`
    gate on this one predicate. If they could disagree, the conversion stage would ask
    cpg_flow for an output this stage never produced, and fail at graph-build time on a
    subset of a cohort.

    Args:
        sequencing_group: The sequencing group to test.

    Returns:
        True for an exome with both a CRAM to call from and a gVCF to supplement.
    """
    return sequencing_group.sequencing_type == EXOME and bool(sequencing_group.cram) and bool(sequencing_group.gvcf)


class PosthocGenotypeOffTargetSites(cpg_flow.stage.SequencingGroupStage):
    """Call blood-group defining sites from an exome CRAM, to fill the capture's blind spots.

    DRAGEN calls an exome against the capture-target BED with no `--vc-target-bed-padding`, so
    gVCF emission hard-stops on the BED edges. A defining coordinate outside the capture gets
    no record at all — not a low-quality one, none — and rbceq2 reads a site missing from its
    input as a *confident homozygous reference call*. Our QC flags that `NOCOV`, so it is at
    least visible, but the system stays unassessable.

    The reads are usually there. Measured over 20 DRAGEN 3.7.8 exomes (5 Twist VCGS, 15
    Agilent CREv2) in 2026-08: ~165 of 1,599 non-HPA defining coordinates are off-target per
    capture design, and ~105 of those sit 1-100bp from a target edge carrying 40-110x MAPQ>=20
    depth in the CRAM, the FY GATA Duffy-null promoter sites among them at 34-210x. Only the
    caller stopped early. Re-calling those sites lifts assessable non-HPA coordinates from
    91.7% to ~98% per design; re-running DRAGEN over a cohort is not affordable.

    This stage calls every assessable defining site rather than only the off-target ones,
    because which sites are off-target is a property of the sample's capture kit and this
    stage has no reason to know it. Nothing is gained by knowing: the whole padded interval
    list is 199 regions over 136kb, so calling all of it costs the same as calling part.
    Deciding what to *keep* is `FilterAndConvertGvcfsForRbceq2`'s job. It keeps a post-hoc
    record only where the DRAGEN gVCF is silent *and* the site lies outside the cohort's
    capture design, which that stage reads from config. Silence is judged per sample from the
    gVCF itself, so a design BED that misdescribes the real footprint still cannot cause a
    DRAGEN call to be overwritten; the design only ever narrows what may be filled.

    Three things to preserve when changing this stage:

    - **The CRAM is streamed, not localised.** The `gs://` path goes to `-I` and GATK reads it
      over NIO, pulling only the blocks `-L` asks for. The interval list is 136kb of genome,
      so this reads megabytes where localising an exome CRAM would move gigabytes per sample.
    - **The reference must be the one the CRAM was aligned against**, the DRAGEN masked
      assembly38 (`references.broad.ref_fasta`). CRAM stores reads differentially against the
      reference, so decoding with a different assembly does not fail loudly, it yields wrong
      bases.
    - **`--dragen-mode`, and no re-alignment.** The primary calls this supplements come from
      DRAGEN 3.7.8, so the supplement is made as close to DRAGEN-equivalent as a re-call can
      be. Reads are used as aligned: no DRAGMAP, no re-alignment. There is no DRAGstr model
      for the sample, so STR genotyping is not fully DRAGEN-equivalent; that is accepted for
      now and worth revisiting if post-hoc indel calls near STRs look discordant.

    The output is an intermediate consumed through the cpg-flow graph, so it is written under
    the tmp category and registers no Analysis of its own. Provenance reaches Metamist through
    the QC flags instead: a site filled from here is flagged `POSTHOC` and never a silent
    `PASS`, and `FlagBloodGroupCallQc` records the caller version in its Analysis meta.
    """

    def expected_outputs(
        self, sequencing_group: cpg_flow.targets.SequencingGroup
    ) -> stage_support.ExpectedOutputs | None:
        # Skipped, not failed, on all three counts, matching how the rest of the pipeline
        # treats a sequencing group it cannot process: a genome (nothing to supplement), no
        # CRAM (nothing to call from), or no gVCF (no primary calls to supplement, so the
        # conversion stage produces nothing for this sequencing group either).
        if not applies_to(sequencing_group):
            return None
        prefix = stage_support.get_sg_output_prefix(sequencing_group, stage_name=self.name, category='tmp')
        return {
            'gvcf': prefix / f'{sequencing_group.id}.posthoc.g.vcf.gz',
            'index': prefix / f'{sequencing_group.id}.posthoc.g.vcf.gz.tbi',
        }

    def queue_jobs(
        self,
        sequencing_group: cpg_flow.targets.SequencingGroup,
        inputs: cpg_flow.stage.StageInput,  # noqa: ARG002
    ) -> cpg_flow.stage.StageOutput | None:
        outputs = self.expected_outputs(sequencing_group)
        if outputs is None:
            return None
        cfg = stage_support.config_section(self)
        cpu = cpg_utils.config.config_retrieve(['workflow', cfg, 'cpu'], 2)
        # Read here rather than left to configure_job's own fallback, so the heap below and the
        # container are sized from one value.
        memory = cpg_utils.config.config_retrieve(['workflow', cfg, 'memory'], MEMORY_TIER)
        genome = cpg_utils.config.genome_build()

        b = cpg_utils.hail_batch.get_batch()
        j = b.new_bash_job(
            f'PosthocGenotypeOffTargetSites/{sequencing_group.id}',
            self.get_job_attrs(sequencing_group) | {'tool': 'gatk'},
        )
        j = stage_support.configure_job(
            j,
            self,
            cpu=cpu,
            memory=memory,
            storage='20Gi',
            image=cpg_utils.config.image_path('gatk', constants.GATK_IMAGE_TAG),
        )

        # The fasta is localised (GATK reads it randomly throughout the run, and needs the .fai
        # and .dict beside it); the CRAM is not (GATK reads a few hundred kb of it, by index).
        fasta_path = cpg_utils.config.reference_path('broad/ref_fasta')
        reference = b.read_input_group(
            base=str(fasta_path),
            fai=f'{fasta_path}.fai',
            dict=str(fasta_path).removesuffix('.fasta') + '.dict',
        )
        padded_bed = b.read_input(stage_support.blood_group_resource(f'bg_defining_sites_padded.{genome}.bed'))
        j.declare_resource_group(out={'g.vcf.gz': '{root}.g.vcf.gz', 'g.vcf.gz.tbi': '{root}.g.vcf.gz.tbi'})
        out = typing.cast('hailtop.batch.resource.ResourceGroup', j.out)

        # -ERC GVCF so the output carries reference blocks, not just variants. A defining site
        # the sample is genuinely hom-ref at has to arrive as a block with a DP and a GQ, or
        # filling the hole would only replace "no record" with "no record", and the QC could
        # not tell a confident reference call from an absent one.
        #
        # The CRAM index is not passed separately: GATK derives the .crai path from the .cram
        # path over NIO the same way it does on a local file.
        j.command(
            f"""
            set -euxo pipefail
            gatk --java-options "-Xms1g -Xmx{_heap_gb(cpu, memory, cfg)}g" HaplotypeCaller \\
                -R {reference.base} \\
                -I {sequencing_group.cram!s} \\
                -L {padded_bed} \\
                -O {out['g.vcf.gz']} \\
                -ERC GVCF \\
                --dragen-mode \\
                --create-output-variant-index true
            if [ ! -s {out['g.vcf.gz.tbi']} ]; then
                echo "ERROR: HaplotypeCaller wrote no index beside its GVCF." >&2
                exit 1
            fi
            """,
        )
        b.write_output(out, str(outputs['gvcf']).removesuffix('.g.vcf.gz'))
        return self.make_outputs(sequencing_group, data=outputs, jobs=[j])
