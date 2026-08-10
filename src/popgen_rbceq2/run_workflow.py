#!/usr/bin/env python3

"""Entry point for the workflow.

It holds no stage list of its own. Which stages run is declared by REQUESTED_STAGES in
stages/pipeline.py, alongside the DAG those stages form; this module hands that list to
cpg_flow and lets it discover the graph.
"""

import argparse

import cpg_flow.workflow

from popgen_rbceq2 import logging_setup
from popgen_rbceq2.stages import pipeline

WORKFLOW_NAME = 'popgen_rbceq2'


def cli_main() -> None:
    """CLI entrypoint — starts up the workflow."""
    logging_setup.setup_logging(force=True)

    parser = argparse.ArgumentParser()
    parser.add_argument('--dry_run', action='store_true', help='Dry run')
    args = parser.parse_args()

    cpg_flow.workflow.run_workflow(
        name=WORKFLOW_NAME,
        stages=pipeline.REQUESTED_STAGES,
        dry_run=args.dry_run,
    )


if __name__ == '__main__':
    cli_main()
