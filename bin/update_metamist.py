#!/usr/bin/env python3

"""
Generic Metamist analysis generation script, must be executable

- takes a path to a primary output, and optionally one or more secondary analyses
  - due to the fixed dictionary structure (type of each secondary is within schema), this requires 'k=v k2=v2' args.
- requires an analysis type to register this as
- can take either SGs or Cohorts to associate the analysis entry with, but not both
- optionally can use the same k=v breaking to register metadata

e.g. CLI

# minimal example
python3 update_metamist.py \
    --project <project> \
    --output <path/to/primary/output> \
    --type cram \
    --cohorts COH1234 \
    --dry

# with optional metaadta
python3 update_metamist.py \
    --project <project> \
    --output <path/to/primary/output> \
    --type cram \
    --secondary index=path/to/index html=/path/to/html \
    --meta stage=name sequencing_type=genome \
    --cohorts COH1234 COH5678
"""

import argparse
import json
import sys

from collections import defaultdict
from typing import TypeAlias

from metamist import graphql


RecursiveDict: TypeAlias = dict[str, 'str | RecursiveDict']


def create_output_block(primary: str, type: str, secondary: list[str]) -> RecursiveDict:
    """Populates the output dict based on the provided arguments."""

    # create the outputs dictionary
    outputs: RecursiveDict = {
        type: {
            'basename': primary,
        },
    }

    # take the list of key=value strings, snap'em, and add each to the secondary files list
    if secondary:
        secondary_kv = defaultdict(dict)
        for keyvaluepair in secondary:
            key, value = keyvaluepair.split('=')
            secondary_kv[key] = {
                'basename': value,
            }

        outputs[type]['secondary_files'] = secondary_kv

    return outputs


def main():
    parser = argparse.ArgumentParser(description='Record MultiQC results to Metamist.')
    parser.add_argument('--project', required=True, help='The Metamist project name.')
    parser.add_argument('--output', required=True, help='Path of the primary output file/dir.')
    parser.add_argument('--type', required=True, help='Analysis type, must be from valid enum.')
    parser.add_argument('--secondary', nargs='+', help='Optional, list of k=v pairs for secondary analyses.', default=[])
    parser.add_argument('--meta', nargs='+', help='List of k=v metadata pairs, optional.', default=[])
    parser.add_argument('--cohorts', nargs='+', help='Metamist cohort ID(s).')
    parser.add_argument('--sgs', nargs='+', help='Metamist Sequencing Group ID(s).')
    parser.add_argument('--dry', action='store_true', help='Dry run, print only.')
    args = parser.parse_args()

    if args.cohorts and args.sgs:
        raise Exception('Cannot specify both --cohorts and --sgs')

    query = """
    mutation updateAnalysis($project: String!, $analysis:AnalysisInput!) {
        analysis {
            createAnalysis(project:$project, analysis:$analysis) {
                id
                type
                status
                output
                active
            }
        }
    }
    """

    outputs = create_output_block(
        primary=args.output,
        type=args.type,
        secondary=args.secondary,
    )

    variables = {
        'project': args.project,
        'analysis': {
            'type': args.type,
            'status': 'COMPLETED',
            'outputs': outputs,
            'meta': {}
        }
    }

    # allow arbitrary values to be passed into meta. Not sure here how best to tolerate numerical values, if we have 'em
    if args.meta:
        meta_kv: dict[str, str] = {}
        for keyvaluepair in args.meta:
            key, value = keyvaluepair.split('=')
            meta_kv[key] = value
        variables['analysis']['meta'] = meta_kv

    if args.cohorts:
        variables['analysis']['cohortIds'] = args.cohorts
    if args.sgs:
        variables['analysis']['sequencingGroupIds'] = args.sgs

    print(json.dumps(variables, indent=2))

    if args.dry:
        sys.exit(0)

    try:
        result = graphql.query(query, variables)
        print(f'Successfully recorded analysis to Metamist: {result}')
    except Exception as e:
        print(f'Error executing GraphQL query: {e}', file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()
