#!/usr/bin/env python3

"""
Queries Metamist for all gVCFs in the cohort(s), writes a CSV samplesheet.
"""

import csv
import argparse
import hashlib
import sys

from metamist import graphql


ANALYSIS_QUERY = """
query AnalysisQuery($project: String!, $cohortFilter: StrGraphQLFilter, $a_type: StrGraphQLFilter) {
  project(name: $project) {
    cohorts(id: $cohortFilter) {
      id
      sequencingGroups {
        id
        analyses(active: {eq: true}, type: $a_type) {
          output
        }
      }
    }
  }
}
"""


def get_sg_hash(sg_ids: list[str]) -> str:
    """Unique hash string from Sequencing Group IDs."""
    h = hashlib.sha256(''.join(sg_ids).encode()).hexdigest()[:38]
    return f'{h}_{len(sg_ids)}'


def try_query_with_variables(query_string: str, variables: dict[str, str | dict[str, str]]):
    """Posts a populated query string with variables."""
    try:
        return graphql.query(query_string, variables)
    except Exception as e:
        print(f"Error executing GraphQL query: {e}", file=sys.stderr)
        sys.exit(1)


def build_samplesheet_rows(project: str, cohort_ids: list[str]) -> tuple[list[dict], str]:
    variables = {
        "project": project,
        "cohortFilter": {'in_': cohort_ids},
        "a_type": {'eq': 'gvcf'},
    }
    result = try_query_with_variables(ANALYSIS_QUERY, variables)

    seen_sgs: set[str] = set()
    rows = []

    for cohort in result['project']['cohorts']:
        cohort_id = cohort['id']
        sg_ids = [sg['id'] for sg in cohort['sequencingGroups']]
        cohort_hash = get_sg_hash(sorted(sg_ids))

        for sg in cohort['sequencingGroups']:
            sg_id = sg['id']

            if sg_id in seen_sgs:
                print(f"Error: SG {sg_id} appears in multiple cohorts", file=sys.stderr)
                sys.exit(1)
            seen_sgs.add(sg_id)

            if not sg['analyses']:
                print(f"Error: SG {sg_id} in cohort {cohort_id} has no gVCF", file=sys.stderr)
                sys.exit(1)

            rows.append({
                'sg_id': sg_id,
                'gvcf': sg['analyses'][0]['output'],
                'project': project,
                'cohort': cohort_id,
                'cohort_hash': cohort_hash,
            })

    if not rows:
        print("Error: no sequencing groups found", file=sys.stderr)
        sys.exit(1)

    return rows, get_sg_hash(sorted(seen_sgs))


def main():
    parser = argparse.ArgumentParser(description="Fetch gVCF paths from Metamist and write a samplesheet CSV.")
    parser.add_argument("--project", required=True, help="The Metamist project name")
    parser.add_argument("--cohorts", required=True, nargs="+", help="Metamist cohort IDs")
    args = parser.parse_args()

    # get all samplesheet rows, and one hash spanning all samples/SG IDs. This is a proxy for output_version in CPG-Flow
    rows, big_hash = build_samplesheet_rows(args.project, args.cohorts)

    fieldnames = ['sg_id', 'gvcf', 'project', 'cohort', 'cohort_hash']
    with open(f'{big_hash}.tsv', 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {big_hash}.tsv")


if __name__ == "__main__":
    main()
