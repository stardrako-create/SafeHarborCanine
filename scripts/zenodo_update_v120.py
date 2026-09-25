#!/usr/bin/env python3
"""
Post a v1.20.0 status update to both Zenodo records (code + dataset) via
the REST API. Token read from the ZENODO_TOKEN environment variable only -
never hardcoded, never written to a file. Creates a new version draft,
appends to the existing description (never overwrites it), and publishes.
"""
import os
import sys

import requests

TOKEN = os.environ.get("ZENODO_TOKEN")
if not TOKEN:
    sys.exit("ZENODO_TOKEN not set in environment")

API = "https://zenodo.org/api/deposit/depositions"
HEADERS = {"Authorization": f"Bearer {TOKEN}"}


def new_version_and_update(record_id, extra_description_html):
    try:
        r = requests.post(f"{API}/{record_id}/actions/newversion", headers=HEADERS, timeout=30)
        r.raise_for_status()
        draft_url = r.json()["links"]["latest_draft"]
        draft = requests.get(draft_url, headers=HEADERS, timeout=30).json()
        draft_id = draft["id"]

        current_desc = draft["metadata"]["description"]
        new_desc = current_desc + extra_description_html

        r = requests.put(
            f"{API}/{draft_id}",
            headers=HEADERS,
            json={"metadata": {**draft["metadata"], "description": new_desc}},
            timeout=30,
        )
        r.raise_for_status()

        r = requests.post(f"{API}/{draft_id}/actions/publish", headers=HEADERS, timeout=30)
        r.raise_for_status()
        return r.json()["doi"], r.json()["doi_url"]
    except requests.exceptions.HTTPError as e:
        # sanitized - never print e (may embed request details) or the token
        raise SystemExit(f"Zenodo API error for record {record_id}: HTTP {e.response.status_code}")


CODE_RECORD_ID = "22055173"     # latest known version of the code deposit
DATASET_RECORD_ID = "22128237"  # latest known version of the dataset deposit

CODE_UPDATE = """
<p><strong>Update 2026-09-10</strong> - fixed a real bug in
score_ship_candidates.py's ATAC accessibility-floor veto: it compared a
window-averaged candidate value against a background built from individual
25bp bins, a statistical mismatch that made 254/461 candidates fail the
floor. Isolated with a 4-condition test and fixed with a dual check (wide
neighborhood window + narrow likely-insertion-point window, each against
its own matched-statistic background). Corrected result: 26/461 candidates
pass every criterion. Also completed v1.20.0, a from-scratch unified
rebuild (76-dog cohort, no post-hoc cohort joining, no raw.bw dependency)
that reproduces the same validated methodology - see the repository's
scripts/score_ship_candidates_v2.py and scripts/build_mother_track_v2.py
docstrings for full details.</p>
"""

DATASET_UPDATE = """
<p><strong>Update 2026-09-10</strong> - tracks rebuilt via v1.20.0 (unified
76-dog cohort for ATAC, no post-hoc cohort joining; RRBS confidence-gate
fix carried over). See the code repository's CHANGELOG/commit history for
the accessibility-floor statistical fix applied downstream in scoring.</p>
"""

if __name__ == "__main__":
    code_doi, code_url = new_version_and_update(CODE_RECORD_ID, CODE_UPDATE)
    print(f"Code record updated: {code_doi} ({code_url})")
    dataset_doi, dataset_url = new_version_and_update(DATASET_RECORD_ID, DATASET_UPDATE)
    print(f"Dataset record updated: {dataset_doi} ({dataset_url})")
