"""Cut the campaign's input out of the ChEMBL web service, once.

Fetches every activity on the two declared targets plus the assay records the
admission rule needs, and writes them as one compact file. The activity fields are
copied verbatim, so the filtering problem the campaign poses is ChEMBL's own.

The output is a ChEMBL derivative and therefore carries CC BY-SA 3.0; the campaign
declares `license_class: share_alike` and ships the notice. See
docs/environment.md for why that obligation is local to the campaign.

  python -m npbench_c.templates.chembl_selectivity.prepare_input \
      <focal_target> <reference_target> <out.json.gz>
"""

from __future__ import annotations

import gzip
import json
import pathlib
import sys
import urllib.parse
import urllib.request

BASE = "https://www.ebi.ac.uk/chembl/api/data"
ACTIVITY_FIELDS = (
    "activity_id,molecule_chembl_id,molecule_pref_name,target_chembl_id,"
    "target_pref_name,target_organism,assay_chembl_id,standard_type,"
    "standard_relation,standard_value,standard_units,pchembl_value,"
    "data_validity_comment,potential_duplicate,document_chembl_id,document_year")
ASSAY_FIELDS = ("assay_chembl_id,confidence_score,confidence_description,"
                "assay_type,assay_organism")
TIMEOUT = 300


def _pages(path: str, params: dict, key: str, cap: int = 20000) -> list[dict]:
    rows: list[dict] = []
    url = f"{BASE}/{path}.json?" + urllib.parse.urlencode(params)
    while url and len(rows) < cap:
        with urllib.request.urlopen(url, timeout=TIMEOUT) as handle:
            payload = json.load(handle)
        rows += payload[key]
        nxt = payload["page_meta"].get("next")
        url = ("https://www.ebi.ac.uk" + nxt) if nxt else None
    return rows


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 3:
        print(__doc__.strip().splitlines()[-1], file=sys.stderr)
        return 2
    focal, reference, out = argv[0], argv[1], pathlib.Path(argv[2])

    with urllib.request.urlopen(f"{BASE}/status.json", timeout=TIMEOUT) as handle:
        status = json.load(handle)

    activities: list[dict] = []
    for target in (focal, reference):
        activities += _pages("activity", {"target_chembl_id": target,
                                          "limit": 1000, "only": ACTIVITY_FIELDS},
                             "activities")
    activities.sort(key=lambda a: a["activity_id"])

    ids = sorted({a["assay_chembl_id"] for a in activities if a.get("assay_chembl_id")})
    assays: dict[str, dict] = {}
    for i in range(0, len(ids), 50):
        for assay in _pages("assay", {"assay_chembl_id__in": ",".join(ids[i:i + 50]),
                                      "limit": 100, "only": ASSAY_FIELDS}, "assays"):
            assays[assay["assay_chembl_id"]] = assay

    body = {
        "description": "ChEMBL activities and assay records for two targets, cut "
                       "verbatim so the filtering problem is ChEMBL's own.",
        "source": "ChEMBL",
        "release": status["chembl_db_version"],
        "release_date": status["chembl_release_date"],
        "license": "CC BY-SA 3.0",
        "focal_target": focal,
        "reference_target": reference,
        "activities": activities,
        "assays": dict(sorted(assays.items())),
    }
    with gzip.open(out, "wt", compresslevel=9) as handle:
        json.dump(body, handle, sort_keys=True)
    print(f"{out}: {status['chembl_db_version']}, {len(activities)} activities, "
          f"{len(assays)} assays, {out.stat().st_size // 1024} KB gzipped")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
