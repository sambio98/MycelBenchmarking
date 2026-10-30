"""Cut the campaign's input out of the full BRENDA release, once.

The release is 710 MB of JSON covering 8,129 EC numbers; this campaign needs the
three linked kinetic fields for one declared EC subclass. The extraction is
verbatim -- the value strings, the protein lists and the reference lists are
copied unchanged, so the agent parses exactly what BRENDA wrote -- and the result
is committed as the campaign's input with the source archive's hash in
provenance.

Run once, by hand, against a downloaded release:

  python -m npbench_c.templates.kinetics_consistency.prepare_input \
      <brenda.json> <ec_prefix> <out.json.gz>
"""

from __future__ import annotations

import gzip
import json
import pathlib
import sys

from npbench_c.templates.kinetics_consistency.core import FIELDS


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 3:
        print(__doc__.strip().splitlines()[-1], file=sys.stderr)
        return 2
    source, prefix, out = pathlib.Path(argv[0]), argv[1], pathlib.Path(argv[2])

    payload = json.loads(source.read_text())
    data = payload["data"]
    selected = {ec: {f: data[ec][f] for f in FIELDS if data[ec].get(f)}
                for ec in sorted(data) if ec.startswith(prefix)}
    selected = {ec: fields for ec, fields in selected.items() if fields}

    body = {
        "description": "Kinetic records cut verbatim from a BRENDA release. The "
                       "value strings are unchanged, so the parsing problem the "
                       "campaign poses is BRENDA's own.",
        "source": "BRENDA",
        "release": payload.get("release") or payload.get("version"),
        "ec_prefix": prefix,
        "fields": list(FIELDS),
        "license": "CC BY 4.0",
        "entries": len(selected),
        "data": selected,
    }
    with gzip.open(out, "wt", compresslevel=9) as handle:
        json.dump(body, handle, sort_keys=True)
    records = sum(len(v) for fields in selected.values() for v in fields.values())
    print(f"{out}: {len(selected)} EC numbers, {records} records, "
          f"{out.stat().st_size // 1024} KB gzipped")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
