# NPBench-C releases

The shippable tree is **generated, never committed**: it is derived from
`campaigns/`, and `npbench_c.release.package` refuses to write it inside this
repository so that the next packaging run cannot package its own output.

```bash
python3 -m npbench_c.release.package /somewhere/outside/npbench-c-1.0.0
```

What that produces, and what it withholds, is recorded here per release:

| file | what it is |
|---|---|
| `MANIFEST-<version>.json` | every shipped file with its sha256, the grader version, the campaigns included, the paths withheld |
| `ATTRIBUTION-<version>.md` | the attribution assembled from each campaign's own `license_notice` |
| `npbench-c-<version>.tar.gz.sha256` | the digest of the packaged tarball |

**Release 1.0.0** — 25 campaigns, 171 files, grader `2.0.0+9fc95658`. Rebuilt
twice on 2026-10-27, for the two output-contract findings: nine `task.yaml`
files moved the first time and twenty-two the second, and nothing else in the
tree either time.

Withheld from every campaign: `gold/`, `oracle/`, `oracle_submission/`,
`grading.json`, `internal_sweep.json`. Excluded trees: `docs/`, `src/`,
`tests/`, `image/`, `audit/`, `oracle/`. A campaign in the release is
`inputs/`, `reference/`, `task.yaml` and an empty `submission/` — the cold
sandbox a system is given, and nothing else.

Verify a tarball against its manifest:

```bash
python3 - <<'PY'
import hashlib, json, pathlib, sys
root = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else ".")
manifest = json.loads((root / "MANIFEST.json").read_text())
bad = []
for relative, digest in manifest["files"].items():
    path = root / relative
    if not path.is_file():
        bad.append(f"missing: {relative}")
        continue
    if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
        bad.append(f"digest differs: {relative}")
print(f"{len(manifest['files'])} files checked;", bad or "all match")
PY
```
