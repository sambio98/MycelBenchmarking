# Environment

Pinned to the exact versions that produced `truth/`. A GECCO cluster inventory
is only valid against the release that generated it, so these are not advisory.

```bash
docker build -t npbench-r-02 environment/
# or, without docker:
pip install --no-deps -r environment/requirements.lock
gecco --version   # must print 0.11.0
```

`requirements.lock` is the resolved dependency closure of `gecco`, recorded from
a working install. To regenerate it after a deliberate version bump, install the
new `gecco` into a clean environment and re-run the closure walk in
`scripts/` — then regenerate `truth/` and re-run `make selfgrade`, because the
old truth is no longer valid.
