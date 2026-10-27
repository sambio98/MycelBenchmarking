# Environment

Pinned to the exact versions that produced `truth/`. A GECCO cluster inventory
is only valid against the release that generated it, so these are not advisory.

```bash
docker build -t npbench-r-02 environment/
# or, without docker:
pip install --no-deps -r environment/requirements.lock
gecco --version   # must print 0.11.0
```

`requirements.lock` is the runtime dependency closure of `gecco`, recorded from
a working install. Regenerate it after a deliberate version bump with:

```bash
python scripts/lock_task_environment.py tasks/02_salinispora_bgc gecco gecco-tool
```

Then regenerate `truth/` and re-run `make selfgrade`: a version bump invalidates
the old truth, because GECCO's cluster counts and type assignments move between
releases.
