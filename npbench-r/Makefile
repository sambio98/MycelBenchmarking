# Command surface from CLAUDE.md. `eval` belongs to build-order step 4 and
# fails loudly rather than doing something approximate.
TASK_DIR = tasks/$(TASK)
JOBS ?= 4
PY = python -m npbench.cli

.PHONY: test validate fetch truth selfgrade grade eval lock calibrate audit
.PHONY: require-task

require-task:
	@test -n "$(TASK)" || { echo "usage: make $(MAKECMDGOALS) TASK=<task dir name>"; exit 2; }

test:
	python -m pytest -q

validate: require-task
	$(PY) validate $(TASK_DIR)

fetch: require-task
	$(PY) fetch $(TASK_DIR)

truth: require-task
	$(PY) truth $(TASK_DIR) --jobs $(JOBS)

selfgrade: require-task
	$(PY) selfgrade $(TASK_DIR) --jobs $(JOBS)

grade: require-task
	@test -n "$(OUTPUTS)" || { echo "usage: make grade TASK=<dir> OUTPUTS=<dir>"; exit 2; }
	$(PY) grade $(TASK_DIR) $(OUTPUTS)

lock:
	@echo "Locks are per task: see tasks/<id>/environment/README.md"
	@exit 2

# The human curation pass. With no TASK it covers every task, which is how a
# reviewer reads it; exits non-zero on a finding nobody has judged yet.
audit:
	$(PY) audit $(if $(TASK),$(TASK_DIR),tasks/*/)

calibrate:
	@test -n "$(RUNS)" || { echo "usage: make calibrate RUNS=<dir of scored run JSON>"; exit 2; }
	$(PY) calibrate $(RUNS)

eval: require-task
	@test -n "$(MODEL)" || { echo "usage: make eval TASK=<dir> MODEL=<provider/model>"; exit 2; }
	inspect eval npbench/harness/task.py -T task_dir=$(TASK_DIR) --model $(MODEL)
