# Command surface from CLAUDE.md. `eval` belongs to build-order step 4 and
# fails loudly rather than doing something approximate.
TASK_DIR = tasks/$(TASK)
JOBS ?= 4
PY = python -m npbench.cli

.PHONY: test validate fetch truth selfgrade grade eval lock
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

eval:
	@echo "'eval' needs the Inspect AI harness (CLAUDE.md build order, step 4)."
	@exit 2
