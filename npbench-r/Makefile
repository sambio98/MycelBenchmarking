# Command surface from CLAUDE.md. Targets for build-order steps that are not
# built yet fail loudly rather than doing something approximate.
TASK_DIR = tasks/$(TASK)

.PHONY: test validate grade lock fetch truth selfgrade eval

test:
	python -m pytest -q

validate:
	@test -n "$(TASK)" || { echo "usage: make validate TASK=01"; exit 2; }
	python -m npbench.cli validate $(TASK_DIR)

grade:
	@test -n "$(TASK)" || { echo "usage: make grade TASK=01 OUTPUTS=path/to/outputs"; exit 2; }
	python -m npbench.cli grade $(TASK_DIR) $(OUTPUTS)

lock fetch truth selfgrade eval:
	@echo "'$@' is not implemented. Build order (CLAUDE.md):"
	@echo "  lock, fetch, truth, selfgrade  -> step 3, one task end to end"
	@echo "  eval                           -> step 4, the Inspect AI harness"
	@exit 2
