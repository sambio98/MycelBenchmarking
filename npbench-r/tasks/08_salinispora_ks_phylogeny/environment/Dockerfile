# Pinned environment for tasks/08_salinispora_ks_phylogeny.
#
# pyrodigal, pyhmmer and scipy all install from wheels, so this task needs no
# conda channel and no compiler. gecco is installed for one reason only: it
# ships the Pfam 35.0 subset that holds the ketosynthase model, so the task
# fetches no database at runtime. The pins are load-bearing -- gene calls and
# domain boundaries both move between releases.
FROM docker.io/library/python:3.13-slim-bookworm@sha256:93498f36dc7e24ef233bfdf6cc17d4b8bd32d94ee737143e906460d5f5a914c5

WORKDIR /task
COPY requirements.lock /task/environment/requirements.lock
RUN pip install --no-cache-dir --no-deps -r /task/environment/requirements.lock \
 && python -c "import pyrodigal, pyhmmer, scipy, Bio; print(pyrodigal.__version__, pyhmmer.__version__, scipy.__version__, Bio.__version__)"
