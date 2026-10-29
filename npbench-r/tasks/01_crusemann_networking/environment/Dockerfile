# Pinned environment for tasks/01_crusemann_networking.
#
# matchms installs from wheels and needs no conda channel and no compiler. The
# pin is load-bearing beyond reproducibility: matchms renamed ModifiedCosine to
# ModifiedCosineGreedy between releases, and cosine scores are
# implementation-dependent, so the reference solution refuses to run against
# any other version rather than quietly producing different truth.
FROM docker.io/library/python:3.13-slim-bookworm@sha256:93498f36dc7e24ef233bfdf6cc17d4b8bd32d94ee737143e906460d5f5a914c5

WORKDIR /task
COPY requirements.lock /task/environment/requirements.lock
RUN pip install --no-cache-dir --no-deps -r /task/environment/requirements.lock \
 && python -c "import matchms, networkx; print(matchms.__version__, networkx.__version__)"
