# Pinned environment for tasks/09_salinispora_paired_omics.
#
# GECCO is pure Python with wheels for pyrodigal and pyhmmer, so this task needs
# no conda channel and no compiler. That is why GECCO was chosen over antiSMASH
# for the first BGC-detection task: PROJECT.md section 4 lists both, and GECCO
# is the one that pins reproducibly without a 9 GB database download.
FROM docker.io/library/python:3.13-slim-bookworm@sha256:93498f36dc7e24ef233bfdf6cc17d4b8bd32d94ee737143e906460d5f5a914c5

WORKDIR /task
COPY requirements.lock /task/environment/requirements.lock
RUN pip install --no-cache-dir --no-deps -r /task/environment/requirements.lock \
 && gecco --version \
 && python -c "import matchms, networkx; from pyteomics import version; print(matchms.__version__, networkx.__version__, version.version)"
