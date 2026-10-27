# Pinned environment for tasks/04_mibig_nrps_pks_modules.
#
# RDKit ships manylinux wheels, so this task needs no conda channel and no
# compiler. Canonical SMILES output, substrate formula computation and
# exact mass are version-sensitive, so the pin is load-bearing:
# truth is only valid against the version named in requirements.lock.
FROM docker.io/library/python:3.13-slim-bookworm@sha256:93498f36dc7e24ef233bfdf6cc17d4b8bd32d94ee737143e906460d5f5a914c5

WORKDIR /task
COPY requirements.lock /task/environment/requirements.lock
RUN pip install --no-cache-dir --no-deps -r /task/environment/requirements.lock \
 && python -c "import rdkit; print(rdkit.__version__)"
