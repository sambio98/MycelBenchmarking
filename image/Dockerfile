# NPBench-C Phase 1 pinned-tool image.
#
# Every tool a build:image campaign needs, at an exact version AND build string,
# resolved once and recorded in environment.lock.json with artifact hashes. A
# score is comparable only against a run of the same image, so the pins are
# requirements rather than suggestions and the lockfile is what makes "the same
# image" checkable.
#
# Ubuntu's archive cannot satisfy two of the pins -- it ships DIAMOND 2.1.9 where
# the determinism controls require >= 2.2.7, and BLAST+ 2.12.0 where the pin is
# 2.17.0 -- so the tools come from conda-forge and bioconda via micromamba.
#
#   docker build -f image/Dockerfile -t npbench-c-tools:phase1 .
#   docker run --rm npbench-c-tools:phase1 \
#       python -m npbench_c.tools.invariance --all
#
# The image deliberately does NOT contain antiSMASH or BiG-SCAPE. Their pins are
# pairs (binary plus database release) that have to be resolved together, and
# they are the next image task; see environment.lock.json#not_in_this_image.
#
# STATUS: this file has not been built end to end. The build host available while
# it was written had the Docker CLI but no daemon. What HAS been run is every step
# it performs, directly: the micromamba download and the tar form below, the
# create with these exact version=build specs, `registry --verify` against the
# resulting prefix, and the full invariance suite -- whose measured verdicts are in
# image/thread_invariance.json. So the pins and the determinism results are real
# and the layer sequence is not yet proven. Build it once before relying on it.

FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive \
    MAMBA_ROOT_PREFIX=/opt/mamba \
    TOOL_PREFIX=/opt/npbench-tools \
    PATH=/opt/npbench-tools/bin:/usr/local/bin:/usr/bin:/bin \
    # Determinism: no tool in this image may pick its own thread count from the
    # host, and no BLAS may reorder a reduction behind our back. The per-tool
    # thread flags in npbench_c.tools.registry are still required and explicit;
    # these are the backstop for anything that reads the environment instead.
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    NUMEXPR_NUM_THREADS=1 \
    PYTHONHASHSEED=0 \
    PYTHONDONTWRITEBYTECODE=1 \
    LC_ALL=C \
    TZ=UTC

RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates curl bzip2 \
 && rm -rf /var/lib/apt/lists/*

# micromamba is pinned too: a different solver can resolve a different graph.
ARG MICROMAMBA_VERSION=2.9.0
RUN curl -sSL \
      "https://micro.mamba.pm/api/micromamba/linux-64/${MICROMAMBA_VERSION}" \
      -o /tmp/micromamba.tar.bz2 \
 && tar -xjf /tmp/micromamba.tar.bz2 -C /usr/local bin/micromamba \
 && rm /tmp/micromamba.tar.bz2 \
 && micromamba --version

# Exact version AND build string for every pin. Resolved 2026-10-02; the full
# 79-package closure with artifact hashes is in image/environment.lock.json.
RUN micromamba create -y -p "${TOOL_PREFIX}" --no-rc \
      -c conda-forge -c bioconda \
      'hmmer=3.4=h7d74f8d_5' \
      'diamond=2.2.8=he361c42_0' \
      'mmseqs2=18.8cc5c=hd6d6fdc_0' \
      'prodigal=2.6.3=h577a1d6_11' \
      'mafft=7.526=ha3fa977_0' \
      'blast=2.17.0=hb02a186_1' \
 && micromamba clean -a -y

# The grading core is stdlib-only by design, so the campaign runtime needs the
# system Python and nothing from PyPI.
RUN apt-get update \
 && apt-get install -y --no-install-recommends python3 python3-yaml \
 && rm -rf /var/lib/apt/lists/* \
 && ln -sf /usr/bin/python3 /usr/local/bin/python

WORKDIR /opt/npbench-c
COPY src/ src/
COPY image/environment.lock.json image/
ENV PYTHONPATH=/opt/npbench-c/src

# Fail the build if a pin drifted, rather than discovering it in a score.
RUN python -m npbench_c.tools.registry --verify

LABEL org.opencontainers.image.title="NPBench-C Phase 1 tools" \
      org.opencontainers.image.description="Pinned bioinformatics tools for \
NPBench-C build:image campaigns, with a thread-invariance suite" \
      org.opencontainers.image.licenses="see image/environment.lock.json"
