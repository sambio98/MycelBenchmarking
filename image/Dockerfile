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
# antiSMASH is here, and its pin is really a PAIR: this binary plus the database
# release downloaded in its own layer below. Both determine which regions come
# out, so neither alone identifies a result. BiG-SCAPE is here too and reuses
# antiSMASH's Pfam library rather than fetching its own. matchms is still absent,
# for a recorded reason; see environment.lock.json.
#
# STATUS: this file has not been built end to end. The build host available while
# it was written had the Docker CLI but no daemon. What HAS been run is every step
# it performs, directly: the micromamba download and the tar form below, the
# create with these exact version=build specs, `registry --verify` against the
# resulting prefix, and the full invariance suite -- whose measured verdicts are in
# image/thread_invariance.json. So the pins and the determinism results are real
# and the layer sequence is not yet proven. Build it once before relying on it.

FROM ubuntu:24.04

# Two things in the block below are deliberate rather than conventional.
#
# PATH puts the tool prefix LAST. antiSMASH and BiG-SCAPE bring their own Python
# into the prefix, and their console scripts carry an absolute shebang to it, so
# the prefix does not need to precede the system interpreter -- while if it does,
# `python3` silently becomes their interpreter and the benchmark's stdlib-only
# grading core starts running on an interpreter nobody chose. Measured the hard
# way: with the prefix first, `python3 -m pytest` could not find pytest.
#
# The thread environment is a backstop. No tool here may pick its own thread
# count from the host, and no BLAS may reorder a reduction behind our back. The
# per-tool thread flags in npbench_c.tools.registry are still required and
# explicit; these catch anything that reads the environment instead.
ENV DEBIAN_FRONTEND=noninteractive \
    MAMBA_ROOT_PREFIX=/opt/mamba \
    TOOL_PREFIX=/opt/npbench-tools \
    PATH=/usr/local/bin:/usr/bin:/bin:/opt/npbench-tools/bin \
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

# Exact version AND build string for every pin. Resolved 2026-10-06; the full
# 185-package closure with artifact hashes is in image/environment.lock.json.
#
# The sqlalchemy pin is an upper bound the upstream recipe does NOT carry.
# BiG-SCAPE 2.0.3 asks for sqlalchemy >= 2.0.2, a fresh solve gives 2.1.x, and it
# then raises ObjectNotExecutableError before reading a single input file.
# Measured by running it, not read off a changelog.
#
# FastTreeMP is deleted immediately after: it arrives inside the fasttree package
# that antiSMASH and BiG-SCAPE both depend on, and it is banned outright because
# thread order affects the neighbour-joining heuristic, so it can never sit on a
# gold-producing path. Deleting the binary is what makes the ban enforceable --
# `registry --verify` below fails if it comes back. Single-threaded FastTree
# stays, because both tools need it.
RUN micromamba create -y -p "${TOOL_PREFIX}" --no-rc \
      -c conda-forge -c bioconda \
      'hmmer=3.4=h7d74f8d_5' \
      'diamond=2.2.8=he361c42_0' \
      'mmseqs2=18.8cc5c=hd6d6fdc_0' \
      'prodigal=2.6.3=h577a1d6_11' \
      'mafft=7.526=ha3fa977_0' \
      'blast=2.17.0=hb02a186_1' \
      'antismash=8.0.4=pyhdfd78af_1' \
      'bigscape=2.0.3=pyhdfd78af_0' \
      'sqlalchemy=2.0.54=py311h194be0f_1' \
 && micromamba clean -a -y \
 && rm -f "${TOOL_PREFIX}/bin/FastTreeMP"

# The antiSMASH reference databases, in their own layer because they are an order
# of magnitude larger than everything above and they version independently of the
# binary. The default location is inside the package, which is where the suite
# expects them; for distribution this layer can be dropped and the same path
# mounted instead, since antiSMASH also takes --databases PATH.
#
# Size, measured: roughly 10 GB on disk once ClusterBlast 4.0 is extracted. That
# is the cost of a BGC-detection campaign and it is better stated here than
# discovered during a build.
RUN download-antismash-databases \
 && find "${TOOL_PREFIX}/lib" -name '*.tar' -delete \
 && find "${TOOL_PREFIX}/lib" -name '*.tar.xz' -delete

# The grading core is stdlib-only by design, so the campaign runtime needs the
# system Python and nothing from PyPI.
RUN apt-get update \
 && apt-get install -y --no-install-recommends python3 python3-yaml \
 && rm -rf /var/lib/apt/lists/* \
 && ln -sf /usr/bin/python3 /usr/local/bin/python \
 && ln -sf /usr/bin/python3 /usr/local/bin/python3

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
