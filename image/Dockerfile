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
# antiSMASH's Pfam library rather than fetching its own.
#
# matchms is requested below and is this image's first LIBRARY pin rather than a
# binary. That distinction is load-bearing: see the comment on the create step.
#
# It is also where the staged build bit. The prefix assembled in stages carried
# matchms WITHOUT `requests`, which matchms declares and imports at import time,
# and `registry --verify` reported it present anyway because the version probe
# ran an import that could also see the surrounding container's user site, where
# an unrelated tool had requests installed. A fresh `create` of the same spec
# links the whole chain -- measured -- so the single create below is the complete
# path and the staged prefix was the incomplete one. The probe now runs with
# `-s`, and `--verify` checks every library pin's declared run dependencies
# against the prefix, so neither form of this can pass again.
#
# STATUS: built, with one layer skipped. 2026-10-27, on a host with 9 GB free
# and an agent proxy in front of the network: every layer below built and
# `registry --verify` passed inside the resulting image -- 9/9 tools at their
# pinned version and build, FastTreeMP absent, matchms' run closure complete.
# What is NOT yet proven is the `download-antismash-databases` layer, which is
# roughly 10 GB and did not fit; it was commented out for that build and is the
# one step still to run before relying on the S1/S2 campaigns. The build also
# needed two host-specific workarounds that are NOT in this file and should not
# be: the proxy tunnels HTTPS and answers plain HTTP with 405, and the base
# image carries no CA store, so apt needed its CaInfo pointed at the proxy
# bundle and the archive URIs moved to https. On a host with direct network
# access neither is necessary.
#
# What that build found is the reason the create step below reads the way it
# does. The prefix had been built in STAGES rather than in one create: six tools
# on 2026-10-02, antiSMASH, BiG-SCAPE and the sqlalchemy bound on 10-06, matchms
# on 10-21, each step recorded with its date in the lock's `build_history`. The
# comment here used to say that a closure assembled over five solves is not
# guaranteed to equal what one solve produces, and that the difference should be
# visible rather than implied. It is now measured: a single fresh solve of the
# eleven pins produced 232 packages against the lock's 231, with cpython,
# curl, brotli, ca-certificates, fonttools, importlib-metadata, bitarray and
# c-ares all newer than the lock records. So this file no longer asks the solver
# anything -- it installs the locked closure by URL and sha256, and a rebuild of
# it was verified to produce those 231 packages exactly.
# The matchms step was dry-run first and moved nothing already installed: 46 new
# packages, no upgrade, no downgrade, no removal.

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
# matchms is the odd one in this list and the reasons are measured rather than
# stylistic. The version is 0.30.2 because that is the version the MS2 fixture
# grounding ran against -- bioconda carries it as a noarch build, so the spec is
# exact like every other line here, and pinning a later one would mean re-running
# the tolerance measurement rather than inheriting it. It is a LIBRARY, so it has
# no binary to invoke and no --version to read: a campaign reaches it by running
# this prefix's python3.11 as a subprocess, which is why PATH above keeps the
# prefix last and `python3` the system interpreter. Importing matchms into the
# campaign runtime instead would put numpy, scipy and numba inside a grading core
# that is stdlib-only by design. And its m/z tolerance is part of gold rather
# than a comparison detail -- 20 of 150 per-query verdicts move across the
# plausible range -- so a campaign pins one value and grades the sweep.
#
# FastTreeMP is deleted immediately after: it arrives inside the fasttree package
# that antiSMASH and BiG-SCAPE both depend on, and it is banned outright because
# thread order affects the neighbour-joining heuristic, so it can never sit on a
# gold-producing path. Deleting the binary is what makes the ban enforceable --
# `registry --verify` below fails if it comes back. Single-threaded FastTree
# stays, because both tools need it.
#
# The create below installs the LOCKED CLOSURE, not the eleven specs above: an
# explicit file of 231 artifact URLs, each with the sha256 the lock recorded,
# generated from environment.lock.json by
# `python -m npbench_c.tools.lockfile --explicit`. The eleven specs stay in the
# lock's `requested` field, which is what they are -- the intent -- and a test
# fails if the explicit file and the lock disagree.
#
# This was a spec list until it was built. A fresh solve of the same eleven pins
# resolved 232 packages against the lock's 231, with cpython 3.11.17 -> .18,
# curl 8.22 -> 8.23, brotli 1.2.0 -> 1.2.1, ca-certificates, fonttools,
# importlib-metadata, bitarray, c-ares and more all newer than the lock. Every
# pinned tool was still at its pinned version and build, so `registry --verify`
# passed and said nothing about it: the pins held and the closure under them
# moved. A score is comparable only against a run of the same image, so the
# closure is the thing that has to be reproducible, and a spec list does not
# reproduce it.
COPY image/environment.explicit.txt /tmp/environment.explicit.txt
RUN micromamba create -y -p "${TOOL_PREFIX}" --no-rc \
      -c conda-forge -c bioconda \
      --file /tmp/environment.explicit.txt \
 && rm -f /tmp/environment.explicit.txt \
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
