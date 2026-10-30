"""Measured difficulty of tasks/09_salinispora_paired_omics.

Three layers, and they cost very different amounts, so everything expensive is
cached on disk and keyed by the configuration that produced it:

  * **The genomic layer** reuses stage 1's GECCO output, which no modelled
    mistake changes -- the mistakes are in how clusters are grouped, not in
    whether they are found. Free.
  * **The metabolomic layer** changes the presence calls, which is ~9 minutes
    of modified-cosine over 111 strains. Cached per presence configuration.
  * **The consensus network** is ~100 seconds and changes only if the spectrum
    filter does. Cached too.
  * **The linking layer** is arithmetic over cached occupancy sets. Free.

Variants are ordered cheap-first so a partial run still says something about
every layer. Same cheap-proxy caveat as the other tasks.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import pickle
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import networkx as nx  # noqa: E402
import numpy as np  # noqa: E402
from matchms.similarity import CosineGreedy, ModifiedCosineGreedy  # noqa: E402

from npbench.grading import load_grading_config, score_task  # noqa: E402

# Grouping the clusters and naming the families: free, reuses cached GECCO.
GENOMIC = {
    "domains_as_multiset": "uses the raw domains field without deduplicating it",
    "keeps_thin_clusters": "skips the fewer-than-3-distinct-domains exclusion",
    "families_numbered_by_rank": "names families GCF_0001.. by size rank, not by member",
    "jaccard_0_30": "joins clusters at Jaccard 0.30 instead of 0.50",
    "pairwise_not_transitive": "treats each linked pair as a family, no single linkage",
}

# Changing the presence calls: ~9 minutes each, cached.
# The top-50 peak cut is NOT modelled here. Changing it means re-extracting
# from the raw mzXML, and this bundle streams those: each file is fetched,
# fingerprinted, reduced to its filtered spectra and deleted, because 111 of
# them do not fit on one container's disk beside everything else. An auditor
# with the full staging can add the lever; the gap is recorded in NOTES.md
# rather than left as a variant that silently measures nothing.
METABOLOMIC = {
    "presence_window_2da": "allows a 2.0 Da precursor difference instead of 0.5",
    "presence_cosine_0_50": "calls a spectrum present at cosine 0.50 instead of 0.70",
    "plain_cosine": "uses CosineGreedy instead of ModifiedCosineGreedy",
}

# Combining and ranking: free.
LINKING = {
    "family_presence_needs_all": "requires every spectrum of a family, not any one",
    "metcalf_symmetric": "scores a cluster-only strain 0 rather than -10",
    "metcalf_no_neither": "drops the +1 for a strain with neither",
    "scores_undetected_families": "scores families present in no strain",
    "rank_not_dense": "ranks 1..500 ordinally, so tied scores take different ranks",
}

MISTAKES = {**GENOMIC, **METABOLOMIC, **LINKING}


def load_solution(task_dir: Path):
    spec = importlib.util.spec_from_file_location(
        "ref_solution_09", task_dir / "solution" / "make_truth.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def presence_for(ref, task_dir: Path, consensus, strains, mistakes: set[str],
                 scratch: Path) -> dict[int, set[str]]:
    """Presence calls, cached per presence configuration."""
    key = (
        "presence_window_2da" in mistakes,
        "presence_cosine_0_50" in mistakes,
        "plain_cosine" in mistakes,
    )
    cache = scratch / f"presence_{'_'.join(str(int(k)) for k in key)}.pkl"
    if cache.is_file():
        return pickle.loads(cache.read_bytes())

    window = 2.0 if key[0] else ref.PRESENCE_WINDOW_DA
    threshold = 0.50 if key[1] else ref.PRESENCE_MIN_COSINE
    similarity = (CosineGreedy(tolerance=ref.FRAGMENT_TOLERANCE) if key[2]
                  else ModifiedCosineGreedy(tolerance=ref.FRAGMENT_TOLERANCE))
    print(f"    (presence calls for {key})", flush=True)

    present: dict[int, set[str]] = defaultdict(set)
    for index, strain in enumerate(strains, start=1):
        accession = strain["accession"]
        spectra, _ = ref.strain_spectra(
            task_dir / "work" / "features" / f"{accession}.npz")
        if not spectra:
            continue
        masses = np.array([float(s.get("precursor_mz")) for s in spectra])
        for spectrum in consensus:
            mass = float(spectrum.get("precursor_mz"))
            lo, hi = np.searchsorted(masses, [mass - window, mass + window])
            for position in range(int(lo), int(hi)):
                if float(similarity.pair(spectrum, spectra[position])["score"]) >= threshold:
                    present[int(spectrum.get("scans"))].add(accession)
                    break
        if index % 25 == 0:
            print(f"      {index}/{len(strains)} strains", flush=True)
    cache.write_bytes(pickle.dumps(dict(present)))
    return dict(present)


def network_for(ref, consensus, scratch: Path):
    cache = scratch / "network.pkl"
    if cache.is_file():
        return pickle.loads(cache.read_bytes())
    print("    (consensus network)", flush=True)
    edges = ref.apply_top_k(ref.candidate_edges(consensus))
    cache.write_bytes(pickle.dumps(edges))
    return edges


def families_of(ref, domains_of, mistakes: set[str]) -> dict[str, str]:
    """Cluster -> family under one genomic configuration."""
    if "domains_as_multiset" in mistakes:
        # Jaccard over a multiset read as a set of (domain, position) pairs
        # inflates the union and so lowers every similarity.
        domains_of = {
            cluster: frozenset(f"{d}#{i}" for i, d in enumerate(sorted(ds)))
            for cluster, ds in domains_of.items()
        }
    threshold = 0.30 if "jaccard_0_30" in mistakes else ref.JACCARD_THRESHOLD
    floor = 0 if "keeps_thin_clusters" in mistakes else ref.MIN_DOMAINS

    eligible = sorted(c for c, d in domains_of.items() if len(d) >= floor)
    if "pairwise_not_transitive" in mistakes:
        # Each linked pair is its own family, named for the smaller member; a
        # cluster in several pairs lands in whichever is lexicographically
        # first, which is the point -- the grouping is no longer a partition.
        assignment = {c: c for c in eligible}
        for i, a in enumerate(eligible):
            for b in eligible[i + 1:]:
                da, db = domains_of[a], domains_of[b]
                union = len(da | db)
                if union and len(da & db) / union >= threshold:
                    assignment[b] = min(assignment[b], a)
        return assignment

    parent = {c: c for c in eligible}

    def find(c):
        while parent[c] != c:
            parent[c] = parent[parent[c]]
            c = parent[c]
        return c

    for i, a in enumerate(eligible):
        for b in eligible[i + 1:]:
            da, db = domains_of[a], domains_of[b]
            union = len(da | db)
            if union and len(da & db) / union >= threshold:
                ra, rb = find(a), find(b)
                if ra != rb:
                    parent[max(ra, rb)] = min(ra, rb)
    assignment = {c: find(c) for c in eligible}

    if "families_numbered_by_rank" in mistakes:
        sizes = Counter(assignment.values())
        order = sorted(sizes, key=lambda f: (-sizes[f], f))
        renamed = {f: f"GCF_{i:04d}" for i, f in enumerate(order, start=1)}
        assignment = {c: renamed[f] for c, f in assignment.items()}
    return assignment


def build(ref, task_dir: Path, out: Path, mistakes: set[str], scratch: Path,
          cached: dict) -> None:
    strains = cached["strains"]
    species = cached["species"]
    all_strains = set(cached["accessions"])
    w = ref.write_csv

    # --- genomic half --------------------------------------------------------
    w(out / "bgc_inventory.csv",
      ["cluster_id", "accession", "product_type", "n_domains"], cached["cluster_rows"])

    family_of = families_of(ref, cached["domains_of"], mistakes)
    strains_of_gcf, clusters_of_gcf = defaultdict(set), defaultdict(list)
    for cluster_id, family_id in family_of.items():
        strains_of_gcf[family_id].add(cached["strain_of_cluster"][cluster_id])
        clusters_of_gcf[family_id].append(cluster_id)
    w(out / "gene_cluster_families.csv",
      ["gcf_id", "n_clusters", "n_strains", "n_species", "representative_cluster"],
      [{"gcf_id": g, "n_clusters": len(clusters_of_gcf[g]),
        "n_strains": len(strains_of_gcf[g]),
        "n_species": len({species[a] for a in strains_of_gcf[g]}),
        "representative_cluster": min(clusters_of_gcf[g])}
       for g in sorted(strains_of_gcf)])

    # --- metabolomic half ----------------------------------------------------
    present = presence_for(ref, task_dir, cached["consensus"], strains, mistakes, scratch)
    strain_rows = []
    for strain in strains:
        accession = strain["accession"]
        found = sum(1 for scan in present if accession in present[scan])
        strain_rows.append({
            "accession": accession, "species": species[accession],
            "strain": strain["strain"],
            "n_ms2_scans": cached["n_ms2"][accession],
            "n_usable_ms2": cached["n_usable"][accession],
            "n_spectra_present": found,
        })
    w(out / "strain_inventory.csv",
      ["accession", "species", "strain", "n_ms2_scans", "n_usable_ms2",
       "n_spectra_present"], strain_rows)
    w(out / "spectrum_occurrence.csv", ["scan", "accession"],
      [{"scan": scan, "accession": a} for scan in sorted(present)
       for a in sorted(present[scan])])

    members_of_family = cached["members_of_family"]
    if "family_presence_needs_all" in mistakes:
        strains_of_family = {
            f: set.intersection(*(present.get(s, set()) for s in members)) if members else set()
            for f, members in members_of_family.items()
        }
    else:
        strains_of_family = {
            f: set().union(*(present.get(s, set()) for s in members)) if members else set()
            for f, members in members_of_family.items()
        }
    w(out / "molecular_families.csv", ["mf_id", "n_spectra", "n_strains", "n_species"],
      [{"mf_id": f, "n_spectra": len(members_of_family[f]),
        "n_strains": len(strains_of_family[f]),
        "n_species": len({species[a] for a in strains_of_family[f]})}
       for f in sorted(members_of_family)])

    # --- links ---------------------------------------------------------------
    both_w = ref.METCALF_BOTH
    gcf_w = 0 if "metcalf_symmetric" in mistakes else ref.METCALF_GCF_ONLY
    neither_w = 0 if "metcalf_no_neither" in mistakes else ref.METCALF_NEITHER

    def score(gcf_strains, family_strains):
        both = len(gcf_strains & family_strains)
        gcf_only = len(gcf_strains - family_strains)
        family_only = len(family_strains - gcf_strains)
        neither = len(all_strains) - both - gcf_only - family_only
        return (both * both_w + gcf_only * gcf_w
                + family_only * ref.METCALF_FAMILY_ONLY + neither * neither_w)

    if "scores_undetected_families" in mistakes:
        linkable = sorted(members_of_family)
    else:
        linkable = sorted(f for f in members_of_family if strains_of_family[f])
    scores = [
        {"gcf_id": g, "mf_id": f,
         "metcalf_score": score(strains_of_gcf[g], strains_of_family[f]),
         "n_shared_strains": len(strains_of_gcf[g] & strains_of_family[f])}
        for g in sorted(strains_of_gcf) for f in linkable
    ]
    w(out / "metcalf_scores.csv",
      ["gcf_id", "mf_id", "metcalf_score", "n_shared_strains"], scores)

    ordered = sorted(scores, key=lambda r: (-r["metcalf_score"], r["gcf_id"], r["mf_id"]))
    top = []
    if "rank_not_dense" in mistakes:
        top = [{**row, "rank": i} for i, row in enumerate(ordered[:ref.TOP_LINKS], start=1)]
    else:
        rank, previous = 0, None
        for row in ordered[:ref.TOP_LINKS]:
            if row["metcalf_score"] != previous:
                rank += 1
                previous = row["metcalf_score"]
            top.append({**row, "rank": rank})
    w(out / "top_links.csv",
      ["rank", "gcf_id", "mf_id", "metcalf_score", "n_shared_strains"], top)


VARIANTS = {
    "reference (no mistakes)": set(),
    "one slip: domains as a multiset": {"domains_as_multiset"},
    "one slip: keeps thin clusters": {"keeps_thin_clusters"},
    "one slip: families numbered by rank": {"families_numbered_by_rank"},
    "one slip: Jaccard 0.30": {"jaccard_0_30"},
    "one slip: pairs, not single linkage": {"pairwise_not_transitive"},
    "one slip: family needs every spectrum": {"family_presence_needs_all"},
    "one slip: symmetric Metcalf": {"metcalf_symmetric"},
    "one slip: no +1 for neither": {"metcalf_no_neither"},
    "one slip: scores undetected families": {"scores_undetected_families"},
    "one slip: ordinal rank, not dense": {"rank_not_dense"},
    "one slip: 2.0 Da presence window": {"presence_window_2da"},
    "one slip: presence at cosine 0.50": {"presence_cosine_0_50"},
    "one slip: plain cosine": {"plain_cosine"},
    "hasty agent (3 slips)": {"families_numbered_by_rank", "metcalf_symmetric",
                              "rank_not_dense"},
    "careless agent (5 slips)": {"keeps_thin_clusters", "jaccard_0_30",
                                 "presence_cosine_0_50", "metcalf_no_neither",
                                 "scores_undetected_families"},
}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--task", type=Path, required=True)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--scratch", type=Path, required=True)
    parser.add_argument("--only", default=None)
    args = parser.parse_args()
    args.scratch.mkdir(parents=True, exist_ok=True)

    config = load_grading_config(args.task / "grading.yaml")
    ref = load_solution(args.task)
    work = args.task / "work"

    strains = ref.load_strains(args.task / "strains.tsv")
    species = {s["accession"]: ref.species_of(s["organism"]) for s in strains}
    domains_of, cluster_rows, strain_of_cluster = {}, [], {}
    for strain in strains:
        accession = strain["accession"]
        for record in ref.read_clusters(work, accession):
            cluster_id = record["cluster_id"]
            domains = ref.domain_set(record.get("domains", ""))
            domains_of[cluster_id] = domains
            strain_of_cluster[cluster_id] = accession
            cluster_rows.append({"cluster_id": cluster_id, "accession": accession,
                                 "product_type": record.get("type") or "Unknown",
                                 "n_domains": len(domains)})
    cluster_rows.sort(key=lambda r: r["cluster_id"])

    consensus = ref.load_consensus(args.task / ref.MGF)
    edges = network_for(ref, consensus, args.scratch)
    graph = nx.Graph()
    graph.add_nodes_from(int(s.get("scans")) for s in consensus)
    graph.add_edges_from(edges)
    members_of_family = {}
    for component in nx.connected_components(graph):
        members = sorted(component)
        members_of_family[f"MF{members[0]}"] = members

    n_ms2, n_usable = {}, {}
    for strain in strains:
        spectra, total = ref.strain_spectra(
            work / "features" / f"{strain['accession']}.npz")
        n_ms2[strain["accession"]] = total
        n_usable[strain["accession"]] = len(spectra)

    cached = {
        "strains": strains, "species": species, "domains_of": domains_of,
        "cluster_rows": cluster_rows, "strain_of_cluster": strain_of_cluster,
        "consensus": consensus, "members_of_family": members_of_family,
        "accessions": [s["accession"] for s in strains],
        "n_ms2": n_ms2, "n_usable": n_usable,
    }

    names = [a.path.removesuffix(".csv") for a in config.artifacts]
    header = f"{'variant':38}{'task':>6}  " + "".join(f"{n[:13]:>15}" for n in names)
    print(header)
    print("-" * len(header))
    passes, measured = Counter(), Counter()
    for label, mistakes in VARIANTS.items():
        if args.only and args.only not in label:
            continue
        out = args.scratch / label.replace(" ", "_").replace(":", "").replace("/", "")
        build(ref, args.task, out, mistakes, args.scratch, cached)
        result = score_task(config, args.truth, out)
        scored = {a.path: a for a in result.artifacts + result.soft_artifacts}
        row = f"{label:38}{result.task_score:6.2f}  "
        for spec in config.artifacts:
            artifact = scored[spec.path]
            row += f"{artifact.score:15.2f}"
            measured[spec.path] += 1
            if artifact.passed:
                passes[spec.path] += 1
        print(row, flush=True)

    total = max(measured.values()) if measured else 0
    print()
    print(f"artifact pass rate across the {total} variants "
          f"(PROJECT.md section 8 wants 0.1 to 0.8):")
    for artifact in config.artifacts:
        n = measured[artifact.path]
        rate = passes[artifact.path] / n if n else 0.0
        verdict = ("soft, out of the headline" if artifact.soft
                   else "ok" if 0.1 <= rate <= 0.8 else "OUTSIDE BAND")
        print(f"  {artifact.path:34} {rate:4.2f}  {verdict}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
