"""Measure task 02's difficulty by grading plausible wrong pipelines.

PROJECT.md §8: difficulty is measured, not designed, and an artifact is only
informative if its pass rate falls between roughly 0.1 and 0.8. Running frontier
models twice per artifact is the real calibration, but it is slow and costs
money, and it cannot be done before the task is built.

This is the cheap proxy. Each variant below is a mistake a competent agent
actually makes, applied to the same GECCO output the reference solution reads.
Grading them says which artifacts discriminate and which are free marks.

    python scripts/measure_difficulty.py --task tasks/02_salinispora_bgc \
        --truth tasks/02_salinispora_bgc/truth
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from npbench.grading import load_grading_config, score_task  # noqa: E402


def load_solution(task_dir: Path):
    spec = importlib.util.spec_from_file_location(
        "ref_solution", task_dir / "solution" / "make_truth.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# Each flag is one mistake. Names are the failure, not the fix.
MISTAKES = {
    "domains_as_multiset": "uses the raw `domains` field without deduplicating it",
    "pairwise_not_transitive": "treats each linked pair as a family, no connected components",
    "families_numbered_by_rank": "names families GCF_0001.. by size rank instead of by member",
    "keeps_thin_clusters": "skips the fewer-than-3-distinct-domains exclusion",
    "species_is_full_organism": "groups by the whole organism string, strain included",
    "length_off_by_one": "computes length as end - start",
    "splits_hybrid_types": "splits `NRP;Polyketide` into separate rows",
    "edge_rule_exact_boundary": "flags a contig edge only at base 1 or the last base",
    "domain_count_includes_duplicates": "reports the domain field's length, not its distinct size",
    "coordinates_zero_based": "re-bases GECCO's coordinates to 0-based half-open",
    "drops_unknown_type": "filters out clusters GECCO could not type",
    "drops_thin_from_inventory": "applies the family eligibility filter to every artifact",
}


_LENGTH_CACHE: dict[str, dict[str, int]] = {}


def build(task_dir: Path, work: Path, out: Path, mistakes: set[str]) -> None:
    ref = load_solution(task_dir)
    with (task_dir / "strains.tsv").open() as handle:
        strains = list(csv.DictReader(handle, delimiter="\t"))
    accessions = sorted(s["accession"] for s in strains)
    species = {
        s["accession"]: (
            s["organism"] if "species_is_full_organism" in mistakes
            else ref.species_of(s["organism"])
        )
        for s in strains
    }

    inventory, edge, domains, assembly_bp = [], [], {}, {}
    for accession in accessions:
        if accession not in _LENGTH_CACHE:
            _LENGTH_CACHE[accession] = ref.contig_lengths(
                task_dir / "data" / "genomes" / f"{accession}.fna.gz"
            )
        lengths = _LENGTH_CACHE[accession]
        assembly_bp[accession] = sum(lengths.values())
        for row in ref.read_clusters(work, accession):
            start, end = int(row["start"]), int(row["end"])
            sequence_id, cluster_id = row["sequence_id"], row["cluster_id"]
            raw = [d for d in row["domains"].split(";") if d]
            if "drops_unknown_type" in mistakes and row["type"] == "Unknown":
                continue
            if "drops_thin_from_inventory" in mistakes and len(set(raw)) < 3:
                continue
            if "coordinates_zero_based" in mistakes:
                start -= 1
            domains[cluster_id] = (
                raw if "domains_as_multiset" in mistakes else frozenset(raw)
            )
            types = (
                row["type"].split(";") if "splits_hybrid_types" in mistakes else [row["type"]]
            )
            for kind in types:
                inventory.append({
                    "cluster_id": cluster_id, "accession": accession,
                    "sequence_id": sequence_id, "start": start, "end": end,
                    "length_bp": end - start if "length_off_by_one" in mistakes
                    else end - start + 1,
                    "n_distinct_domains": len(raw)
                    if "domain_count_includes_duplicates" in mistakes else len(set(raw)),
                    "type": kind,
                })
            contig = lengths[sequence_id]
            if "edge_rule_exact_boundary" in mistakes:
                at_start, at_end = start == 1, end == contig
            else:
                at_start = start <= ref.EDGE_DISTANCE
                at_end = end >= contig - ref.EDGE_DISTANCE + 1
            if at_start or at_end:
                edge.append({
                    "cluster_id": cluster_id, "accession": accession,
                    "sequence_id": sequence_id, "contig_length_bp": contig,
                    "boundary": "both" if at_start and at_end
                    else ("start" if at_start else "end"),
                })

    inventory.sort(key=lambda r: r["cluster_id"])
    edge.sort(key=lambda r: r["cluster_id"])
    families = _families(domains, mistakes)

    size = Counter(families.values())
    accession_of = {r["cluster_id"]: r["accession"] for r in inventory}
    members: dict[str, list[str]] = {}
    for cluster, family in families.items():
        members.setdefault(family, []).append(cluster)
    per_type = Counter(r["type"] for r in inventory)
    genomes_per_type: dict[str, set[str]] = {}
    for row in inventory:
        genomes_per_type.setdefault(row["type"], set()).add(row["accession"])
    per_genome = Counter(r["accession"] for r in inventory)
    per_species_type = Counter((species[r["accession"]], r["type"]) for r in inventory)
    edge_per_genome = Counter(r["accession"] for r in edge)
    cluster_bp = Counter()
    for row in inventory:
        cluster_bp[row["accession"]] += row["length_bp"]
    species_per_type: dict[str, set[str]] = {}
    lengths_per_type: dict[str, list[int]] = {}
    for row in inventory:
        species_per_type.setdefault(row["type"], set()).add(species[row["accession"]])
        lengths_per_type.setdefault(row["type"], []).append(row["length_bp"])

    w = ref.write_csv
    w(out / "bgc_inventory.csv",
      ["cluster_id", "accession", "sequence_id", "start", "end", "length_bp",
       "n_distinct_domains", "type"], inventory)
    w(out / "contig_edge_clusters.csv",
      ["cluster_id", "accession", "sequence_id", "contig_length_bp", "boundary"], edge)
    w(out / "clusters_per_genome.csv",
      ["accession", "n_clusters", "n_edge_clusters", "assembly_bp", "pct_bp_in_clusters"],
      [{"accession": a, "n_clusters": per_genome.get(a, 0),
        "n_edge_clusters": edge_per_genome.get(a, 0), "assembly_bp": assembly_bp[a],
        "pct_bp_in_clusters": round(100 * cluster_bp.get(a, 0) / assembly_bp[a], 4)}
       for a in accessions])
    w(out / "type_counts.csv",
      ["type", "n_clusters", "n_genomes", "n_species", "median_length_bp"],
      [{"type": t, "n_clusters": per_type[t], "n_genomes": len(genomes_per_type[t]),
        "n_species": len(species_per_type[t]),
        "median_length_bp": ref._median(lengths_per_type[t])}
       for t in sorted(per_type)])
    w(out / "species_type_counts.csv", ["species", "type", "n_clusters"],
      [{"species": s, "type": t, "n_clusters": n} for (s, t), n in sorted(per_species_type.items())])
    w(out / "gene_cluster_families.csv", ["cluster_id", "family_id", "family_size"],
      [{"cluster_id": c, "family_id": f, "family_size": size[f]}
       for c, f in sorted(families.items())])
    w(out / "gcf_summary.csv", ["family_id", "n_clusters", "n_genomes", "n_species"],
      [{"family_id": f, "n_clusters": len(g),
        "n_genomes": len({accession_of[c] for c in g}),
        "n_species": len({species[accession_of[c]] for c in g})}
       for f, g in sorted(members.items())])


def _similarity(x, y) -> float:
    """Jaccard. Over multisets when the variant kept the duplicates."""
    if isinstance(x, frozenset):
        return len(x & y) / len(x | y) if (x | y) else 0.0
    a, b = Counter(x), Counter(y)
    union = sum((a | b).values())
    return sum((a & b).values()) / union if union else 0.0


def _families(domains: dict, mistakes: set[str]) -> dict[str, str]:
    counts = {c: (len(set(d)) if not isinstance(d, frozenset) else len(d))
              for c, d in domains.items()}
    eligible = sorted(
        domains if "keeps_thin_clusters" in mistakes
        else {c: d for c, d in domains.items() if counts[c] >= 3}
    )
    jaccard = {}
    for i, a in enumerate(eligible):
        for b in eligible[i + 1:]:
            if _similarity(domains[a], domains[b]) >= 0.50:
                jaccard[(a, b)] = True

    if "pairwise_not_transitive" in mistakes:
        # Name each cluster after its best single partner: no components at all.
        best = {c: c for c in eligible}
        for a, b in jaccard:
            best[b] = min(best[b], a)
        assignment = best
    else:
        parent = {c: c for c in eligible}

        def find(c):
            while parent[c] != c:
                parent[c] = parent[parent[c]]
                c = parent[c]
            return c

        for a, b in jaccard:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[max(ra, rb)] = min(ra, rb)
        assignment = {c: find(c) for c in eligible}

    if "families_numbered_by_rank" in mistakes:
        counts = Counter(assignment.values())
        order = sorted(counts, key=lambda f: (-counts[f], f))
        rename = {f: f"GCF_{i + 1:04d}" for i, f in enumerate(order)}
        assignment = {c: rename[f] for c, f in assignment.items()}
    return assignment


VARIANTS = {
    "reference (no mistakes)": set(),
    "one slip: domains as multiset": {"domains_as_multiset"},
    "one slip: families numbered by rank": {"families_numbered_by_rank"},
    "one slip: no transitive closure": {"pairwise_not_transitive"},
    "one slip: species keeps the strain": {"species_is_full_organism"},
    "one slip: length off by one": {"length_off_by_one"},
    "one slip: edge rule too strict": {"edge_rule_exact_boundary"},
    "one slip: hybrid types split": {"splits_hybrid_types"},
    "one slip: thin clusters kept": {"keeps_thin_clusters"},
    "one slip: domain count counts dupes": {"domain_count_includes_duplicates"},
    "one slip: coordinates 0-based": {"coordinates_zero_based"},
    "one slip: drops Unknown-type clusters": {"drops_unknown_type"},
    "one slip: thin clusters dropped everywhere": {"drops_thin_from_inventory"},
    "hasty agent (3 slips)": {"domains_as_multiset", "species_is_full_organism",
                              "edge_rule_exact_boundary"},
    "careless agent (5 slips)": {"domains_as_multiset", "families_numbered_by_rank",
                                 "species_is_full_organism", "length_off_by_one",
                                 "splits_hybrid_types"},
}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--task", type=Path, required=True)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--work", type=Path, default=None)
    parser.add_argument("--scratch", type=Path, required=True)
    args = parser.parse_args()
    work = args.work or args.task / "work"
    config = load_grading_config(args.task / "grading.yaml")

    names = [a.path.removesuffix(".csv") for a in config.artifacts]
    header = f"{'variant':34}{'task':>6}  " + "".join(f"{n[:13]:>15}" for n in names)
    print(header)
    print("-" * len(header))
    passes = Counter()
    for label, mistakes in VARIANTS.items():
        out = args.scratch / label.replace(" ", "_").replace(":", "")
        build(args.task, work, out, mistakes)
        result = score_task(config, args.truth, out)
        row = f"{label:34}{result.task_score:6.2f}  "
        for artifact in result.artifacts:
            row += f"{artifact.score:15.2f}"
            if artifact.passed:
                passes[artifact.path] += 1
        print(row)

    print()
    print("artifact pass rate across the", len(VARIANTS), "variants "
          "(PROJECT.md section 8 wants 0.1 to 0.8):")
    for artifact in config.artifacts:
        rate = passes[artifact.path] / len(VARIANTS)
        flag = "ok" if 0.1 <= rate <= 0.8 else "OUTSIDE BAND"
        print(f"  {artifact.path:32} {rate:4.2f}  {flag}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
