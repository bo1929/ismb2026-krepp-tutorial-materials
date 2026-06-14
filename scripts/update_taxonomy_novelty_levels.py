#!/usr/bin/env python3
"""
Assign novelty_level fields from taxonomic LCA between queries (query_taxonomy.tsv)
and references (NCBI lineages). Updates data/query_info.tsv and data/reference_info.tsv.

Usage:
  python3 scripts/update_taxonomy_novelty_levels.py
  python3 scripts/update_taxonomy_novelty_levels.py --dry-run
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
QUERY_TAXONOMY = ROOT / "data" / "query_taxonomy.tsv"
QUERY_INFO = ROOT / "data" / "query_info.tsv"
REFERENCE_INFO = ROOT / "data" / "reference_info.tsv"
REF_TAXONOMY_CACHE = ROOT / "data" / "reference_taxonomy.tsv"

ENTREZ = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

# NCBI ranks we emit as novelty_level (deepest shared rank name).
RANK_PRIORITY = {
    "species": 7,
    "subspecies": 7,
    "strain": 7,
    "genus": 6,
    "family": 5,
    "order": 4,
    "class": 3,
    "phylum": 2,
    "kingdom": 1,
    "domain": 1,
    "superkingdom": 0,
}


def entrez_sleep():
    time.sleep(0.34)


def parse_taxpath(taxpath: str) -> list[str]:
    return [t for t in taxpath.split("|") if t.strip()]


def load_query_taxpaths(path: Path) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            taxid, taxpath = parts[0], parts[2]
            if taxid.isdigit():
                out[taxid] = parse_taxpath(taxpath)
    return out


def fetch_lineage(tax_id: str, path_cache: dict[str, list[str]], rank_cache: dict[str, str]) -> list[str]:
    if tax_id in path_cache:
        return path_cache[tax_id]
    url = f"{ENTREZ}efetch.fcgi?db=taxonomy&id={tax_id}&retmode=xml"
    entrez_sleep()
    try:
        with urllib.request.urlopen(url, timeout=120) as r:
            root = ET.fromstring(r.read())
    except urllib.error.HTTPError as exc:
        sys.stderr.write(f"warning: taxonomy fetch failed for {tax_id}: {exc}\n")
        path_cache[tax_id] = [tax_id]
        rank_cache[tax_id] = "species"
        return path_cache[tax_id]
    taxon = root.find(".//Taxon")
    if taxon is None:
        path_cache[tax_id] = [tax_id]
        rank_cache[tax_id] = "species"
        return path_cache[tax_id]
    ids: list[str] = []
    for node in taxon.findall("LineageEx/Taxon"):
        tid = (node.findtext("TaxId") or "").strip()
        rk = (node.findtext("Rank") or "").strip().lower()
        if tid:
            ids.append(tid)
            rank_cache[tid] = rk
    tid = (taxon.findtext("TaxId") or tax_id).strip()
    rk = (taxon.findtext("Rank") or "species").strip().lower()
    if not ids or ids[-1] != tid:
        ids.append(tid)
    rank_cache[tid] = rk
    path_cache[tax_id] = ids
    return ids


def load_ref_taxonomy_cache(path: Path) -> dict[str, list[str]]:
    if not path.is_file():
        return {}
    out: dict[str, list[str]] = {}
    with path.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2 or not parts[0].isdigit():
                continue
            out[parts[0]] = parse_taxpath(parts[1])
    return out


def write_ref_taxonomy_cache(path: Path, cache: dict[str, list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        fh.write("taxid\ttaxpath\n")
        for tid in sorted(cache, key=int):
            fh.write(f"{tid}\t{'|'.join(cache[tid])}\n")


def lca_taxid(path_a: list[str], path_b: list[str]) -> str | None:
    """Deepest taxon shared between paths (Kraken paths may omit NCBI root nodes)."""
    in_b = set(path_b)
    for tid in reversed(path_a):
        if tid in in_b:
            return tid
    for tid in reversed(path_b):
        if tid in set(path_a):
            return tid
    return None


def novelty_level_from_lca(lca_id: str | None, rank_cache: dict[str, str]) -> str:
    if not lca_id:
        return "kingdom"
    rk = rank_cache.get(lca_id, "")
    if rk in ("domain", "superkingdom"):
        return "kingdom"
    if rk in RANK_PRIORITY:
        if rk in ("subspecies", "strain"):
            return "species"
        return rk
    return "kingdom"


def lca_depth(path_a: list[str], path_b: list[str], lca: str | None) -> int:
    if not lca:
        return 0
    depth = 0
    if lca in path_a:
        depth = max(depth, path_a.index(lca) + 1)
    if lca in path_b:
        depth = max(depth, path_b.index(lca) + 1)
    return depth


def best_match(
    source_id: str,
    source_path: list[str],
    targets: dict[str, list[str]],
    rank_cache: dict[str, str],
) -> tuple[str, str, str]:
    """Return (target_id, novelty_level, lca_taxid) for closest target."""
    best_tid = ""
    best_level = "kingdom"
    best_lca = ""
    best_depth = -1
    best_prio = -1
    for tid, tpath in targets.items():
        lca = lca_taxid(source_path, tpath)
        depth = lca_depth(source_path, tpath, lca)
        level = novelty_level_from_lca(lca, rank_cache)
        prio = RANK_PRIORITY.get(level, 0)
        if depth > best_depth or (depth == best_depth and prio > best_prio):
            best_depth = depth
            best_prio = prio
            best_tid = tid
            best_level = level
            best_lca = lca or ""
    return best_tid, best_level, best_lca


def read_reference_info(path: Path) -> tuple[list[str], list[dict]]:
    with path.open(encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        rows = []
        for parts in reader:
            if not parts or not parts[0].strip():
                continue
            rows.append({
                "i": parts[0],
                "label": parts[1],
                "novelty_level": parts[2],
                "profile_taxid": parts[3] if len(parts) > 3 else "",
                "profile_genus": parts[4] if len(parts) > 4 else "",
                "ref_accession": parts[5] if len(parts) > 5 else "",
                "ref_taxid": parts[6] if len(parts) > 6 else "",
                "ref_organism": parts[7] if len(parts) > 7 else "",
                "ref_class_taxid": parts[8] if len(parts) > 8 else "",
                "ref_class_name": parts[9] if len(parts) > 9 else "",
                "query_accession": parts[10] if len(parts) > 10 else "",
            })
    return header, rows


def read_query_info(path: Path) -> tuple[list[str], list[dict]]:
    with path.open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []
        rows = [dict(r) for r in reader]
    return list(header), rows


def write_reference_info(path: Path, header: list[str], rows: list[dict]) -> None:
    out_header = [
        "i", "label", "novelty_level", "profile_taxid", "profile_genus",
        "ref_accession", "ref_taxid", "ref_organism",
        "ref_class_taxid", "ref_class_name", "query_accession",
    ]
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(out_header)
        for r in rows:
            w.writerow([
                r["i"], r["label"], r["novelty_level"], r["profile_taxid"], r["profile_genus"],
                r["ref_accession"], r["ref_taxid"], r["ref_organism"],
                r["ref_class_taxid"], r["ref_class_name"], r["query_accession"],
            ])


def write_query_info(path: Path, header: list[str], rows: list[dict]) -> None:
    out_header = ["taxid", "taxon", "novelty_level", "accession", "organism"]
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=out_header, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in out_header})


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    query_paths = load_query_taxpaths(QUERY_TAXONOMY)
    if not query_paths:
        sys.stderr.write(f"No query paths in {QUERY_TAXONOMY}\n")
        return 1

    _, ref_rows = read_reference_info(REFERENCE_INFO)
    ref_taxids = sorted({r["ref_taxid"] for r in ref_rows if r["ref_taxid"].isdigit()}, key=int)

    ref_paths = load_ref_taxonomy_cache(REF_TAXONOMY_CACHE)
    rank_cache: dict[str, str] = {}
    missing = [t for t in ref_taxids if t not in ref_paths]
    if missing:
        sys.stderr.write(f"Fetching {len(missing)} reference lineages from NCBI...\n")
    for tid in ref_taxids:
        fetch_lineage(tid, ref_paths, rank_cache)

    rank_fetch_buf: dict[str, list[str]] = {}
    for path in query_paths.values():
        for tid in path:
            if tid not in rank_cache:
                fetch_lineage(tid, rank_fetch_buf, rank_cache)

    if not args.dry_run:
        write_ref_taxonomy_cache(REF_TAXONOMY_CACHE, ref_paths)

    qi_header, qi_rows = read_query_info(QUERY_INFO)
    query_targets = {tid: path for tid, path in query_paths.items()}

    ref_targets = {tid: ref_paths[tid] for tid in ref_taxids if tid in ref_paths}

    for row in qi_rows:
        qtid = (row.get("taxid") or "").strip()
        qpath = query_paths.get(qtid, [])
        if not qpath:
            continue
        _, level, _ = best_match(qtid, qpath, ref_targets, rank_cache)
        row["novelty_level"] = level

    for row in ref_rows:
        rtid = row["ref_taxid"].strip()
        rpath = ref_paths.get(rtid, [])
        if not rpath:
            continue
        _, level, _ = best_match(rtid, rpath, query_targets, rank_cache)
        row["novelty_level"] = level

    sys.stderr.write("Query novelty levels:\n")
    for row in sorted(qi_rows, key=lambda r: r.get("taxid", "")):
        sys.stderr.write(
            f"  {row.get('taxid')}\t{row.get('organism', '')[:40]}\t{row.get('novelty_level')}\n"
        )

    sys.stderr.write("\nReference novelty levels:\n")
    for row in ref_rows:
        sys.stderr.write(f"  {row['i']}\t{row['label'][:28]}\t{row['novelty_level']}\n")

    if args.dry_run:
        sys.stderr.write("\n(dry run; files not written)\n")
        return 0

    write_query_info(QUERY_INFO, qi_header, qi_rows)
    write_reference_info(REFERENCE_INFO, [], ref_rows)
    sys.stderr.write(f"\nWrote {QUERY_INFO} and {REFERENCE_INFO}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
