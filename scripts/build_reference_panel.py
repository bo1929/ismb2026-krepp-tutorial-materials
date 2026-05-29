#!/usr/bin/env python3
"""
Build 40-genome reference panel per reference_set plan:
  20 profile-linked + 10 genus-pool extras + 10 class-disjoint outgroups.

Resolves accessions via NCBI Entrez (esearch/esummary/efetch). Downloads FASTA
via NCBI Datasets HTTP API (same as fetch_genomes.py).

Usage:
  python3 scripts/build_reference_panel.py
  python3 scripts/build_reference_panel.py --relabel-only
"""

from __future__ import annotations

import argparse
import csv
import re
import gzip
import json
import random
import shutil
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

ENTREZ = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
NCBI_GENOME_DOWNLOAD_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/download"

# Fixed partition (20 profile species).
# NCBI has only one assembly for several species (and some monotypic genera), so
# strict 15/3/2 is infeasible. Implemented: 11 same-species, 6 same-genus,
# 3 cross-genus (still one cross row is Archaeon: Hyperthermus).
CROSS_GENUS = [
    ("54248", "Hyperthermus"),
    ("186192", "Marinithermus"),
    ("80679", "Thioflavicoccus"),
]
SAME_GENUS = [
    ("1263979", "Candidatus Endolissoclinum"),
    ("1298608", "Psychrobacter"),
    ("971279", "Palaeococcus"),
    ("133539", "Nitrosococcus"),
    ("1335757", "Spiribacter"),
    ("387093", "Sulfurovum"),
]
SAME_SPECIES = [
    "28108",
    "1229",
    "314275",
    "39960",
    "936476",
    "1335746",
    "80852",
    "40269",
    "55601",
    "187493",
    "467094",
]

RNG_SEED = 42
EXTRA_COUNT = 10
OUTGROUP_COUNT = 10

# Canonical column order for all panel TSV outputs (header + every row).
DESIGN_FIELDNAMES = [
    "short_id",
    "genome_id",
    "role",
    "profile_taxid",
    "profile_genus",
    "ref_accession",
    "ref_taxid",
    "ref_organism",
    "ref_class_taxid",
    "ref_class_name",
    "query_accession",
]


def entrez_sleep():
    time.sleep(0.34)


def http_json(url: str):
    entrez_sleep()
    with urllib.request.urlopen(url, timeout=180) as r:
        return json.loads(r.read().decode())


def http_bytes(url: str, data: bytes | None = None, headers: dict | None = None):
    entrez_sleep()
    req = urllib.request.Request(url, data=data, headers=headers or {}, method="POST" if data else "GET")
    with urllib.request.urlopen(req, timeout=300) as r:
        return r.read()


def taxonomy_primary_scientific_name(tax_id: str, cache: dict[str, str]) -> str:
    if tax_id in cache:
        return cache[tax_id]
    url = f"{ENTREZ}efetch.fcgi?db=taxonomy&id={tax_id}&retmode=xml"
    entrez_sleep()
    with urllib.request.urlopen(url, timeout=120) as r:
        root = ET.fromstring(r.read())
    name = ""
    for taxon in root.findall(".//Taxon"):
        if taxon.findtext("Rank") == "species":
            name = taxon.findtext("ScientificName") or ""
            break
    if not name:
        taxon = root.find(".//Taxon")
        if taxon is not None:
            name = taxon.findtext("ScientificName") or ""
    out = name.strip() or f"taxid_{tax_id}"
    cache[tax_id] = out
    return out


def taxonomy_class(tax_id: int) -> tuple[int | None, str | None]:
    url = f"{ENTREZ}efetch.fcgi?db=taxonomy&id={tax_id}&retmode=xml"
    entrez_sleep()
    with urllib.request.urlopen(url, timeout=120) as r:
        xml = r.read()
    root = ET.fromstring(xml)
    for taxon in root.findall(".//Taxon"):
        if taxon.findtext("Rank") == "class":
            tid = taxon.findtext("TaxId")
            name = taxon.findtext("ScientificName")
            return (int(tid), name)
    return (None, None)


def load_query_accessions(path: Path) -> dict[str, str]:
    out = {}
    with open(path, newline="") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for rec in rdr:
            tid = (rec.get("taxid") or "").strip()
            acc = (rec.get("assembly_accession") or "").strip()
            if tid.isdigit() and acc:
                out[tid] = acc
    return out


def assembly_esummary(uids: list[str]) -> dict:
    if not uids:
        return {}
    idstr = ",".join(uids)
    url = f"{ENTREZ}esummary.fcgi?db=assembly&id={idstr}&retmode=json"
    return http_json(url)


def assembly_search(term: str, retmax: int = 100) -> list[str]:
    q = urllib.parse.quote(term)
    url = f"{ENTREZ}esearch.fcgi?db=assembly&term={q}&retmax={retmax}&retmode=json&sort=significance"
    jr = http_json(url)
    return jr.get("esearchresult", {}).get("idlist", [])


def score_doc(doc: dict) -> int:
    s = 0
    if doc.get("ftppath_refseq"):
        s += 40
    ref = str(doc.get("refseq_category", "")).lower()
    if "reference" in ref:
        s += 50
    st = doc.get("assemblystatus", "")
    if st == "Complete Genome":
        s += 30
    elif st == "Chromosome":
        s += 20
    elif st == "Scaffold":
        s += 10
    ex = doc.get("exclfromrefseq", [])
    if ex:
        s -= 25
    return s


def species_taxid_for_node(tax_id: str, cache: dict[str, str]) -> str:
    if tax_id in cache:
        return cache[tax_id]
    url = f"{ENTREZ}efetch.fcgi?db=taxonomy&id={tax_id}&retmode=xml"
    entrez_sleep()
    with urllib.request.urlopen(url, timeout=120) as r:
        xml = r.read()
    root = ET.fromstring(xml)
    # Prefer explicit species rank in lineage
    sp = None
    for taxon in root.findall(".//Taxon"):
        if taxon.findtext("Rank") == "species":
            sp = taxon.findtext("TaxId")
            break
    if sp is None:
        for taxon in root.findall(".//Taxon"):
            if taxon.findtext("Rank") in ("species", "no rank") and taxon.findtext("ScientificName"):
                cand = taxon.findtext("TaxId")
                if cand:
                    sp = cand
                    break
    out = sp or tax_id
    cache[tax_id] = out
    return out


def pick_same_species_alt(
    tax_id: str,
    exclude_acc: str,
    profile_taxa: set[str],
    lineage_cache: dict[str, str],
    forbidden_acc: set[str],
) -> tuple[str, str, str]:
    passes = (
        (True, True),
        (False, True),
        (True, False),
        (False, False),
    )
    best_global = None
    ids = assembly_search(f"txid{tax_id}[Organism]", retmax=80)
    if not ids:
        raise RuntimeError(f"no assemblies for taxid {tax_id}")
    summ = assembly_esummary(ids)
    result = summ.get("result", {})
    for req_ftp, need_good_level in passes:
        best = None
        for uid in result.get("uids", []):
            doc = result[uid]
            acc = doc.get("assemblyaccession")
            if not acc or acc == exclude_acc or acc in forbidden_acc:
                continue
            if req_ftp and not doc.get("ftppath_refseq"):
                continue
            atid = str(doc.get("taxid", ""))
            stid = species_taxid_for_node(atid, lineage_cache)
            if stid != tax_id:
                continue
            st = doc.get("assemblystatus", "")
            if need_good_level and st not in ("Complete Genome", "Chromosome", "Scaffold"):
                continue
            rec = (score_doc(doc), acc, doc.get("organism_name", ""))
            if best is None or rec[0] > best[0]:
                best = rec
        if best:
            best_global = best
            break
    if not best_global:
        raise RuntimeError(f"no alternate assembly for taxid {tax_id} != {exclude_acc}")
    return best_global[1], tax_id, best_global[2]


def taxid_from_assembly_accession(acc: str) -> str:
    ids = assembly_search(acc, retmax=5)
    if not ids:
        raise RuntimeError(f"assembly not found {acc}")
    summ = assembly_esummary(ids[:1])
    uid = summ["result"]["uids"][0]
    doc = summ["result"][uid]
    return str(doc.get("taxid", ""))


def pick_genus_sibling(
    genus_term: str,
    profile_species_taxids: set[str],
    lineage_cache: dict[str, str],
    forbidden_acc: set[str],
    prefer_complete: bool = True,
) -> tuple[str, str, str]:
    def scan(retmax: int, pref: bool, require_ftp_refseq: bool) -> tuple[str, str, str] | None:
        ids = assembly_search(f"{genus_term}[Organism]", retmax=retmax)
        summ = assembly_esummary(ids)
        result = summ.get("result", {})
        best = None
        for uid in result.get("uids", []):
            doc = result[uid]
            acc = doc.get("assemblyaccession")
            if not acc or acc in forbidden_acc:
                continue
            if require_ftp_refseq and not doc.get("ftppath_refseq"):
                continue
            atid = str(doc.get("taxid", ""))
            stid = species_taxid_for_node(atid, lineage_cache)
            if stid in profile_species_taxids:
                continue
            if pref and doc.get("assemblystatus") not in (
                "Complete Genome", "Chromosome", "Scaffold",
            ):
                continue
            rec = (score_doc(doc), acc, stid, doc.get("organism_name", ""))
            if best is None or rec[0] > best[0]:
                best = rec
        if not best:
            return None
        return best[1], best[2], best[3]

    out = scan(120, prefer_complete, True)
    if out is None:
        out = scan(250, False, True)
    if out is None:
        out = scan(250, False, False)
    if out is None:
        raise RuntimeError(f"no genus sibling for {genus_term}")
    return out


def pick_cross_genus_hyperthermus() -> tuple[str, str, str]:
    acc = "GCF_000007085.1"
    tid = taxid_from_assembly_accession(acc)
    return acc, tid, "Pyrodictium occultum DSM 2709"


def pick_cross_genus_marinithermus() -> tuple[str, str, str]:
    acc = "GCF_000006785.2"
    tid = taxid_from_assembly_accession(acc)
    return acc, tid, "Thermus thermophilus HB8"


def pick_cross_genus_thioflavicoccus() -> tuple[str, str, str]:
    acc = "GCF_900119735.1"
    tid = taxid_from_assembly_accession(acc)
    return acc, tid, "Allochromatium vinosum DSM 180"


def pick_cross_genus_bacterium(_profile_genus: str) -> tuple[str, str, str]:
    # Escherichia coli str. K-12 - different genus; common lab ref
    acc = "GCF_000005845.2"
    tid = taxid_from_assembly_accession(acc)
    return acc, tid, "Escherichia coli str. K-12 substr. MG1655"


def download_via_ncbi_ftp_https(accession: str, out_fa: Path) -> bool:
    ids = assembly_search(accession, retmax=5)
    if not ids:
        return False
    summ = assembly_esummary(ids[:1])
    uid = summ["result"]["uids"][0]
    doc = summ["result"][uid]
    ft = doc.get("ftppath_refseq") or doc.get("ftppath_genbank")
    if not ft:
        return False
    base = ft.replace("ftp://ftp.ncbi.nlm.nih.gov", "https://ftp.ncbi.nlm.nih.gov").rstrip("/")
    leaf = base.split("/")[-1]
    gz_url = f"{base}/{leaf}_genomic.fna.gz"
    entrez_sleep()
    try:
        with urllib.request.urlopen(gz_url, timeout=300) as resp:
            raw_gz = resp.read()
        raw = gzip.decompress(raw_gz)
    except Exception:
        return False
    with open(out_fa, "wb") as out:
        out.write(raw)
    return True


def download_fasta(accession: str, out_fa: Path) -> None:
    if out_fa.exists() and out_fa.stat().st_size > 0:
        return
    body = json.dumps({
        "accessions": [accession],
        "include_annotation_type": ["GENOME_FASTA"],
    }).encode("utf-8")
    last_err = None
    for attempt in range(8):
        req = urllib.request.Request(
            NCBI_GENOME_DOWNLOAD_URL,
            data=body,
            headers={"Content-Type": "application/json"},
            method="POST",
        )
        try:
            with tempfile.TemporaryDirectory() as td:
                zp = Path(td) / "g.zip"
                with urllib.request.urlopen(req, timeout=300) as resp:
                    with open(zp, "wb") as out:
                        shutil.copyfileobj(resp, out)
                head = zp.read_bytes()[:4]
                if head[:2] != b"PK":
                    raise zipfile.BadZipFile("not a zip response")
                with zipfile.ZipFile(zp) as zf:
                    for name in zf.namelist():
                        if name.endswith(".fna"):
                            with zf.open(name) as src, open(out_fa, "wb") as dst:
                                shutil.copyfileobj(src, dst)
                            return
                raise zipfile.BadZipFile(f"no fna in zip for {accession}")
        except (zipfile.BadZipFile, urllib.error.HTTPError, OSError) as e:
            last_err = e
            time.sleep(3.0 * (attempt + 1))
    if download_via_ncbi_ftp_https(accession, out_fa):
        return
    raise RuntimeError(f"download failed for {accession}: {last_err}") from last_err


def fill_missing_ref_organism(rows: list[dict], cache: dict[str, str]) -> None:
    for r in rows:
        if (r.get("ref_organism") or "").strip():
            continue
        tid = str(r.get("ref_taxid") or "").strip()
        if not tid:
            continue
        r["ref_organism"] = taxonomy_primary_scientific_name(tid, cache)


def species_slug_base(organism_label: str) -> str:
    words = organism_label.strip().split()
    if not words:
        return "unknown"
    if words[0] == "Candidatus" and len(words) >= 3:
        tokens = words[:3]
    elif len(words) >= 2:
        tokens = words[:2]
    else:
        tokens = words[:1]
    raw = "_".join(tokens)
    raw = re.sub(r"[^A-Za-z0-9_]+", "_", raw)
    raw = re.sub(r"_+", "_", raw).strip("_")
    return raw or "unknown"


def assign_genome_ids(rows: list[dict]) -> None:
    counts: dict[str, int] = {}
    used_final: set[str] = set()
    for r in rows:
        base = species_slug_base(r.get("ref_organism") or "")
        counts[base] = counts.get(base, 0) + 1
        n = counts[base]
        gid = base if n == 1 else f"{base}_{n}"
        if len(gid) > 120:
            gid = gid[:120].rstrip("_")
        stem = gid
        bump = 0
        while gid in used_final:
            bump += 1
            suffix = f"_{bump}"
            head = stem[: max(1, 120 - len(suffix))].rstrip("_")
            gid = (head + suffix)[:120]
        used_final.add(gid)
        r["genome_id"] = gid


def write_panel_tables(
    rows: list[dict],
    design_path: Path,
    manifest_path: Path,
    input_map_path: Path,
) -> None:
    design_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    input_map_path.parent.mkdir(parents=True, exist_ok=True)
    for path in (design_path, manifest_path):
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(
                fh,
                fieldnames=DESIGN_FIELDNAMES,
                delimiter="\t",
                extrasaction="ignore",
                lineterminator="\n",
            )
            w.writeheader()
            for r in rows:
                w.writerow({k: r.get(k, "") for k in DESIGN_FIELDNAMES})
    with open(input_map_path, "w", newline="", encoding="utf-8") as fh:
        for r in rows:
            gid = r["genome_id"]
            fh.write(f"{gid}\t{gid}.fna\n")


def load_panel_rows_tsv(path: Path) -> list[dict]:
    rows: list[dict] = []
    with open(path, newline="", encoding="utf-8") as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for raw in rdr:
            row = {fn: "" for fn in DESIGN_FIELDNAMES}
            for key in rdr.fieldnames or []:
                if key in DESIGN_FIELDNAMES:
                    val = raw.get(key)
                    row[key] = val.strip() if isinstance(val, str) else (val or "")
            rows.append(row)
    return rows


def relabel_panel_fastas(
    rows: list[dict],
    outdir: Path,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    for r in rows:
        sid = r.get("short_id") or ""
        gid = r.get("genome_id") or ""
        if not sid or not gid:
            continue
        old_p = outdir / f"{sid}.fna"
        new_p = outdir / f"{gid}.fna"
        if old_p == new_p:
            continue
        if new_p.exists() and not old_p.exists():
            continue
        if old_p.exists():
            if new_p.exists():
                sys.stderr.write(
                    f"error: both exist, refusing overwrite: {old_p} and {new_p}\n"
                )
                sys.exit(4)
            shutil.move(str(old_p), str(new_p))


def genus_name_for_taxid(tax_id: str) -> str:
    url = f"{ENTREZ}efetch.fcgi?db=taxonomy&id={tax_id}&retmode=xml"
    entrez_sleep()
    with urllib.request.urlopen(url, timeout=120) as r:
        xml = r.read()
    root = ET.fromstring(xml)
    for taxon in root.findall(".//Taxon"):
        if taxon.findtext("Rank") == "genus":
            return taxon.findtext("ScientificName") or ""
    return ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parent.parent)
    ap.add_argument("--profile-accessions", type=Path, default=None)
    ap.add_argument("--outdir", type=Path, default=None)
    ap.add_argument(
        "--relabel-only",
        action="store_true",
        help="Read design/manifest, assign genome_id names, rename FASTAs, rewrite TSVs (no download).",
    )
    args = ap.parse_args()
    root = args.repo_root
    acc_tsv = args.profile_accessions or (root / "profile_species_accessions.tsv")
    outdir = args.outdir or (root / "reference_panel")
    manifest_path = outdir / "reference_panel_manifest.tsv"
    design_path = root / "data" / "reference_set_design.tsv"
    input_map_path = outdir / "input_map.tsv"

    if args.relabel_only:
        src = design_path if design_path.exists() else manifest_path
        if not src.exists():
            sys.stderr.write(f"error: need {design_path} or {manifest_path}\n")
            sys.exit(2)
        rows = load_panel_rows_tsv(src)
        tax_cache: dict[str, str] = {}
        fill_missing_ref_organism(rows, tax_cache)
        assign_genome_ids(rows)
        relabel_panel_fastas(rows, outdir)
        write_panel_tables(rows, design_path, manifest_path, input_map_path)
        sys.stderr.write(
            f"Relabeled {len(rows)} genomes; updated {design_path}, {manifest_path}, {input_map_path}\n"
        )
        return

    query_map = load_query_accessions(acc_tsv)
    profile_taxa = set(query_map.keys())

    same_genus_taxa = {p[0] for p in SAME_GENUS}
    cross_taxa = {p[0] for p in CROSS_GENUS}
    same_species_set = set(SAME_SPECIES)
    assert len(profile_taxa) == 20
    assert same_genus_taxa | cross_taxa | same_species_set == profile_taxa

    class_by_species: dict[str, tuple[int | None, str | None]] = {}
    query_classes: set[int] = set()
    for tid in profile_taxa:
        c_id, c_name = taxonomy_class(int(tid))
        class_by_species[tid] = (c_id, c_name)
        if c_id is not None:
            query_classes.add(c_id)

    lineage_cache: dict[str, str] = {}

    rows: list[dict] = []
    rid = 0
    used_acc: set[str] = set(query_map.values())

    # --- profile linked ---
    for ptid, pgen in CROSS_GENUS:
        exclude = query_map[ptid]
        if ptid == "54248":
            acc, rtid, org = pick_cross_genus_hyperthermus()
        elif ptid == "186192":
            acc, rtid, org = pick_cross_genus_marinithermus()
        elif ptid == "80679":
            acc, rtid, org = pick_cross_genus_thioflavicoccus()
        else:
            acc, rtid, org = pick_cross_genus_bacterium(pgen)
        if acc in used_acc:
            raise RuntimeError(f"cross-genus acc already used: {acc}")
        used_acc.add(acc)
        rid += 1
        short = f"ref_{rid:02d}"
        cls_id, cls_name = taxonomy_class(int(rtid))
        rows.append({
            "short_id": short,
            "role": "profile_cross_genus",
            "profile_taxid": ptid,
            "profile_genus": pgen,
            "ref_accession": acc,
            "ref_taxid": rtid,
            "ref_organism": org,
            "ref_class_taxid": cls_id or "",
            "ref_class_name": cls_name or "",
            "query_accession": exclude,
        })

    for ptid, pgen in SAME_GENUS:
        exclude = query_map[ptid]
        if ptid == "971279":
            term = "Palaeococcus"
        elif ptid in ("133539", "1335757", "387093"):
            term = pgen
        else:
            term = pgen.replace("Candidatus ", "").split()[0]
            if ptid == "1263979":
                term = "Endolissoclinum"
        acc, rtid, org = pick_genus_sibling(term, profile_taxa, lineage_cache, used_acc)
        used_acc.add(acc)
        rid += 1
        short = f"ref_{rid:02d}"
        cls_id, cls_name = taxonomy_class(int(rtid))
        rows.append({
            "short_id": short,
            "role": "profile_same_genus",
            "profile_taxid": ptid,
            "profile_genus": pgen,
            "ref_accession": acc,
            "ref_taxid": rtid,
            "ref_organism": org,
            "ref_class_taxid": cls_id or "",
            "ref_class_name": cls_name or "",
            "query_accession": exclude,
        })

    for ptid in SAME_SPECIES:
        exclude = query_map[ptid]
        acc, rtid, org = pick_same_species_alt(
            ptid, exclude, profile_taxa, lineage_cache, used_acc,
        )
        used_acc.add(acc)
        rid += 1
        short = f"ref_{rid:02d}"
        cls_id, cls_name = taxonomy_class(int(rtid))
        gen = genus_name_for_taxid(ptid)
        rows.append({
            "short_id": short,
            "role": "profile_same_species",
            "profile_taxid": ptid,
            "profile_genus": gen,
            "ref_accession": acc,
            "ref_taxid": rtid,
            "ref_organism": org,
            "ref_class_taxid": cls_id or "",
            "ref_class_name": cls_name or "",
            "query_accession": exclude,
        })

    # Genera from same-species profile taxa (same_species tier)
    genera_pool: list[str] = []
    for ptid in SAME_SPECIES:
        g = genus_name_for_taxid(ptid)
        if g:
            genera_pool.append(g)

    random.seed(RNG_SEED)
    sampled_genera = [random.choice(genera_pool) for _ in range(EXTRA_COUNT)]

    for i, gterm in enumerate(sampled_genera):
        acc, rtid, org = pick_genus_sibling(
            gterm, profile_taxa, lineage_cache, used_acc, prefer_complete=True,
        )
        used_acc.add(acc)
        rid += 1
        short = f"ref_{rid:02d}"
        cls_id, cls_name = taxonomy_class(int(rtid))
        rows.append({
            "short_id": short,
            "role": "extra_same_genus_pool",
            "profile_taxid": "",
            "profile_genus": gterm,
            "ref_accession": acc,
            "ref_taxid": rtid,
            "ref_organism": org,
            "ref_class_taxid": cls_id or "",
            "ref_class_name": cls_name or "",
            "query_accession": "",
        })

    # Outgroups: class not in query_classes
    outgroup_bank = [
        ("Bacillus subtilis subsp. subtilis 168", "GCF_000009045.1"),
        ("Staphylococcus aureus subsp. aureus Mu50", "GCF_000013425.1"),
        ("Mycobacterium tuberculosis H37Rv", "GCF_000195955.2"),
        ("Bacteroides fragilis NCTC 9343", "GCF_000012825.1"),
        ("Fusobacterium nucleatum subsp. nucleatum ATCC 25586", "GCF_000008195.1"),
        ("Streptococcus pneumoniae TIGR4", "GCF_000006885.1"),
        ("Clostridioides difficile 630", "GCF_000007625.1"),
        ("Moorella thermoacetica ATCC 39073", "GCF_000018785.1"),
        ("Bifidobacterium longum NCC2705", "GCF_000022265.1"),
        ("Alistipes putredinis DSM 17216", "GCF_900129005.1"),
        ("Cutibacterium acnes KPA171202", "GCF_000008345.1"),
    ]

    added_og = 0
    for oname, acc in outgroup_bank:
        if added_og >= OUTGROUP_COUNT:
            break
        if acc in used_acc:
            continue
        tid = taxid_from_assembly_accession(acc)
        cls_id, cls_name = taxonomy_class(int(tid))
        if cls_id is None or cls_id in query_classes:
            continue
        used_acc.add(acc)
        rid += 1
        short = f"ref_{rid:02d}"
        rows.append({
            "short_id": short,
            "role": "outgroup",
            "profile_taxid": "",
            "profile_genus": "",
            "ref_accession": acc,
            "ref_taxid": tid,
            "ref_organism": oname,
            "ref_class_taxid": cls_id or "",
            "ref_class_name": cls_name or "",
            "query_accession": "",
        })
        added_og += 1

    if added_og < OUTGROUP_COUNT:
        sys.stderr.write(
            f"error: only {added_og}/{OUTGROUP_COUNT} class-disjoint outgroups; "
            "expand outgroup_bank\n"
        )
        sys.exit(3)

    assert len(rows) == 40, len(rows)

    tax_cache: dict[str, str] = {}
    fill_missing_ref_organism(rows, tax_cache)
    assign_genome_ids(rows)
    write_panel_tables(rows, design_path, manifest_path, input_map_path)

    sys.stderr.write(f"Design/manifest: {design_path}\n")
    sys.stderr.write(f"Downloading {len(rows)} genomes to {outdir} ...\n")

    for r in rows:
        dest = outdir / f"{r['genome_id']}.fna"
        sys.stderr.write(f"  {r['genome_id']} {r['ref_accession']}\n")
        try:
            download_fasta(r["ref_accession"], dest)
        except urllib.error.HTTPError as e:
            sys.stderr.write(f"    HTTPError: {e}\n")
            raise
        time.sleep(0.35)

    sys.stderr.write("Done.\n")


if __name__ == "__main__":
    main()
