#!/usr/bin/env python3
"""
fetch_genomes.py - download a curated set of marine reference genomes from NCBI.

Inputs (one mode):
    --species-list  TSV: <short_id>\\t<species_or_strain_name>\\t<genome_kind>
                    genome_kind is one of: reference | holdout
    --accession-tsv TSV with header: taxid, species, assembly_accession, ...
                    (e.g. profile_species_accessions.tsv). Writes <taxid>.fna.

    --outdir        directory to populate with <short_id>.fna files

Implementation:
    Prefer NCBI Datasets CLI (`datasets`) when installed; otherwise download
    via the public NCBI Datasets REST API. For --species-list, resolve the
    reference or complete assembly per taxon name.
"""

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.request
import zipfile
from pathlib import Path

NCBI_GENOME_DOWNLOAD_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/download"


def run(cmd, **kw):
    return subprocess.run(cmd, check=False, capture_output=True, text=True, **kw)


def best_accession(taxon):
    for flags in (["--reference"], ["--assembly-source", "RefSeq", "--assembly-level", "complete"],
                  ["--assembly-level", "complete"], []):
        cp = run(["datasets", "summary", "genome", "taxon", taxon,
                  "--as-json-lines", *flags])
        if cp.returncode != 0 or not cp.stdout.strip():
            continue
        for line in cp.stdout.splitlines():
            try:
                rec = json.loads(line)
            except Exception:
                continue
            acc = rec.get("accession")
            if acc and acc.startswith(("GCF_", "GCA_")):
                return acc, rec.get("organism", {}).get("organism_name", "")
    return None, None


def _extract_first_fna_from_zip(zip_path: Path, out_fa: Path) -> bool:
    with zipfile.ZipFile(zip_path) as zf:
        for name in zf.namelist():
            if name.endswith(".fna"):
                with zf.open(name) as src, open(out_fa, "wb") as dst:
                    shutil.copyfileobj(src, dst)
                return True
    return False


def download_one_datasets_cli(accession, out_fa: Path) -> bool:
    with tempfile.TemporaryDirectory() as td:
        zip_path = Path(td) / "g.zip"
        cp = run(["datasets", "download", "genome", "accession", accession,
                  "--include", "genome", "--filename", str(zip_path)])
        if cp.returncode != 0 or not zip_path.exists():
            sys.stderr.write(f"  download failed for {accession}: {cp.stderr.strip()[:200]}\n")
            return False
        if not _extract_first_fna_from_zip(zip_path, out_fa):
            sys.stderr.write(f"  no .fna inside zip for {accession}\n")
            return False
    return True


def download_one_http(accession, out_fa: Path) -> bool:
    body = json.dumps({
        "accessions": [accession],
        "include_annotation_type": ["GENOME_FASTA"],
    }).encode("utf-8")
    req = urllib.request.Request(
        NCBI_GENOME_DOWNLOAD_URL,
        data=body,
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    with tempfile.TemporaryDirectory() as td:
        zip_path = Path(td) / "g.zip"
        try:
            with urllib.request.urlopen(req, timeout=300) as resp:
                if resp.status != 200:
                    sys.stderr.write(f"  HTTP {resp.status} for {accession}\n")
                    return False
                with open(zip_path, "wb") as out:
                    shutil.copyfileobj(resp, out)
        except urllib.error.HTTPError as e:
            sys.stderr.write(f"  HTTP error for {accession}: {e.code} {e.reason}\n")
            return False
        except urllib.error.URLError as e:
            sys.stderr.write(f"  network error for {accession}: {e.reason}\n")
            return False
        if not _extract_first_fna_from_zip(zip_path, out_fa):
            sys.stderr.write(f"  no .fna inside zip for {accession}\n")
            return False
    return True


def download_one(accession, short_id, outdir):
    out_fa = outdir / f"{short_id}.fna"
    if out_fa.exists() and out_fa.stat().st_size > 0:
        return True
    if shutil.which("datasets"):
        return download_one_datasets_cli(accession, out_fa)
    return download_one_http(accession, out_fa)


def main():
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--species-list")
    g.add_argument("--accession-tsv",
                   help="TSV with header containing taxid, species, assembly_accession")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--manifest", default=None,
                    help="Output TSV recording id, accession, taxon, genome_kind.")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    manifest_path = Path(args.manifest) if args.manifest else outdir / "manifest.tsv"

    rows = []
    if args.accession_tsv:
        import csv
        with open(args.accession_tsv, newline="") as fh:
            rdr = csv.DictReader(fh, delimiter="\t")
            lower = {k.lower(): k for k in rdr.fieldnames or []}
            k_tid = lower.get("taxid")
            k_acc = lower.get("assembly_accession") or lower.get("accession")
            k_spec = lower.get("species") or lower.get("organism")
            if not k_tid or not k_acc:
                sys.stderr.write("accession-tsv needs taxid and assembly_accession columns\n")
                sys.exit(2)
            for rec in rdr:
                tid = (rec.get(k_tid) or "").strip()
                acc = (rec.get(k_acc) or "").strip()
                spec = (rec.get(k_spec) or "").strip() if k_spec else ""
                if not tid or not acc:
                    continue
                if not tid.isdigit():
                    continue
                short_id = tid
                taxon = spec or tid
                kind = "reference"
                rows.append((short_id, taxon, kind, acc))
    else:
        with open(args.species_list) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line or line.startswith("#") or line.startswith("@"):
                    continue
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                short_id, taxon, kind = parts[0], parts[1], parts[2]
                rows.append((short_id, taxon, kind, None))

    sys.stderr.write(f"Resolving and downloading {len(rows)} genomes...\n")
    success = []
    failed = []
    for row in rows:
        short_id, taxon, kind, fixed_acc = row
        if fixed_acc:
            sys.stderr.write(f"  {short_id} <- {fixed_acc} ({taxon})\n")
            acc, org_name = fixed_acc, taxon
        else:
            sys.stderr.write(f"  {short_id} <- {taxon} ({kind})\n")
            acc, org_name = best_accession(taxon)
        if not acc:
            sys.stderr.write(f"    no accession found for {taxon}\n")
            failed.append((short_id, taxon, kind, "", ""))
            continue
        ok = download_one(acc, short_id, outdir)
        if not shutil.which("datasets"):
            time.sleep(0.35)
        if ok:
            success.append((short_id, taxon, kind, acc, org_name))
            sys.stderr.write(f"    {acc} ({org_name})\n")
        else:
            failed.append((short_id, taxon, kind, acc, org_name))

    with open(manifest_path, "w") as fh:
        fh.write("short_id\ttaxon\tgenome_kind\taccession\torganism\n")
        for r in success + failed:
            fh.write("\t".join(str(x) for x in r) + "\n")

    sys.stderr.write(f"\nSuccess: {len(success)}/{len(rows)}; failed: {len(failed)}\n")
    if failed:
        sys.stderr.write("Failed entries:\n")
        for r in failed:
            sys.stderr.write("  " + "\t".join(str(x) for x in r) + "\n")
        sys.exit(2 if not success else 0)


if __name__ == "__main__":
    main()
