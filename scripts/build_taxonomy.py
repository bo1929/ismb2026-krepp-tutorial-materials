#!/usr/bin/env python3
"""
build_taxonomy.py - assemble metadata.tsv and lineages.tsv from a manifest.

For tutorial use we need stable, fully-named GTDB-style lineages so that the
taxonomic placement step does not encounter unclassified bacteria. We map each
short_id to a hand-curated lineage based on the well-known taxonomy of marine
type strains (no reliance on NCBI taxdump).
"""

import argparse
import sys

LINEAGE = {
    # Alphaproteobacteria - SAR11 / Pelagibacterales
    "SAR11_HTCC1062": ("Bacteria", "Pseudomonadota", "Alphaproteobacteria",
                       "Pelagibacterales", "Pelagibacteraceae", "Pelagibacter",
                       "Pelagibacter ubique HTCC1062"),
    "SAR11_HTCC7211": ("Bacteria", "Pseudomonadota", "Alphaproteobacteria",
                       "Pelagibacterales", "Pelagibacteraceae", "Pelagibacter",
                       "Pelagibacter sp. HTCC7211"),
    "SAR11_IMCC9063": ("Bacteria", "Pseudomonadota", "Alphaproteobacteria",
                       "Pelagibacterales", "Pelagibacteraceae", "Pelagibacter",
                       "Pelagibacter sp. IMCC9063"),
    # Alphaproteobacteria - Roseobacter clade (Rhodobacterales)
    "ROSE_DSS3": ("Bacteria", "Pseudomonadota", "Alphaproteobacteria",
                  "Rhodobacterales", "Rhodobacteraceae", "Ruegeria",
                  "Ruegeria pomeroyi DSS-3"),
    "ROSE_OCh114": ("Bacteria", "Pseudomonadota", "Alphaproteobacteria",
                    "Rhodobacterales", "Rhodobacteraceae", "Roseobacter",
                    "Roseobacter denitrificans OCh114"),
    "ROSE_DG874": ("Bacteria", "Pseudomonadota", "Alphaproteobacteria",
                   "Rhodobacterales", "Rhodobacteraceae", "Roseovarius",
                   "Roseovarius nubinhibens ISM"),
    "ROSE_PHAEO": ("Bacteria", "Pseudomonadota", "Alphaproteobacteria",
                   "Rhodobacterales", "Rhodobacteraceae", "Phaeobacter",
                   "Phaeobacter inhibens DSM 17395"),
    # Cyanobacteria - Prochlorococcus
    "PROC_MED4": ("Bacteria", "Cyanobacteriota", "Cyanophyceae",
                  "Synechococcales", "Prochlorococcaceae", "Prochlorococcus",
                  "Prochlorococcus marinus MED4"),
    "PROC_MIT9301": ("Bacteria", "Cyanobacteriota", "Cyanophyceae",
                     "Synechococcales", "Prochlorococcaceae", "Prochlorococcus",
                     "Prochlorococcus marinus MIT9301"),
    "PROC_MIT9313": ("Bacteria", "Cyanobacteriota", "Cyanophyceae",
                     "Synechococcales", "Prochlorococcaceae", "Prochlorococcus",
                     "Prochlorococcus marinus MIT9313"),
    "PROC_NATL2A": ("Bacteria", "Cyanobacteriota", "Cyanophyceae",
                    "Synechococcales", "Prochlorococcaceae", "Prochlorococcus",
                    "Prochlorococcus marinus NATL2A"),
    "PROC_AS9601": ("Bacteria", "Cyanobacteriota", "Cyanophyceae",
                    "Synechococcales", "Prochlorococcaceae", "Prochlorococcus",
                    "Prochlorococcus marinus AS9601"),
    # Cyanobacteria - Synechococcus
    "SYNE_WH8102": ("Bacteria", "Cyanobacteriota", "Cyanophyceae",
                    "Synechococcales", "Synechococcaceae", "Synechococcus",
                    "Synechococcus sp. WH8102"),
    "SYNE_CC9311": ("Bacteria", "Cyanobacteriota", "Cyanophyceae",
                    "Synechococcales", "Synechococcaceae", "Synechococcus",
                    "Synechococcus sp. CC9311"),
    "SYNE_WH7803": ("Bacteria", "Cyanobacteriota", "Cyanophyceae",
                    "Synechococcales", "Synechococcaceae", "Synechococcus",
                    "Synechococcus sp. WH7803"),
    # Bacteroidota - Flavobacteriaceae
    "FLAV_PSYCHRO": ("Bacteria", "Bacteroidota", "Flavobacteriia",
                     "Flavobacteriales", "Flavobacteriaceae", "Flavobacterium",
                     "Flavobacterium psychrophilum JIP02/86"),
    "FLAV_JOHNS": ("Bacteria", "Bacteroidota", "Flavobacteriia",
                   "Flavobacteriales", "Flavobacteriaceae", "Flavobacterium",
                   "Flavobacterium johnsoniae UW101"),
    "FLAV_POLAR23P": ("Bacteria", "Bacteroidota", "Flavobacteriia",
                      "Flavobacteriales", "Flavobacteriaceae", "Polaribacter",
                      "Polaribacter irgensii 23-P"),
    "FLAV_POLAR_MED152": ("Bacteria", "Bacteroidota", "Flavobacteriia",
                          "Flavobacteriales", "Flavobacteriaceae", "Polaribacter",
                          "Polaribacter sp. MED152"),
    "FLAV_DOKDONIA_MED134": ("Bacteria", "Bacteroidota", "Flavobacteriia",
                             "Flavobacteriales", "Flavobacteriaceae", "Dokdonia",
                             "Dokdonia sp. MED134"),
    "FLAV_GRAMELLA": ("Bacteria", "Bacteroidota", "Flavobacteriia",
                      "Flavobacteriales", "Flavobacteriaceae", "Christiangramia",
                      "Christiangramia forsetii KT0803"),
    "FLAV_KORDIA": ("Bacteria", "Bacteroidota", "Flavobacteriia",
                    "Flavobacteriales", "Flavobacteriaceae", "Kordia",
                    "Kordia algicida OT-1"),
    # Gammaproteobacteria
    "GAMMA_VIBRIO": ("Bacteria", "Pseudomonadota", "Gammaproteobacteria",
                     "Vibrionales", "Vibrionaceae", "Vibrio",
                     "Vibrio cholerae N16961"),
    "GAMMA_ALTERO": ("Bacteria", "Pseudomonadota", "Gammaproteobacteria",
                     "Alteromonadales", "Alteromonadaceae", "Alteromonas",
                     "Alteromonas macleodii ATCC 27126"),
    "ALTERO_ATL": ("Bacteria", "Pseudomonadota", "Gammaproteobacteria",
                   "Alteromonadales", "Alteromonadaceae", "Alteromonas",
                   "Alteromonas mediterranea"),
    "GAMMA_PSEUDO": ("Bacteria", "Pseudomonadota", "Gammaproteobacteria",
                     "Alteromonadales", "Pseudoalteromonadaceae", "Pseudoalteromonas",
                     "Pseudoalteromonas haloplanktis"),
    "GAMMA_MARINOBACTER": ("Bacteria", "Pseudomonadota", "Gammaproteobacteria",
                           "Pseudomonadales", "Marinobacteraceae", "Marinobacter",
                           "Marinobacter nauticus ATCC 49840"),
    # Thaumarchaeota
    "THAUM_NITROSO": ("Archaea", "Thermoproteota", "Nitrososphaeria",
                      "Nitrosopumilales", "Nitrosopumilaceae", "Nitrosopumilus",
                      "Nitrosopumilus maritimus SCM1"),
    "THAUM_NITROSOARCH": ("Archaea", "Thermoproteota", "Nitrososphaeria",
                          "Nitrosopumilales", "Nitrosopumilaceae", "Nitrosarchaeum",
                          "Nitrosarchaeum koreense MY1"),
    "THAUM_NITROSARCH": ("Archaea", "Thermoproteota", "Nitrososphaeria",
                         "Nitrososphaerales", "Nitrososphaeraceae", "Nitrososphaera",
                         "Nitrososphaera viennensis EN76"),
    # Outgroups (non-marine, distantly related)
    "OUT_ECOLI": ("Bacteria", "Pseudomonadota", "Gammaproteobacteria",
                  "Enterobacterales", "Enterobacteriaceae", "Escherichia",
                  "Escherichia coli K-12 MG1655"),
    "OUT_BSUBTILIS": ("Bacteria", "Bacillota", "Bacilli",
                      "Bacillales", "Bacillaceae", "Bacillus",
                      "Bacillus subtilis subsp. subtilis 168"),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--lineages", required=True)
    args = ap.parse_args()

    rows = []
    with open(args.manifest) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5 or not parts[3]:
                continue
            row = dict(zip(header, parts))
            if row["short_id"] not in LINEAGE:
                sys.stderr.write("WARN: no lineage curated for {}\n".format(row["short_id"]))
                continue
            rows.append(row)

    rows.sort(key=lambda r: r["short_id"])

    with open(args.metadata, "w") as fh:
        fh.write("genome\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\trole\taccession\n")
        for r in rows:
            lin = LINEAGE[r["short_id"]]
            fh.write("\t".join([r["short_id"], *lin, r["role"], r["accession"]]) + "\n")

    with open(args.lineages, "w") as fh:
        for r in rows:
            k, p, c, o, fam, g, s = LINEAGE[r["short_id"]]
            label = "k__{}; p__{}; c__{}; o__{}; f__{}; g__{}; s__{}".format(k, p, c, o, fam, g, s)
            fh.write("{}\t{}\n".format(r["short_id"], label))

    sys.stderr.write("Wrote {} rows to {} and {}\n".format(len(rows), args.metadata, args.lineages))


if __name__ == "__main__":
    main()
