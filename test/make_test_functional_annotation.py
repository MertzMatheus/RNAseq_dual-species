#!/usr/bin/env python3
"""
Write fake functional-annotation downloads (same formats as GO/CGD/UniProt/KEGG) for the simulated project,
so steps 08-09 can be tested offline. Plants a KEGG pathway made of the S. aureus genes simulated as
UP with C. albicans (truth_Saur_dual_effect.txt) and a decoy pathway of random genes.
Usage: make_test_functional_annotation.py <project_dir>
"""
import gzip, os, random, sys
root = sys.argv[1]; A = f"{root}/analysis"; raw = f"{A}/reference/functional/raw"
os.makedirs(raw, exist_ok=True); random.seed(3)
ann = [l.rstrip("\n").split("\t") for l in open(f"{A}/reference/gene_annotation.tsv")]
H = ann[0]; ann = [dict(zip(H, r)) for r in ann[1:]]
calb = [r for r in ann if r["species"] == "Calb"]; saur = [r for r in ann if r["species"] == "Saur"]
truth = [float(x) for x in open(f"{root}/truth_Saur_dual_effect.txt")]

# --- ontology: small tree under the three roots, with an alt_id and an obsolete term
obo = ["format-version: 1.2", "data-version: releases/2026-10-01", ""]
def term(i, name, ns, parents=(), extra=()):
    obo.extend(["[Term]", f"id: {i}", f"name: {name}", f"namespace: {ns}", *[f"is_a: {p} ! x" for p in parents], *extra, ""])
term("GO:0008150", "biological_process", "biological_process")
term("GO:0003674", "molecular_function", "molecular_function")
term("GO:0005575", "cellular_component", "cellular_component")
term("GO:0050896", "response to stimulus", "biological_process", ["GO:0008150"])
term("GO:0009607", "response to biotic stimulus", "biological_process", ["GO:0050896"], ["alt_id: GO:0009999"])
term("GO:0008152", "metabolic process", "biological_process", ["GO:0008150"])
bp = []
for k in range(25):
    i = f"GO:10{k:05d}"; term(i, f"simulated process {k}", "biological_process", [random.choice(["GO:0009607", "GO:0008152"])]); bp.append(i)
mf = []
for k in range(10):
    i = f"GO:20{k:05d}"; term(i, f"simulated activity {k}", "molecular_function", ["GO:0003674"]); mf.append(i)
term("GO:0005737", "cytoplasm", "cellular_component", ["GO:0005575"])
term("GO:0000000", "obsolete thing", "biological_process", extra=["is_obsolete: true", "replaced_by: GO:0009607"])
open(f"{raw}/go-basic.obo", "w").write("\n".join(obo) + "\n")

# --- C. albicans GAF 2.2 (CGD style): mapped via symbol (gene name) or via a synonym holding the locus tag
with gzip.open(f"{raw}/Calb_GO.gaf.gz", "wt") as fo:
    fo.write("!gaf-version: 2.2\n")
    for n, r in enumerate(calb):
        sym = r["name"] or f"orf19.{1000+n}"
        syn = f"orf19.{1000+n}|{r['gene_id']}" if not r["name"] else f"orf19.{1000+n}"
        for go in random.sample(bp, 3) + random.sample(mf, 1) + ["GO:0005737"] + (["GO:0009999"] if n % 7 == 0 else []):
            q = "NOT|involved_in" if n % 50 == 1 else "involved_in"
            fo.write("\t".join(["CGD", f"CAL{n:07d}", sym, q, go, "CGD_REF:1", "IEA", "", "P", "", syn, "gene",
                                "taxon:237561", "20261001", "CGD", "", ""]) + "\n")

# --- UniProt TSV (new header labels); S. aureus mapped via ordered locus name = old locus tag
hdr = "Entry\tGene Names (ordered locus)\tGene Names (ORF)\tGene Names\tGene Ontology (biological process)\tGene Ontology (molecular function)\tGene Ontology (cellular component)\n"
with gzip.open(f"{raw}/Saur_uniprot.tsv.gz", "wt") as fo:
    fo.write(hdr)
    for n, r in enumerate(saur):
        gos = random.sample(bp, 2)
        fo.write(f"Q{n:05d}\t{r['old_locus_tag']}\t\t{r['name']}\t" + "; ".join(f"x [{g}]" for g in gos) +
                 f"\tbinding [{random.choice(mf)}]\tcytoplasm [GO:0005737]\n")
with gzip.open(f"{raw}/Calb_uniprot.tsv.gz", "wt") as fo:
    fo.write(hdr)
    for n, r in enumerate(calb[:20]):     # mapped only through the KEGG conv bridge (accession)
        fo.write(f"A0A{n:05d}\t\t\t\tx [{random.choice(bp)}]\t\t\n")

# --- KEGG REST files
up = [r["old_locus_tag"] for r, t in zip(saur, truth) if t > 0]
rnd = random.sample([r["old_locus_tag"] for r in saur], 25)
with open(f"{raw}/Saur_kegg_link_pathway.txt", "w") as fo:
    for g in up[: max(6, int(len(up) * .8))] + random.sample([r["old_locus_tag"] for r in saur], 6):
        fo.write(f"saa:{g}\tpath:saa99901\n")
    for g in rnd: fo.write(f"saa:{g}\tpath:saa99902\n")
    for g in rnd[:10]: fo.write(f"saa:{g}\tpath:saa01100\n")          # global map -> must be dropped
sfx = " - Staphylococcus aureus subsp. aureus USA300_FPR3757 (CA-MRSA)"
open(f"{raw}/Saur_kegg_list_pathway.txt", "w").write(
    f"saa99901\tSimulated C. albicans-induced pathway{sfx}\nsaa99902\tDecoy random pathway{sfx}\nsaa01100\tMetabolic pathways{sfx}\n")
open(f"{raw}/Saur_kegg_conv_uniprot.txt", "w").write("".join(f"saa:{r['old_locus_tag']}\tup:Q{n:05d}\n" for n, r in enumerate(saur)))
with open(f"{raw}/Calb_kegg_link_pathway.txt", "w") as fo:
    for r in random.sample(calb, 30): fo.write(f"cal:{r['gene_id']}\tpath:cal00010\n")
open(f"{raw}/Calb_kegg_list_pathway.txt", "w").write("cal00010\tGlycolysis / Gluconeogenesis - Candida albicans\n")
open(f"{raw}/Calb_kegg_conv_uniprot.txt", "w").write("".join(f"cal:{r['gene_id']}\tup:A0A{n:05d}\n" for n, r in enumerate(calb[:20])))
print(f"planted pathway: {len(up)} true S. aureus up genes; files in {raw}")
