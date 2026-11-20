#!/usr/bin/env python3
"""
Build gene -> functional term tables (GO BP/MF/CC and KEGG pathways) for the pipeline's gene IDs.

Inputs (downloaded by 08_functional_annotation.sh into <rawdir>):
  go-basic.obo                      GO ontology (is_a + part_of, used to propagate annotations to ancestors)
  Calb_GO.gaf.gz                    C. albicans GO annotations (CGD GAF; GO Consortium CANAL-mod or CGD site)
  <sp>_uniprot.tsv.gz               UniProt: accession, ordered locus, ORF names, gene names, GO P/F/C
  <sp>_kegg_link_pathway.txt        KEGG gene -> pathway
  <sp>_kegg_list_pathway.txt        KEGG pathway names
  <sp>_kegg_conv_uniprot.txt        KEGG gene -> UniProt accession (identifier bridge)
Gene identifiers are taken from gene_annotation.tsv and the GFF3 (locus_tag, old_locus_tag, Name, gene,
gene_synonym, Dbxref). For C. albicans the Assembly-22 systematic name (C1_00010W_A) is derived from the
NCBI locus tag (CAALFM_C100010WA) and CaO19.N / orf19.N aliases are unified.

Output: <outdir>/<sp>_term2gene.tsv.gz  (database, term_id, term_name, gene_id)
        <outdir>/annotation_coverage.tsv, <outdir>/versions.txt
"""
import argparse, collections, gzip, os, re, sys
from urllib.parse import unquote

SPECIES_GO_SOURCE = {"Calb": "gaf+uniprot", "Saur": "uniprot"}
ROOTS = {"GO:0008150", "GO:0003674", "GO:0005575"}
NS = {"biological_process": "GO_BP", "molecular_function": "GO_MF", "cellular_component": "GO_CC"}
ASPECT = {"P": "GO_BP", "F": "GO_MF", "C": "GO_CC"}


def op(p):
    return gzip.open(p, "rt", errors="replace") if p.endswith(".gz") else open(p, errors="replace")


def norm(tok):
    t = tok.strip().strip('"').strip()
    if not t:
        return None
    t = re.sub(r"^CaO19\.", "orf19.", t, flags=re.I)
    return t.upper()


# ---------------------------------------------------------------- ontology
def parse_obo(path):
    terms, alt, cur, version = {}, {}, None, "NA"
    for line in op(path):
        line = line.rstrip("\n")
        if line.startswith("data-version:"):
            version = line.split(":", 1)[1].strip()
        if line == "[Term]":
            cur = {"parents": set(), "obsolete": False, "replaced_by": None}
            continue
        if line.startswith("[") and line != "[Term]":
            cur = None; continue
        if cur is None or ":" not in line:
            continue
        k, v = line.split(":", 1); v = v.strip()
        if k == "id":
            cur["id"] = v; terms[v] = cur
        elif k == "name": cur["name"] = v
        elif k == "namespace": cur["ns"] = v
        elif k == "alt_id": alt[v] = cur.get("id")
        elif k == "is_a": cur["parents"].add(v.split()[0])
        elif k == "relationship" and v.startswith("part_of "): cur["parents"].add(v.split()[1])
        elif k == "is_obsolete" and v == "true": cur["obsolete"] = True
        elif k == "replaced_by": cur["replaced_by"] = v.split()[0]
    return terms, alt, version


def make_resolver(terms, alt):
    anc_cache = {}

    def primary(t):
        t = alt.get(t, t)
        x = terms.get(t)
        if x is None:
            return None
        if x["obsolete"]:
            return primary(x["replaced_by"]) if x.get("replaced_by") else None
        return t

    def ancestors(t):
        if t in anc_cache:
            return anc_cache[t]
        out = {t}
        for p in terms[t]["parents"]:
            if p in terms:
                out |= ancestors(p)
        anc_cache[t] = out
        return out
    return primary, ancestors


# ---------------------------------------------------------------- gene identifiers
def gene_tokens(annot_tsv, gff, species):
    """return token -> set(gene_id) for one species"""
    genes = {}
    with open(annot_tsv) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(hdr, line.rstrip("\n").split("\t")))
            if r["species"] != species:
                continue
            toks = {r["gene_id"], r.get("name", ""), r.get("old_locus_tag", "")}
            toks |= set(re.split(r"[,;| ]+", r.get("alias", "")))
            genes[r["gene_id"]] = toks
    # extra identifiers from the GFF (Dbxref, gene_synonym) keyed by locus_tag
    if gff and os.path.exists(gff):
        for line in op(gff):
            if line.startswith("#") or "\t" not in line:
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("gene", "pseudogene", "mRNA", "CDS"):
                continue
            a = dict(kv.split("=", 1) for kv in f[8].strip(";").split(";") if "=" in kv)
            lt = a.get("locus_tag")
            if lt not in genes:
                continue
            for key in ("Dbxref", "gene_synonym", "Alias", "Name", "gene", "old_locus_tag", "protein_id"):
                for v in unquote(a.get(key, "")).split(","):
                    genes[lt].add(v); genes[lt].add(v.split(":", 1)[-1])
    tok2gene = collections.defaultdict(set)
    for g, toks in genes.items():
        if species == "Calb":
            m = re.match(r"^CAALFM_C([1-7R])(\d{5})([WC])([AB])$", g)
            if m:
                toks.add(f"C{m.group(1)}_{m.group(2)}{m.group(3)}_{m.group(4)}")
        for t in toks:
            n = norm(t)
            if n and n not in ("NA", "-"):
                tok2gene[n].add(g)
    # keep only unambiguous tokens
    return {t: next(iter(gs)) for t, gs in tok2gene.items() if len(gs) == 1}, set(genes)


def map_tokens(cands, tok2gene):
    for c in cands:
        n = norm(c)
        if n and n in tok2gene:
            return tok2gene[n]
    return None


# ---------------------------------------------------------------- sources
def read_kegg(raw, sp, tok2gene):
    link, lst, conv = (f"{raw}/{sp}_kegg_link_pathway.txt", f"{raw}/{sp}_kegg_list_pathway.txt",
                       f"{raw}/{sp}_kegg_conv_uniprot.txt")
    up2gene, names, out = {}, {}, collections.defaultdict(set)
    if os.path.exists(conv):
        for line in op(conv):
            p = line.rstrip("\n").split("\t")
            if len(p) == 2:
                g = map_tokens([p[0].split(":", 1)[-1]], tok2gene)
                if g: up2gene[p[1].split(":", 1)[-1]] = g
    if not os.path.exists(link):
        return out, names, up2gene
    rows = [l.rstrip("\n").split("\t") for l in op(lst)] if os.path.exists(lst) else []
    suf = collections.Counter(r[1].rsplit(" - ", 1)[1] for r in rows if len(r) == 2 and " - " in r[1])
    common = suf.most_common(1)[0][0] if suf else None
    for r in rows:
        if len(r) != 2: continue
        pid = r[0].replace("path:", "")
        nm = r[1][: -len(" - " + common)] if common and r[1].endswith(" - " + common) else r[1]
        names[pid] = nm
    for line in op(link):
        p = line.rstrip("\n").split("\t")
        if len(p) != 2: continue
        g = map_tokens([p[0].split(":", 1)[-1]], tok2gene)
        pid = p[1].replace("path:", "")
        num = re.search(r"(\d{5})$", pid)
        if g and not (num and 1100 <= int(num.group(1)) <= 1299):   # drop global/overview maps (01100-01299)
            out[pid].add(g)
    return out, names, up2gene


def read_uniprot(raw, sp, tok2gene, up2gene, primary):
    path = f"{raw}/{sp}_uniprot.tsv.gz"
    ann = collections.defaultdict(set)
    if not os.path.exists(path):
        return ann, 0
    fh = op(path); hdr = fh.readline().rstrip("\n").split("\t")
    idx = {h: i for i, h in enumerate(hdr)}
    n_mapped = 0
    for line in fh:
        f = line.rstrip("\n").split("\t")
        acc = f[idx.get("Entry", 0)]
        g = up2gene.get(acc) or map_tokens([acc], tok2gene)
        if not g:
            cands = []
            # header labels differ between UniProt versions ("Ordered locus name" / "Gene Names (ordered locus)" ...)
            idcols = sorted((c for c in hdr if re.search(r"ordered locus|orf|^gene names", c, re.I)),
                            key=lambda c: (0 if re.search("locus", c, re.I) else 1 if re.search("orf", c, re.I) else 2))
            for col in idcols:
                if idx[col] < len(f):
                    cands += re.split(r"[ ;,/]+", f[idx[col]])
            g = map_tokens(cands, tok2gene)
        if not g:
            continue
        n_mapped += 1
        for col in hdr:
            if col.startswith("Gene Ontology (") and idx[col] < len(f):
                for go in re.findall(r"GO:\d{7}", f[idx[col]]):
                    p = primary(go)
                    if p: ann[g].add(p)
    return ann, n_mapped


def read_gaf(path, tok2gene, primary, taxon="237561"):
    ann = collections.defaultdict(set); n = 0
    if not os.path.exists(path):
        return ann, 0
    for line in op(path):
        if line.startswith("!"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 13 or "NOT" in f[3].upper().split("|"):
            continue
        if taxon and taxon not in f[12]:
            continue
        cands = [f[1], f[2]] + f[10].split("|")
        g = map_tokens(cands, tok2gene)
        if not g: continue
        p = primary(f[4])
        if p: ann[g].add(p); n += 1
    return ann, n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True); ap.add_argument("--annot", required=True)
    ap.add_argument("--calb-gff"); ap.add_argument("--saur-gff"); ap.add_argument("--outdir", required=True)
    a = ap.parse_args(); os.makedirs(a.outdir, exist_ok=True)
    terms, alt, version = parse_obo(f"{a.raw}/go-basic.obo")
    primary, ancestors = make_resolver(terms, alt)
    cov = []
    for sp, gff in (("Calb", a.calb_gff), ("Saur", a.saur_gff)):
        tok2gene, allgenes = gene_tokens(a.annot, gff, sp)
        kegg, knames, up2gene = read_kegg(a.raw, sp, tok2gene)
        go_direct = collections.defaultdict(set); src = []
        if sp == "Calb":
            gaf, n = read_gaf(f"{a.raw}/Calb_GO.gaf.gz", tok2gene, primary)
            for g, s in gaf.items(): go_direct[g] |= s
            src.append(f"CGD GAF ({n} annotations mapped)")
        up, n = read_uniprot(a.raw, sp, tok2gene, up2gene, primary)
        for g, s in up.items(): go_direct[g] |= s
        src.append(f"UniProt ({n} entries mapped)")
        rows = []
        per_db = collections.defaultdict(set)
        for g, gos in go_direct.items():
            full = set()
            for t in gos: full |= ancestors(t)
            for t in full - ROOTS:
                db = NS.get(terms[t].get("ns"))
                if db:
                    rows.append((db, t, terms[t].get("name", t), g)); per_db[db].add(g)
        for pid, gs in kegg.items():
            for g in gs:
                rows.append(("KEGG", pid, knames.get(pid, pid), g)); per_db["KEGG"].add(g)
        rows = sorted(set(rows))
        with gzip.open(f"{a.outdir}/{sp}_term2gene.tsv.gz", "wt") as fo:
            fo.write("database\tterm_id\tterm_name\tgene_id\n")
            for r in rows: fo.write("\t".join(r) + "\n")
        for db in ("GO_BP", "GO_MF", "GO_CC", "KEGG"):
            nterm = len({r[1] for r in rows if r[0] == db})
            cov.append((sp, db, len(allgenes), len(per_db[db]), round(100 * len(per_db[db]) / max(len(allgenes), 1), 1), nterm))
        print(f"[{sp}] GO sources: {'; '.join(src)}; KEGG pathways: {len(kegg)}", file=sys.stderr)
    with open(f"{a.outdir}/annotation_coverage.tsv", "w") as fo:
        fo.write("species\tdatabase\tgenes_total\tgenes_annotated\tpct_annotated\tterms\n")
        for c in cov: fo.write("\t".join(map(str, c)) + "\n")
    with open(f"{a.outdir}/versions.txt", "a") as fo:
        fo.write(f"GO ontology data-version: {version}\n")
    for c in cov: print("\t".join(map(str, c)), file=sys.stderr)


if __name__ == "__main__":
    main()
