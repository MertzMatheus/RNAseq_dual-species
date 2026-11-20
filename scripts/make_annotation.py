#!/usr/bin/env python3
"""
Build a combined, species-tagged annotation for dual RNA-seq from two GFF3 files.

Outputs (in --outdir):
  combined.saf            featureCounts SAF (GeneID Chr Start End Strand); eukaryote = exons, prokaryote = gene span
  gene_annotation.tsv     one row per gene: gene_id, species, name, old_locus_tag, alias, biotype, product, chr, start, end, strand, length
  calb_exons.gtf          exon GTF of the eukaryote (used by STAR for splice-junction annotation)
  contig_species.tsv      contig -> species (from the FASTA files)

Works with NCBI RefSeq/GenBank GFF3 and CGD (Candida Genome Database) GFF3.
The gene identifier is the locus_tag when present (e.g. CAALFM_C100010WA, SAUSA300_RS00005),
otherwise the GFF ID with "gene-" prefix removed.

Usage:
  make_annotation.py --euk-gff calb.gff --euk-fna calb.fna --euk-name Calb \
                     --prok-gff saur.gff --prok-fna saur.fna --prok-name Saur --outdir ref/
"""
import argparse, collections, gzip, os, sys
from urllib.parse import unquote

RNA_TYPES = {"mRNA", "transcript", "primary_transcript", "ncRNA", "rRNA", "tRNA", "snRNA", "snoRNA",
             "lnc_RNA", "lncRNA", "antisense_RNA", "RNase_P_RNA", "RNase_MRP_RNA", "SRP_RNA", "tmRNA",
             "telomerase_RNA", "pseudogenic_transcript", "misc_RNA", "guide_RNA"}
GENE_TYPES = {"gene", "pseudogene", "ncRNA_gene", "rRNA_gene", "tRNA_gene", "snRNA_gene", "snoRNA_gene",
              "transposable_element_gene", "blocked_reading_frame", "ORF"}
SKIP_TYPES = {"region", "chromosome", "contig", "supercontig", "databank_entry", "source", "repeat_region",
              "centromere", "long_terminal_repeat", "sequence_feature", "riboswitch", "binding_site",
              "origin_of_replication", "direct_repeat", "mobile_genetic_element", "intron",
              "five_prime_UTR_intron", "polyA_signal_sequence", "polyA_site", "uORF", "telomere",
              "repeat_family", "retrotransposon", "nucleotide_match", "regulatory_region"}


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def fasta_ids(path):
    ids = []
    with opener(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def parse_attrs(s):
    d = {}
    for part in s.strip().strip(";").split(";"):
        if not part.strip():
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            d[k.strip()] = unquote(v.strip())
    return d


def parse_gff(path):
    feats = {}          # id -> feature dict
    anon = []           # features without ID (exons/CDS often have no ID in some GFFs)
    n = 0
    with opener(path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9:
                continue
            seqid, src, ftype, start, end, score, strand, phase, attrs = f
            if ftype in SKIP_TYPES:
                continue
            a = parse_attrs(attrs)
            n += 1
            rec = dict(seqid=seqid, type=ftype, start=int(start), end=int(end),
                       strand=strand if strand in "+-" else "+", attrs=a,
                       parents=[p for p in a.get("Parent", "").split(",") if p])
            fid = a.get("ID")
            if fid and fid not in feats:
                rec["id"] = fid
                feats[fid] = rec
            else:  # discontinuous features (CDS split across lines share an ID) or no ID
                rec["id"] = None
                anon.append(rec)
    return feats, anon


def top_ancestor(fid, feats, depth=0):
    rec = feats.get(fid)
    if rec is None or depth > 10:
        return fid
    if not rec["parents"]:
        return fid
    return top_ancestor(rec["parents"][0], feats, depth + 1)


def clean_id(x):
    for pre in ("gene-", "gene:", "rna-", "cds-"):
        if x.startswith(pre):
            return x[len(pre):]
    return x


def merge_intervals(iv):
    iv = sorted(iv)
    out = []
    for s, e in iv:
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def build(gff, species, mode, valid_contigs):
    feats, anon = parse_gff(gff)
    allrecs = list(feats.values()) + anon
    genes = collections.OrderedDict()   # top_id -> info

    def gene_info(top):
        if top not in genes:
            g = feats.get(top)
            if g is None:
                return None
            genes[top] = dict(top=g, exons=[], cds=[], children=[])
        return genes[top]

    # register top-level genes
    for rec in feats.values():
        if not rec["parents"] and (rec["type"] in GENE_TYPES or rec["type"] in RNA_TYPES or rec["type"] == "CDS"):
            gene_info(rec["id"])
    # attach children
    for rec in allrecs:
        if not rec["parents"]:
            continue
        top = top_ancestor(rec["parents"][0], feats)
        gi = gene_info(top)
        if gi is None:
            continue
        if rec["type"] == "exon":
            gi["exons"].append(rec)
        elif rec["type"] == "CDS":
            gi["cds"].append(rec)
        elif rec["type"] in RNA_TYPES:
            gi["children"].append(rec)

    rows, saf, gtf = [], [], []
    seen = set()
    for top, gi in genes.items():
        g = gi["top"]
        if g["seqid"] not in valid_contigs:
            continue
        a = g["attrs"]
        gid = a.get("locus_tag") or clean_id(a.get("ID", top))
        if gid in seen:
            gid = f"{gid}_{g['start']}"
        seen.add(gid)
        # biotype
        biotype = a.get("gene_biotype") or a.get("biotype")
        child_types = [c["type"] for c in gi["children"]]
        if not biotype:
            if g["type"] == "pseudogene":
                biotype = "pseudogene"
            elif child_types:
                ct = child_types[0]
                biotype = "protein_coding" if ct == "mRNA" else ct
            elif gi["cds"] or g["type"] in ("CDS", "ORF", "blocked_reading_frame"):
                biotype = "protein_coding"
            elif g["type"] in RNA_TYPES:
                biotype = g["type"]
            else:
                biotype = g["type"]
        # product / description
        product = a.get("product", "")
        if not product:
            for c in gi["children"] + gi["cds"]:
                if c["attrs"].get("product"):
                    product = c["attrs"]["product"]; break
        if not product:
            product = a.get("description", "") or a.get("Note", "")
        name = a.get("gene") or a.get("Name") or ""
        if name == gid or name == a.get("ID"):
            name = a.get("gene", "") if a.get("gene") != gid else ""
        alias = a.get("Alias", "") or a.get("gene_synonym", "")
        old = a.get("old_locus_tag", "")
        # intervals
        if mode == "euk":
            parts = gi["exons"] or gi["cds"]
            iv = merge_intervals([(p["start"], p["end"]) for p in parts]) if parts else [[g["start"], g["end"]]]
        else:
            iv = [[g["start"], g["end"]]]
        length = sum(e - s + 1 for s, e in iv)
        for s, e in iv:
            saf.append((gid, g["seqid"], s, e, g["strand"]))
        if mode == "euk":
            # GTF for STAR: one transcript per child RNA (exons grouped by parent)
            by_tx = collections.defaultdict(list)
            for ex in gi["exons"]:
                by_tx[ex["parents"][0]].append(ex)
            for tx, exs in by_tx.items():
                for ex in sorted(exs, key=lambda r: r["start"]):
                    gtf.append(f'{ex["seqid"]}\tGFF\texon\t{ex["start"]}\t{ex["end"]}\t.\t{ex["strand"]}\t.\t'
                               f'gene_id "{gid}"; transcript_id "{clean_id(tx)}";')
        rows.append(dict(gene_id=gid, species=species, name=name, old_locus_tag=old, alias=alias,
                         biotype=biotype, product=product, chr=g["seqid"], start=g["start"], end=g["end"],
                         strand=g["strand"], length=length))
    return rows, saf, gtf


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--euk-gff", required=True); ap.add_argument("--euk-fna", required=True)
    ap.add_argument("--euk-name", default="Calb")
    ap.add_argument("--prok-gff", required=True); ap.add_argument("--prok-fna", required=True)
    ap.add_argument("--prok-name", default="Saur")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    euk_ids, prok_ids = fasta_ids(a.euk_fna), fasta_ids(a.prok_fna)
    dup = set(euk_ids) & set(prok_ids)
    if dup:
        sys.exit(f"[ERROR] contig names shared by both genomes: {sorted(dup)[:5]}")
    with open(os.path.join(a.outdir, "contig_species.tsv"), "w") as fh:
        for c in euk_ids:
            fh.write(f"{c}\t{a.euk_name}\n")
        for c in prok_ids:
            fh.write(f"{c}\t{a.prok_name}\n")

    all_rows, all_saf, gtf = [], [], []
    for gff, fna_ids, sp, mode in ((a.euk_gff, euk_ids, a.euk_name, "euk"), (a.prok_gff, prok_ids, a.prok_name, "prok")):
        gff_seqids = set()
        with opener(gff) as fh:
            for line in fh:
                if line.startswith("##FASTA"):
                    break
                if not line.startswith("#") and "\t" in line:
                    gff_seqids.add(line.split("\t", 1)[0])
        shared = gff_seqids & set(fna_ids)
        if not shared:
            sys.exit(f"[ERROR] {sp}: no sequence names in common between GFF ({sorted(gff_seqids)[:3]}) "
                     f"and FASTA ({fna_ids[:3]}). Use a GFF from the same source/assembly as the FASTA "
                     f"(00_get_annotations.sh downloads a matched pair from NCBI).")
        missing = set(fna_ids) - gff_seqids
        if missing:
            print(f"[WARN] {sp}: {len(missing)} FASTA contigs have no annotation: {sorted(missing)[:5]}", file=sys.stderr)
        rows, saf, g = build(gff, sp, mode, set(fna_ids))
        print(f"[INFO] {sp}: {len(rows)} genes, {len(saf)} SAF intervals "
              f"({collections.Counter(r['biotype'] for r in rows).most_common(6)})", file=sys.stderr)
        all_rows += rows; all_saf += saf; gtf += g

    ids = collections.Counter(r["gene_id"] for r in all_rows)
    clash = [k for k, v in ids.items() if v > 1]
    if clash:
        sys.exit(f"[ERROR] gene IDs shared between species: {clash[:5]}")

    with open(os.path.join(a.outdir, "combined.saf"), "w") as fh:
        fh.write("GeneID\tChr\tStart\tEnd\tStrand\n")
        for r in all_saf:
            fh.write("\t".join(map(str, r)) + "\n")
    cols = ["gene_id", "species", "name", "old_locus_tag", "alias", "biotype", "product", "chr", "start", "end",
            "strand", "length"]
    with open(os.path.join(a.outdir, "gene_annotation.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in all_rows:
            fh.write("\t".join(str(r[c]).replace("\t", " ") for c in cols) + "\n")
    with open(os.path.join(a.outdir, "calb_exons.gtf"), "w") as fh:
        fh.write("\n".join(gtf) + ("\n" if gtf else ""))


if __name__ == "__main__":
    main()
