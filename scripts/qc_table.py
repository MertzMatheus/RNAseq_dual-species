#!/usr/bin/env python3
"""
Aggregate QC metrics of every pipeline step into one table (one row per sample).

Steps covered: file integrity (MD5) -> raw reads (fastp 'before') -> FastQC -> trimming (fastp 'after')
-> alignment (STAR) -> species split (samtools) -> counting (featureCounts) -> gene detection.

Usage: qc_table.py --outdir /path/to/analysis --samplesheet samplesheet.tsv
Writes: <outdir>/qc/QC_summary.xlsx and <outdir>/qc/QC_per_sample.tsv
"""
import argparse, json, os, re, zipfile
import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import PatternFill, Font, Alignment
from openpyxl.utils import get_column_letter

# ---- thresholds used for the automatic flags (edit if needed) ----
MIN_RETAINED = 80        # % read pairs kept after trimming
MIN_UNIQUE_WARN = 70     # % uniquely mapped (STAR)
MIN_UNIQUE_FAIL = 50
CONTAM_WARN = 1.0        # % of uniquely mapped pairs on the UNEXPECTED species (single-species biofilms)
CONTAM_FAIL = 5.0
MIN_ASSIGNED = {"Calb": 5e6, "Saur": 2e6}   # assigned non-rRNA pairs recommended for DE
MAX_RRNA = 5             # % of assigned pairs on annotated rRNA genes AFTER SortMeRNA (residual)
MAX_RRNA_TOTAL = 30      # % of read pairs removed as rRNA by SortMeRNA (depletion efficiency)
SPECIES = ["Calb", "Saur"]
SPNAME = {"Calb": "C. albicans", "Saur": "S. aureus"}


def read_star(path):
    d = {}
    if not os.path.exists(path):
        return d
    for line in open(path):
        if "|" in line:
            k, v = [x.strip() for x in line.split("|", 1)]
            d[k] = v.rstrip("%")
    return d


def fastqc_fails(zips):
    fails, warns = [], []
    for z in zips:
        if not os.path.exists(z):
            continue
        with zipfile.ZipFile(z) as zf:
            name = [n for n in zf.namelist() if n.endswith("summary.txt")][0]
            for line in zf.read(name).decode().splitlines():
                st, mod, _ = line.split("\t")
                tag = ("R1" if re.search(r"(_1|_R1)[._]", os.path.basename(z)) else "R2") + ":" + mod
                (fails if st == "FAIL" else warns if st == "WARN" else []).append(tag)
    return "; ".join(sorted(set(fails))) or "none", len(set(warns))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--min-count", type=int, default=10)
    a = ap.parse_args()
    O = a.outdir
    ss = pd.read_csv(a.samplesheet, sep="\t", dtype=str)
    md5 = pd.read_csv(f"{O}/qc/md5_check.tsv", sep="\t") if os.path.exists(f"{O}/qc/md5_check.tsv") else None
    ann = pd.read_csv(f"{O}/reference/gene_annotation.tsv", sep="\t", dtype=str)
    counts = pd.read_csv(f"{O}/05_counts/counts_matrix.tsv", sep="\t", index_col=0)
    fcs = pd.read_csv(f"{O}/05_counts/featureCounts_counts.txt.summary", sep="\t", index_col=0)
    fcs.columns = [re.sub(r"\.Aligned\.sortedByCoord\.out\.bam$", "", os.path.basename(c)) for c in fcs.columns]
    strand = pd.read_csv(f"{O}/05_counts/strandedness.tsv", sep="\t", comment="#", index_col=0) \
        if os.path.exists(f"{O}/05_counts/strandedness.tsv") else None
    gene_sp = ann.set_index("gene_id")["species"]
    is_rrna = ann.set_index("gene_id")["biotype"].str.contains("rRNA", case=False, na=False)

    rows = []
    for _, r in ss.iterrows():
        s = r["sample"]; R = {"Sample": s, "Description": r["description"], "Group": r["group"],
                              "Expected species": r["species"]}
        # --- MD5 ---
        if md5 is not None:
            m = md5[md5["sample"] == s]
            R["MD5"] = "OK" if len(m) and (m["status"] == "OK").all() else (",".join(sorted(set(m["status"]))) or "NA")
            files = list(m["file"])
        else:
            R["MD5"] = "not run"; files = []
        # --- fastp ---
        jf = f"{O}/02_trimmed/{s}.fastp.json"
        if os.path.exists(jf):
            j = json.load(open(jf)); b = j["summary"]["before_filtering"]; aft = j["summary"]["after_filtering"]
            R["Raw read pairs"] = b["total_reads"] // 2
            R["Raw bases (Gb)"] = round(b["total_bases"] / 1e9, 2)
            R["Raw read length"] = f'{b["read1_mean_length"]}/{b["read2_mean_length"]}'
            R["Raw Q30 (%)"] = round(100 * b["q30_rate"], 2)
            R["Raw GC (%)"] = round(100 * b["gc_content"], 2)
            R["Duplication (%)"] = round(100 * j.get("duplication", {}).get("rate", float("nan")), 2)
            R["Adapter-trimmed reads (%)"] = round(100 * j.get("adapter_cutting", {}).get("adapter_trimmed_reads", 0)
                                                   / max(b["total_reads"], 1), 2)
            R["Clean read pairs"] = aft["total_reads"] // 2
            R["Pairs retained (%)"] = round(100 * aft["total_reads"] / max(b["total_reads"], 1), 2)
            R["Clean Q30 (%)"] = round(100 * aft["q30_rate"], 2)
            R["Clean mean length R1/R2"] = f'{aft["read1_mean_length"]}/{aft["read2_mean_length"]}'
        # --- SortMeRNA (rRNA removal) ---
        sl = f"{O}/02b_rRNA_removed/{s}.sortmerna.log"
        if os.path.exists(sl):
            txt = open(sl).read()
            tot = re.search(r"Total reads = (\d+)", txt); pas = re.search(r"Total reads passing E-value threshold = (\d+)", txt)
            if tot and pas:
                t, p_ = int(tot.group(1)), int(pas.group(1))
                R["rRNA pairs removed (%)"] = round(100 * p_ / max(t, 1), 2)
                R["Non-rRNA read pairs"] = (t - p_) // 2
        # --- FastQC ---
        zr = [f"{O}/01_fastqc_raw/" + re.sub(r"\.f(ast)?q\.gz$", "_fastqc.zip", f) for f in files] + \
             [f"{O}/01_fastqc_raw/{s}_R{i}.raw_fastqc.zip" for i in (1, 2)]
        zt = [f"{O}/03_fastqc_trimmed/{s}_R{i}.trim_fastqc.zip" for i in (1, 2)]
        R["FastQC raw: FAIL modules"], _ = fastqc_fails(zr)
        R["FastQC clean: FAIL modules"], _ = fastqc_fails(zt)
        # --- STAR ---
        st = read_star(f"{O}/04_alignment/{s}.Log.final.out")
        if st:
            R["STAR input pairs"] = int(st["Number of input reads"])
            R["Uniquely mapped (%)"] = float(st["Uniquely mapped reads %"])
            R["Multi-mapped (%)"] = float(st["% of reads mapped to multiple loci"])
            R["Too many loci (%)"] = float(st["% of reads mapped to too many loci"])
            R["Unmapped: too short (%)"] = float(st["% of reads unmapped: too short"])
            R["Unmapped: other/mismatch (%)"] = round(float(st["% of reads unmapped: other"]) +
                                                       float(st["% of reads unmapped: too many mismatches"]), 2)
            R["Mismatch rate per base (%)"] = float(st["Mismatch rate per base, %"])
            R["Spliced reads"] = int(st["Number of splices: Total"])
        # --- species split ---
        spf = f"{O}/04_alignment/{s}.species_counts.tsv"
        if os.path.exists(spf):
            sp = pd.read_csv(spf, sep="\t").set_index("species")
            tot = sp["pairs_unique"].sum()
            for k in SPECIES:
                R[f"{k} unique pairs"] = int(sp.loc[k, "pairs_unique"])
                R[f"{k} (% of unique)"] = round(100 * sp.loc[k, "pairs_unique"] / max(tot, 1), 2)
            for k in SPECIES:
                v = sp.loc[k, "mismatch_rate"]
                R[f"{k} mismatch rate (%)"] = round(100 * float(v), 3) if str(v) not in ("NA", "nan") else None
        # --- featureCounts ---
        if s in fcs.columns:
            c = fcs[s]; tot = c.sum()
            R["fC assigned (%)"] = round(100 * c.get("Assigned", 0) / max(tot, 1), 2)
            R["fC no feature (%)"] = round(100 * c.get("Unassigned_NoFeatures", 0) / max(tot, 1), 2)
            R["fC ambiguous (%)"] = round(100 * c.get("Unassigned_Ambiguity", 0) / max(tot, 1), 2)
            R["fC multi-map/low MAPQ (%)"] = round(100 * (c.get("Unassigned_MultiMapping", 0) +
                                                           c.get("Unassigned_MappingQuality", 0)) / max(tot, 1), 2)
        if s in counts.columns:
            col = counts[s]
            for k in SPECIES:
                idx = gene_sp.index[gene_sp == k]
                cs = col.reindex(idx).fillna(0)
                rr = cs[is_rrna.reindex(idx).fillna(False).values].sum()
                R[f"{k} assigned pairs"] = int(cs.sum())
                R[f"{k} residual rRNA (% of assigned)"] = round(100 * rr / max(cs.sum(), 1), 2)
                R[f"{k} genes detected (>= {a.min_count})"] = int((cs >= a.min_count).sum())
                R[f"{k} genes total"] = len(idx)
        if strand is not None and s in strand.index:
            R["Forward-strand fraction"] = strand.loc[s, "frac_forward"]

        # --- flags ---
        flags, level = [], 0
        def flag(msg, lv):
            nonlocal level; flags.append(msg); level = max(level, lv)
        if R.get("MD5") not in ("OK", "not run"): flag("MD5 problem", 2)
        if R.get("Pairs retained (%)", 100) < MIN_RETAINED: flag(f"<{MIN_RETAINED}% pairs kept after trimming", 1)
        u = R.get("Uniquely mapped (%)")
        if u is not None and u < MIN_UNIQUE_FAIL: flag(f"unique mapping {u}%", 2)
        elif u is not None and u < MIN_UNIQUE_WARN: flag(f"unique mapping {u}%", 1)
        exp = r["species"].split("+")
        for k in SPECIES:
            pc = R.get(f"{k} (% of unique)")
            if k not in exp and pc is not None:
                if pc > CONTAM_FAIL: flag(f"{SPNAME[k]} reads {pc}% (unexpected species)", 2)
                elif pc > CONTAM_WARN: flag(f"{SPNAME[k]} reads {pc}% (unexpected species)", 1)
            if k in exp and f"{k} assigned pairs" in R:
                nonr = R[f"{k} assigned pairs"] * (1 - R[f"{k} residual rRNA (% of assigned)"] / 100)
                if nonr < MIN_ASSIGNED[k]:
                    flag(f"{SPNAME[k]}: {nonr/1e6:.1f} M non-rRNA assigned pairs (<{MIN_ASSIGNED[k]/1e6:.0f} M)", 1)
                if R[f"{k} residual rRNA (% of assigned)"] > MAX_RRNA:
                    flag(f"{SPNAME[k]} residual rRNA {R[f'{k} residual rRNA (% of assigned)']}%", 1)
        if R.get("rRNA pairs removed (%)", 0) > MAX_RRNA_TOTAL:
            flag(f"rRNA {R['rRNA pairs removed (%)']}% of pairs (removed; check depletion)", 1)
        R["QC status"] = ["PASS", "WARN", "FAIL"][level]
        R["QC notes"] = "; ".join(flags)
        rows.append(R)

    df = pd.DataFrame(rows)
    os.makedirs(f"{O}/qc", exist_ok=True)
    df.to_csv(f"{O}/qc/QC_per_sample.tsv", sep="\t", index=False)

    # group summary (mean ± sd of numeric columns)
    num = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
    g = df.groupby("Group", sort=False)[num]
    mean, sd = g.mean(), g.std()
    gs = mean.copy().astype(object)
    for c in num:
        big = mean[c].abs().max() >= 1000
        gs[c] = [(f"{m/1e6:.2f} ± {s_/1e6:.2f} M" if big else f"{m:.2f} ± {s_:.2f}") for m, s_ in zip(mean[c], sd[c])]
    gs.insert(0, "n", df.groupby("Group", sort=False).size())
    gs = gs.reset_index()

    desc = pd.DataFrame([
        ("MD5", "Integrity check of delivered FASTQ files against MD5.txt"),
        ("Raw read pairs / bases / Q30 / GC", "fastp, before filtering"),
        ("rRNA pairs removed (%) / Non-rRNA read pairs", "SortMeRNA (SILVA + Rfam rRNA databases): pairs in which either mate matched rRNA were removed before alignment. Measures rRNA-depletion efficiency, including multi-copy rRNA operons."),
        ("Duplication (%)", "fastp estimate (high values are expected for highly expressed bacterial genes)"),
        ("Clean read pairs / Pairs retained", "fastp after adapter, poly-G, quality (4-bp window, Q20) and length (36 bp) filtering"),
        ("FastQC FAIL modules", "FastQC modules flagged FAIL (per-base sequence content and duplication FAILs are typical of RNA-seq)"),
        ("Uniquely mapped / Multi-mapped / Unmapped", "STAR, % of input read pairs, combined C. albicans + S. aureus genome"),
        ("Mismatch rate per base (%)", "STAR, all species together"),
        ("Calb / Saur unique pairs, % of unique", "Uniquely mapped read pairs (MAPQ 255) on each species' contigs"),
        ("Calb / Saur mismatch rate (%)", "samtools stats error rate on each genome; for S. aureus it reflects divergence of the isolate from USA300 FPR3757"),
        ("fC assigned / no feature / ambiguous", "featureCounts read-pair assignment (% of pairs presented)"),
        ("Calb / Saur assigned pairs", "Read pairs counted on genes of each species"),
        ("residual rRNA (% of assigned)", "Counted pairs still falling on annotated rRNA genes after SortMeRNA (should be ~0)"),
        ("genes detected", f"Genes with >= {a.min_count} counted pairs"),
        ("Forward-strand fraction", "Assigned(forward)/(forward+reverse): ~0.5 unstranded, ~0 reverse-stranded (dUTP), ~1 forward"),
        ("QC status", f"FAIL: MD5 problem, unique mapping <{MIN_UNIQUE_FAIL}%, or >{CONTAM_FAIL}% reads from the unexpected species. "
                      f"WARN: unique mapping <{MIN_UNIQUE_WARN}%, <{MIN_RETAINED}% pairs retained, >{CONTAM_WARN}% unexpected species, "
                      f"rRNA removed >{MAX_RRNA_TOTAL}%, residual rRNA >{MAX_RRNA}%, or fewer non-rRNA assigned pairs than recommended (C. albicans {MIN_ASSIGNED['Calb']/1e6:.0f} M, S. aureus {MIN_ASSIGNED['Saur']/1e6:.0f} M)."),
    ], columns=["Metric", "Definition"])

    xl = f"{O}/qc/QC_summary.xlsx"
    with pd.ExcelWriter(xl, engine="openpyxl") as w:
        df.to_excel(w, sheet_name="Per_sample_QC", index=False)
        gs.to_excel(w, sheet_name="Group_summary", index=False)
        desc.to_excel(w, sheet_name="Definitions", index=False)
        if md5 is not None: md5.to_excel(w, sheet_name="MD5", index=False)
        if strand is not None: strand.reset_index().to_excel(w, sheet_name="Strandedness", index=False)

    # styling
    wb = load_workbook(xl)
    step_col = [("MD5", "D9D9D9"), ("Raw", "DDEBF7"), ("Duplication", "DDEBF7"), ("Adapter", "DDEBF7"),
                ("Clean", "E2EFDA"), ("Pairs retained", "E2EFDA"), ("rRNA pairs", "D9E1F2"), ("Non-rRNA", "D9E1F2"), ("FastQC", "FFF2CC"), ("STAR", "FCE4D6"),
                ("Uniquely", "FCE4D6"), ("Multi", "FCE4D6"), ("Too many", "FCE4D6"), ("Unmapped", "FCE4D6"),
                ("Mismatch", "FCE4D6"), ("Spliced", "FCE4D6"), ("Calb", "EDE1F5"), ("Saur", "D0F0EF"),
                ("fC", "F2F2F2"), ("Forward", "F2F2F2"), ("QC", "FFFFFF")]
    fills = {"PASS": "C6EFCE", "WARN": "FFEB9C", "FAIL": "FFC7CE"}
    for ws in wb.worksheets:
        for cell in ws[1]:
            cell.font = Font(bold=True); cell.alignment = Alignment(wrap_text=True, vertical="top")
            for key, colr in step_col:
                if str(cell.value).startswith(key):
                    cell.fill = PatternFill("solid", fgColor=colr); break
        ws.freeze_panes = "B2"
        for i, col in enumerate(ws.columns, 1):
            ws.column_dimensions[get_column_letter(i)].width = 14 if ws.title != "Definitions" else (36 if i == 1 else 120)
        ws.row_dimensions[1].height = 45
    ws = wb["Per_sample_QC"]
    hdr = [c.value for c in ws[1]]
    if "QC status" in hdr:
        ci = hdr.index("QC status") + 1
        for row in ws.iter_rows(min_row=2):
            v = row[ci - 1].value
            if v in fills:
                row[ci - 1].fill = PatternFill("solid", fgColor=fills[v])
        ws.column_dimensions[get_column_letter(ci + 1)].width = 80
    wb.save(xl)
    print(f"[qc_table] wrote {xl}")
    print(df[["Sample", "Group", "QC status", "QC notes"]].to_string(index=False))


if __name__ == "__main__":
    main()
