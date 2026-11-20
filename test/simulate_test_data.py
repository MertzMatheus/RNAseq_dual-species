#!/usr/bin/env python3
"""
Simulate a tiny dual RNA-seq project (same layout as a real run) to test the pipeline end-to-end.
  - mini C. albicans genome (2 chromosomes, NCBI-style GFF3 gene>mRNA>exon/CDS, some introns, rRNA/tRNA)
  - mini S. aureus genome (chromosome + plasmid, NCBI prokaryote GFF3 gene>CDS with old_locus_tag, rRNA)
  - 27 samples, reverse-stranded (dUTP) 150-bp paired reads, known fold changes
  - isolate SNP divergence from the reference, one gene absent in CI5, LR5 split in two lanes,
    LR2 lightly contaminated with S. aureus reads
Usage: simulate_test_data.py <project_dir> [pairs_per_sample] [rRNA.fasta]
  rRNA.fasta (optional): records named calb_lsu and saur_ssu; the simulated rRNA genes then carry
  real rRNA sequence, so SortMeRNA can detect them.
"""
import gzip, hashlib, os, random, sys
import numpy as np

root = sys.argv[1]; NPAIRS = int(sys.argv[2]) if len(sys.argv) > 2 else 60000
REAL_RRNA = {}
if len(sys.argv) > 3:
    nm = None
    for line in open(sys.argv[3]):
        if line.startswith(">"): nm = line[1:].split()[0]; REAL_RRNA[nm] = ""
        else: REAL_RRNA[nm] += line.strip().upper().replace("U", "T")
rng = np.random.default_rng(7); random.seed(7)
G = os.path.join(root, "genomes"); R = os.path.join(root, "raw_data")
os.makedirs(G, exist_ok=True); os.makedirs(R, exist_ok=True)
COMP = str.maketrans("ACGT", "TGCA")
rc = lambda s: s.translate(COMP)[::-1]
rand_seq = lambda n, gc=0.35: "".join(rng.choice(list("ACGT"), n, p=[(1-gc)/2, gc/2, gc/2, (1-gc)/2]))


def write_fa(path, recs):
    with open(path, "w") as f:
        for n, s in recs:
            f.write(f">{n} simulated\n" + "\n".join(s[i:i+80] for i in range(0, len(s), 80)) + "\n")


# ---------------- C. albicans ----------------
calb_chr = {"NC_032089.1": rand_seq(220000, .33), "NC_032090.1": rand_seq(180000, .33)}
calb_genes, gff = [], ["##gff-version 3"]
k = 0
for ch, seq in calb_chr.items():
    gff.append(f"{ch}\tRefSeq\tregion\t1\t{len(seq)}\t.\t+\t.\tID={ch}:1..{len(seq)};Dbxref=taxon:237561")
    pos = 500
    while pos < len(seq) - 4000:
        k += 1; lt = f"CAALFM_C{k:05d}0WA"; strand = "+" if rng.random() < .5 else "-"
        kind = "rRNA" if k in (5, 6) else ("tRNA" if k in (10, 11, 12) else "mRNA")
        L = 1500 if kind == "rRNA" else (72 if kind == "tRNA" else int(rng.integers(600, 2500)))
        if kind == "rRNA" and "calb_lsu" in REAL_RRNA:
            r = REAL_RRNA["calb_lsu"][:L + 1]
            seq = seq[:pos - 1] + r + seq[pos - 1 + len(r):]; calb_chr[ch] = seq
        intron = kind == "mRNA" and rng.random() < .15
        exons = [(pos, pos + L // 3), (pos + L // 3 + 80, pos + L + 80)] if intron else [(pos, pos + L)]
        gs, ge = exons[0][0], exons[-1][1]
        name = f"GEN{k}" if rng.random() < .5 else ""
        nm = f";Name={name};gene={name}" if name else f";Name={lt}"
        bt = {"mRNA": "protein_coding", "rRNA": "rRNA", "tRNA": "tRNA"}[kind]
        gff.append(f"{ch}\tRefSeq\tgene\t{gs}\t{ge}\t.\t{strand}\t.\tID=gene-{lt};Dbxref=GeneID:{k}{nm};gbkey=Gene;gene_biotype={bt};locus_tag={lt}")
        tid = f"rna-XM_{k:06d}.1"
        gff.append(f"{ch}\tRefSeq\t{kind}\t{gs}\t{ge}\t.\t{strand}\t.\tID={tid};Parent=gene-{lt};gbkey={kind};locus_tag={lt};product=hypothetical protein {k}%2C putative")
        for i, (s, e) in enumerate(exons):
            gff.append(f"{ch}\tRefSeq\texon\t{s}\t{e}\t.\t{strand}\t.\tID=exon-XM_{k:06d}.1-{i+1};Parent={tid};gbkey={kind};locus_tag={lt}")
            if kind == "mRNA":
                gff.append(f"{ch}\tRefSeq\tCDS\t{s}\t{e}\t.\t{strand}\t0\tID=cds-XP_{k:06d}.1;Parent={tid};product=hypothetical protein {k}%2C putative;locus_tag={lt}")
        tx = "".join(seq[s-1:e] for s, e in exons)
        calb_genes.append(dict(id=lt, chr=ch, strand=strand, tx=tx if strand == "+" else rc(tx), kind=kind))
        pos = ge + int(rng.integers(300, 1500))
write_fa(f"{G}/C_albicans_SC5314.fna", calb_chr.items())
open(f"{G}/C_albicans_SC5314.gff", "w").write("\n".join(gff) + "\n")

# ---------------- S. aureus ----------------
saur_chr = {"NC_007793.1": rand_seq(300000, .33), "NC_007792.1": rand_seq(4000, .30)}
saur_genes, gff = [], ["##gff-version 3"]
k = 0
for ch, seq in saur_chr.items():
    gff.append(f"{ch}\tRefSeq\tregion\t1\t{len(seq)}\t.\t+\t.\tID={ch}:1..{len(seq)}")
    pos = 200
    while pos < len(seq) - 1500:
        k += 1; lt = f"SAUSA300_RS{k*5:05d}"; old = f"SAUSA300_{k:04d}"
        strand = "+" if rng.random() < .6 else "-"
        kind = "rRNA" if k in (3, 4) else "CDS"
        L = 1500 if kind == "rRNA" else int(rng.integers(300, 1500))
        bt = "rRNA" if kind == "rRNA" else "protein_coding"
        if kind == "rRNA" and "saur_ssu" in REAL_RRNA:
            r = REAL_RRNA["saur_ssu"][:L + 1]; L = len(r) - 1
            seq = seq[:pos - 1] + r + seq[pos - 1 + len(r):]; saur_chr[ch] = seq
        name = f"sau{k}" if rng.random() < .4 else ""
        nm = f";Name={name};gene={name}" if name else f";Name={lt}"
        gff.append(f"{ch}\tRefSeq\tgene\t{pos}\t{pos+L}\t.\t{strand}\t.\tID=gene-{lt}{nm};gbkey=Gene;gene_biotype={bt};locus_tag={lt};old_locus_tag={old}")
        gff.append(f"{ch}\tProtein Homology\t{kind}\t{pos}\t{pos+L}\t.\t{strand}\t0\tID=cds-WP_{k:06d}.1;Parent=gene-{lt};product=protein {k};locus_tag={lt}")
        tx = seq[pos-1:pos+L]
        saur_genes.append(dict(id=lt, chr=ch, strand=strand, tx=tx if strand == "+" else rc(tx), kind=kind))
        pos += L + int(rng.integers(50, 400))
write_fa(f"{G}/S_aureus_USA300_FPR3757.fna", saur_chr.items())
open(f"{G}/S_aureus_USA300_FPR3757.gff", "w").write("\n".join(gff) + "\n")
print(f"C. albicans genes {len(calb_genes)}, S. aureus genes {len(saur_genes)}")

# ---------------- expression model ----------------
nC, nS = len(calb_genes), len(saur_genes)
baseC = rng.lognormal(3, 1.3, nC); baseS = rng.lognormal(3, 1.3, nS)
for i, g in enumerate(calb_genes):
    if g["kind"] == "rRNA": baseC[i] *= 30
for i, g in enumerate(saur_genes):
    if g["kind"] == "rRNA": baseS[i] *= 30
iso = ["JE2", "CI19", "S5", "CI5"]
# S. aureus: isolate effects (mono), C. albicans effect (dual), isolate-specific C. albicans effect
fS_iso = {i: np.where(rng.random(nS) < .12, rng.choice([-2, 2], nS), 0.0) for i in iso}; fS_iso["JE2"][:] = 0
fS_dual = np.where(rng.random(nS) < .10, rng.choice([-2, 2], nS), 0.0)
fS_dual_ci19 = np.where(rng.random(nS) < .05, rng.choice([-2.5, 2.5], nS), 0.0)
absent_in_CI5 = 20                                   # gene index absent from CI5
divergence = {"JE2": .001, "CI19": .01, "S5": .015, "CI5": .02}
# C. albicans: generic S. aureus effect + isolate-specific effects
fC_dual = np.where(rng.random(nC) < .10, rng.choice([-2, 2], nC), 0.0)
fC_iso = {i: np.where(rng.random(nC) < .05, rng.choice([-1.8, 1.8], nC), 0.0) for i in iso}; fC_iso["JE2"][:] = 0

def mutate(s, rate, seed):
    r = np.random.default_rng(seed); a = np.frombuffer(s.encode(), dtype="S1").copy()
    idx = np.where(r.random(len(a)) < rate)[0]
    for j in idx: a[j] = r.choice([b for b in b"ACGT" if b != a[j][0]]).tobytes()
    return a.tobytes().decode()
saur_tx = {i: [mutate(g["tx"], divergence[i], 100 + n) for n, g in enumerate(saur_genes)] for i in iso}

ADAPT = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
def frag_reads(tx, n, out1, out2, tag, rs):
    for _ in range(n):
        L = len(tx); fl = int(min(L, max(60, rs.normal(260, 70))))
        st = int(rs.integers(0, L - fl + 1)); f = tx[st:st+fl]
        r1 = rc(f)[:150]; r2 = f[:150]              # dUTP: read 1 antisense
        if len(r1) < 150: r1 = (r1 + ADAPT + "A" * 150)[:150]; r2 = (r2 + ADAPT.replace("CACAC", "TCGTA") + "A" * 150)[:150]
        q = "I" * 120 + "".join(rs.choice(list("?5+#"), 30))
        out1.append(f"@{tag}\n{r1}\n+\n{q}\n"); out2.append(f"@{tag}\n{r2}\n+\n{q}\n")

groups = [("none", "mono", "C"), *[(i, "mono", "S") for i in iso], *[(i, "dual", "CS") for i in iso]]
n = 0
for isol, cul, sp in groups:
    for rep in range(3):
        n += 1; s = f"LR{n}"; rs = np.random.default_rng(1000 + n)
        wC, wS = [], []
        if "C" in sp:
            lf = fC_dual + fC_iso[isol] if cul == "dual" else 0
            wC = baseC * 2.0 ** lf * rs.gamma(20, 1/20, nC)
        if "S" in sp:
            lf = fS_iso[isol] + (fS_dual + (fS_dual_ci19 if isol == "CI19" else 0) if cul == "dual" else 0)
            wS = baseS * 2.0 ** lf * rs.gamma(20, 1/20, nS)
            if isol == "CI5": wS[absent_in_CI5] = 0
        if s == "LR2":   # light contamination
            wS = baseS * 0.02
        fracC = .0 if "C" not in sp else (1.0 if "S" not in sp else .55)
        if s == "LR2": fracC = .98
        out1, out2 = [], []; c = 0
        for frac, w, genes, sel in ((fracC, wC, calb_genes, "C"), (1 - fracC, wS, saur_genes, "S")):
            if frac <= 0 or len(w) == 0: continue
            lens = np.array([len(g["tx"]) for g in genes]); p = w * lens; p /= p.sum()
            cnt = rs.multinomial(int(NPAIRS * frac), p)
            for gi, m in enumerate(cnt):
                if m == 0: continue
                tx = genes[gi]["tx"] if sel == "C" else saur_tx[isol if isol != "none" else "JE2"][gi]
                frag_reads(tx, m, out1, out2, f"{s}:{sel}:{gi}:{c}", rs); c += m
        perm = rs.permutation(len(out1)); out1 = [out1[i] for i in perm]; out2 = [out2[i] for i in perm]
        # add 2% random junk pairs (unmappable)
        for j in range(int(.02 * NPAIRS)):
            a = rand_seq(150, .5); out1.append(f"@junk{j}\n{a}\n+\n{'I'*150}\n"); out2.append(f"@junk{j}\n{rc(a)}\n+\n{'I'*150}\n")
        d = f"{R}/{s}"; os.makedirs(d, exist_ok=True)
        lanes = 2 if s == "LR5" else 1; md5 = []
        for ln in range(lanes):
            sl = slice(ln * len(out1) // lanes, (ln + 1) * len(out1) // lanes)
            for rd, data in ((1, out1), (2, out2)):
                fn = f"{s}_CKDL2300{n:02d}-1A_HYK3LDSX{ln+5}_L{ln+1}_{rd}.fq.gz"
                with gzip.open(f"{d}/{fn}", "wt", compresslevel=1) as fh: fh.writelines(data[sl])
                md5.append(f"{hashlib.md5(open(f'{d}/{fn}', 'rb').read()).hexdigest()}  {fn}")
        open(f"{d}/MD5.txt", "w").write("\n".join(md5) + "\n")
print("simulated", n, "samples in", R)
np.savetxt(f"{root}/truth_Saur_dual_effect.txt", fS_dual)
