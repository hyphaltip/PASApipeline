#!/usr/bin/env python3.12
"""Run ORIGINAL funannotate bam2gff3 (source pulled from library.py by ast)
and the FIXED converter on the same BAM/SAM; compare to truth; emulate the
PASA import + validate checks."""
import ast, json, os, subprocess, sys, types
SAMTOOLS = "/rhome/jstajich/projects/funannotate/funannotate-live/.pixi/envs/default/bin/samtools"
LIB = os.path.expanduser("~/projects/funannotate/funannotate-live/funannotate/library.py")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bam2gff3_fixed as fixed

# ---- load original function bodies verbatim
src = open(LIB).read()
tree = ast.parse(src)
ns = {"os": os}
def execute(cmd):
    cmd = [SAMTOOLS if c == "samtools" else c for c in cmd]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
    for l in p.stdout:
        yield l
    p.wait()
ns["execute"] = execute
for node in tree.body:
    if isinstance(node, ast.FunctionDef) and node.name in ("tokenizeString", "bam2gff3"):
        exec(compile(ast.Module([node], []), LIB, "exec"), ns)
orig_bam2gff3 = ns["bam2gff3"]

COMP = str.maketrans("ACGTacgt", "TGCAtgca")
def rc(s): return s.translate(COMP)[::-1]
def readfa(fn):
    d = {}; k = None
    for l in open(fn):
        l = l.strip()
        if l.startswith(">"): k = l[1:].split()[0]; d[k] = []
        else: d[k].append(l)
    return {k: "".join(v) for k, v in d.items()}

def parse_gff(fn):
    """return {ID: [(contig, s, e, pid, strand, ts, te)]}"""
    out = {}
    for l in open(fn):
        if l.startswith("#") or not l.strip(): continue
        c = l.rstrip("\n").split("\t")
        attrs = c[8]
        mid = attrs.split("ID=")[1].split(";")[0]
        tgt = attrs.split("Target=")[1].split(";")[0].split()
        out.setdefault(mid, []).append((c[0], int(c[3]), int(c[4]), float(c[5]), c[6], int(tgt[1]), int(tgt[2])))
    return out

def pasa_emulate(segs, genome, txseq, min_pct=90, min_id=95, nbp=3):
    """Emulate import_spliced_alignments.dbi + validate_alignments_in_db.dbi.
    segs: list of (contig, s, e, pid, strand, tl, tr) ; returns (ok, reason)"""
    orient = segs[0][4]
    contig = segs[0][0]
    G = genome[contig]
    T = txseq.upper()
    # CDNA_alignment::determine_alignment_attributes
    nts = sum(abs(tr - tl) + 1 for *_, tl, tr in segs)
    avg = sum(pid * (abs(tr - tl) + 1) for (_, _, _, pid, _, tl, tr) in segs) / nts
    pct = nts / len(T) * 100
    reasons = []
    # splice consensus, both orientations (identify_splice_junctions)
    segs_sorted = sorted(segs, key=lambda x: x[1])
    if len(segs) > 1:
        ok_any = False
        for o in ("+", "-"):
            good = True
            for a, b in zip(segs_sorted, segs_sorted[1:]):
                il, ir = a[2] + 1, b[1] - 1
                d, ac = G[il - 1:il + 1].upper(), G[ir - 2:ir].upper()
                pairs = {("GT", "AG"), ("GC", "AG"), ("AT", "AC")} if o == "+" else {("CT", "AC"), ("CT", "GC"), ("GT", "AT")}
                if (d, ac) not in pairs: good = False
            if good: ok_any = True
        if not ok_any: reasons.append("splice_consensus")
        # NUM_BP_PERFECT_SPLICE_BOUNDARY (uses aligned orient = gff col7)
        for i, (c, gl, gr, pid, st, tl, tr) in enumerate(segs_sorted):
            ml, mr = sorted((tl, tr))
            if i != 0:
                for j in range(nbp):
                    if orient == "+":
                        tc = T[ml + j - 1]; gc = G[gl + j - 1].upper()
                    else:
                        tc = rc(T[mr - j - 1]); gc = G[gl + j - 1].upper()
                    if tc != gc: reasons.append(f"left_boundary_seg{i}"); break
            if i != len(segs_sorted) - 1:
                for j in range(nbp):
                    if orient == "+":
                        tc = T[mr - j - 1]; gc = G[gr - j - 1].upper()
                    else:
                        tc = rc(T[ml + j - 1]); gc = G[gr - j - 1].upper()
                    if tc != gc: reasons.append(f"right_boundary_seg{i}"); break
    if pct < min_pct: reasons.append(f"pct_aligned={pct:.1f}<{min_pct}")
    if avg < min_id: reasons.append(f"avg_id={avg:.1f}<{min_id}")
    return (not reasons), ";".join(reasons) or "OK", pct, avg

def compare(label, gff, truth, genome, tx, orient_expect=None):
    got = parse_gff(gff)
    print(f"\n=== {label}: {gff}")
    for tid in tx:
        segs = got.get(tid)
        if not segs:
            print(f"  {tid:10s} NOT EMITTED by converter"); continue
        contig = segs[0][0]
        gex = [(s, e) for _, s, e, *_ in sorted(segs, key=lambda x: x[1])]
        truth_g = [tuple(x) for x in truth["genome_exons"][contig]]
        kind = "clean" if tid.startswith("clean") else "err"
        tc = [tuple(x) for x in truth["tx_coords_fwd"][kind]]
        if tid.endswith("_rc"):
            L = truth["tx_len"][tid]
            tc = [(L - b + 1, L - a + 1) for a, b in tc]
        # transcript coords must be paired with genome exons; on contigB genome
        # order is reversed relative to transcript order
        pairs = list(zip(sorted(truth_g), tc if contig == "contigA" else tc[::-1]))
        got_pairs = [(( s, e), (tl, te)) for _, s, e, _, _, tl, te in sorted(segs, key=lambda x: x[1])]
        okg = gex == sorted(truth_g)
        okt = [p[1] for p in got_pairs] == [p[1] for p in pairs]
        ok, why, pct, avg = pasa_emulate(segs, genome, tx[tid])
        print(f"  {tid:10s} strand={segs[0][4]} genome_exons_ok={okg} target_ok={okt} pid={segs[0][3]:.1f} pasa={why} (pct={pct:.1f},avg={avg:.1f})")
        if not okg or not okt:
            print("     got   :", got_pairs)
            print("     truth :", pairs)

if __name__ == "__main__":
    bam = sys.argv[1]
    truth = json.load(open("truth.json"))
    genome = readfa("genome.fa"); tx = readfa("tx.fa")
    n1 = orig_bam2gff3(bam, "orig.gff3")
    n2 = fixed.bam2gff3(bam, "fixed.gff3", samtools=SAMTOOLS)
    print("records emitted: original", n1, "fixed", n2)
    compare("ORIGINAL", "orig.gff3", truth, genome, tx)
    compare("FIXED", "fixed.gff3", truth, genome, tx)
