#!/usr/bin/env python3.12
"""Build a synthetic genome + transcripts with known truth.

Gene: 4 exons, 3 introns (GT-AG, GC-AG, GT-AG).
contigA = gene on + strand; contigB = revcomp(contigA) so the gene is on - strand.
Transcripts:
  clean_fwd   : exact exon concat
  err_fwd     : 3 substitutions, one 2bp insertion (extra in transcript),
                one 2bp deletion (missing from transcript), 30bp soft-clip tail
  clean_rc / err_rc : reverse complements of the above
Truth is written as JSON: genome exon coords per contig and transcript coords
per exon (in the transcript's own forward orientation, 1-based inclusive).
"""
import json, random, sys
random.seed(7)
COMP = str.maketrans("ACGTacgt", "TGCAtgca")
def rc(s): return s.translate(COMP)[::-1]
def rnd(n): return "".join(random.choice("ACGT") for _ in range(n))

# exons lengths and introns (motif embedded)
exlen = [180, 220, 160, 240]
introns = [("GT", 300, "AG"), ("GC", 250, "AG"), ("GT", 400, "AG")]
left = rnd(200)
parts = [left]
pos = len(left)
exons_A = []
for i, L in enumerate(exlen):
    ex = rnd(L)
    parts.append(ex)
    exons_A.append((pos + 1, pos + L))
    pos += L
    if i < len(introns):
        d, n, a = introns[i]
        body = d + rnd(n - 4) + a
        # avoid accidental splice motif near ends: fine for a test
        parts.append(body)
        pos += n
parts.append(rnd(200))
contigA = "".join(parts)
tx_clean = "".join(contigA[s - 1:e] for s, e in exons_A)
tlen = len(tx_clean)

# transcript coords of exons (clean)
tcoords_clean = []
p = 0
for L in exlen:
    tcoords_clean.append((p + 1, p + L))
    p += L

# introduce errors into the transcript; record how transcript coords shift
# substitutions at tx positions (1-based) 50, 400, 700 (all inside exons, away from boundaries)
# insertion of 2bp after tx pos 250 (exon2), deletion of tx 600-601 (exon3 region: 401-560? no)
# exon layout: e1 1-180, e2 181-400, e3 401-560, e4 561-800
tx = list(tx_clean)
def sub(i):
    b = tx[i - 1]; tx[i - 1] = {"A": "C", "C": "G", "G": "T", "T": "A"}[b]
for s in (50, 300, 700): sub(s)
# deletion: remove tx 480-481 (exon3)
del tx[479:481]
# insertion: insert "TT" after original tx pos 250 (exon2) -> after deletion index unchanged since 250<480
tx.insert(250, "TT")
tx_err = "".join(tx)
tail = rnd(30)
tx_err_clip = tx_err + tail
# truth coords in err transcript: exon1 1-180 ; exon2 181-402 (+2 ins); exon3 403-560 (-2 del); exon4 561-800 ; total 800, clip 801-830
tcoords_err = [(1, 180), (181, 402), (403, 560), (561, 800)]
assert len(tx_err) == 800

contigB = rc(contigA)
G = len(contigA)
exons_B = [(G - e + 1, G - s + 1) for s, e in exons_A][::-1]  # ascending on contigB

with open("genome.fa", "w") as f:
    f.write(f">contigA\n{contigA}\n>contigB\n{contigB}\n")
with open("tx.fa", "w") as f:
    f.write(f">clean_fwd\n{tx_clean}\n>clean_rc\n{rc(tx_clean)}\n")
    f.write(f">err_fwd\n{tx_err_clip}\n>err_rc\n{rc(tx_err_clip)}\n")
truth = {
    "genome_exons": {"contigA": exons_A, "contigB": exons_B},
    "tx_len": {"clean_fwd": tlen, "clean_rc": tlen, "err_fwd": 830, "err_rc": 830},
    # transcript coords per exon, ascending in genome order on contigA, in the
    # transcript's OWN forward coordinates
    "tx_coords_fwd": {"clean": tcoords_clean, "err": tcoords_err},
}
json.dump(truth, open("truth.json", "w"), indent=1)
print("genome", G, "tx", tlen, "exonsA", exons_A, "exonsB", exons_B)
print("tx_err_coords", tcoords_err)
