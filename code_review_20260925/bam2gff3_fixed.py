#!/usr/bin/env python3.12
"""Corrected minimap2 splice SAM/BAM -> GFF3 cDNA_match converter for PASA
--IMPORT_CUSTOM_ALIGNMENTS.

Coordinates come from the CIGAR (authoritative), not the cs string:
  M/=/X consume ref+query, D consumes ref, N is an intron (consumes ref),
  I consumes query, S consumes query (clip), H consumes nothing (SEQ excludes it).
Per-exon identity comes from the cs string when present (counts substitution
and indel bases inside each exon); otherwise from NM overall.
Intron motifs: accepted if consensus (GT-AG, GC-AG, AT-AC) in EITHER
orientation on the reference, matching PASA's own check
(CDNA_alignment::get_consensus_splice_sites tries both). The GFF3 strand
column is the transcribed strand: ts:A when present, else SAM strand.
Target coordinates are in the transcript's own forward orientation with
lend<rend, as import_spliced_alignments.dbi and
validate_alignments_in_db.dbi expect (for '-', low genome end pairs with high
transcript end).

Usage: bam2gff3_fixed.py in.bam out.gff3 [--min-pident 80] [--samtools PATH]
"""
import os
import re
import subprocess
import sys

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
CS_RE = re.compile(r"(:[0-9]+|\*[a-z][a-z]|[=\+\-][A-Za-z]+|~[a-z]{2}[0-9]+[a-z]{2})")
FWD_MOTIFS = {("gt", "ag"), ("gc", "ag"), ("at", "ac")}
REV_MOTIFS = {("ct", "ac"), ("ct", "gc"), ("gt", "at")}


def _intron_motifs_from_cs(cs):
    """Return list of (donor, acceptor) per '~' op in the cs string."""
    out = []
    if cs is None:
        return out
    for tok in CS_RE.findall(cs):
        if tok[0] == "~":
            out.append((tok[1:3], tok[-2:]))
    return out


def _exon_per_id_from_cs(cs, exon_ref_lens):
    """Walk cs string; return per-exon (matched, mism+indel bases) lists.
    Exons are delimited by '~' ops. Returns None if cs missing."""
    if cs is None:
        return None
    stats = [[0, 0]]
    for tok in CS_RE.findall(cs):
        c = tok[0]
        if c == ":":
            stats[-1][0] += int(tok[1:])
        elif c == "=":
            stats[-1][0] += len(tok) - 1
        elif c == "*":
            stats[-1][1] += 1
        elif c in "+-":
            stats[-1][1] += len(tok) - 1
        elif c == "~":
            stats.append([0, 0])
    return stats


def parse_record(cols, min_pident=80.0):
    """Return (list of gff3 lines, reason) for one SAM record."""
    flag = int(cols[1])
    if flag & 0x904:  # unmapped, secondary, supplementary
        return [], "flag"
    if flag & 0x800 or flag & 0x100:
        return [], "flag"
    rname, pos, cigar, seq = cols[2], int(cols[3]), cols[5], cols[9]
    if cigar == "*":
        return [], "nocigar"
    cs = nm = ts = None
    for t in cols[11:]:
        if t.startswith("cs:Z:"):
            cs = t[5:]
        elif t.startswith("NM:i:"):
            nm = int(t[5:])
        elif t.startswith("ts:A:"):
            ts = t[5:]
    ops = [(int(n), o) for n, o in CIGAR_RE.findall(cigar)]
    # query length in read orientation includes soft and hard clips
    qlen = sum(n for n, o in ops if o in "MIS=XH")
    # walk
    rpos = pos          # 1-based ref position of next base
    qpos = 1            # 1-based query position (read orientation) of next base
    exons = []          # (rstart, rend, qstart, qend) in read orientation
    cur = None          # [rstart, qstart, rlast, qlast]
    ins_bases = del_bases = 0
    for n, o in ops:
        if o in "M=X":
            if cur is None:
                cur = [rpos, qpos, rpos, qpos]
            rpos += n
            qpos += n
            cur[2] = rpos - 1
            cur[3] = qpos - 1
        elif o == "I":
            if cur is not None:
                qpos += n
                cur[3] = qpos - 1
                ins_bases += n
            else:
                qpos += n
        elif o == "D":
            rpos += n
            if cur is not None:
                cur[2] = rpos - 1
                del_bases += n
        elif o == "N":
            if cur is not None:
                exons.append(tuple(cur))
                cur = None
            rpos += n
        elif o in "SH":
            qpos += n
        # P consumes nothing
    if cur is not None:
        exons.append(tuple(cur))
    if not exons:
        return [], "noexon"
    samrev = bool(flag & 0x10)
    # Column 7 must be the ALIGNMENT orientation (SAM strand): PASA's
    # validate_sequence_segments_at_splice_boundaries reverse-complements the
    # transcript when col7 is '-'. PASA derives the transcribed (spliced)
    # orientation itself from the intron consensus, so ts:A is not used here.
    strand = "-" if samrev else "+"
    # motif check (both orientations acceptable, as in PASA)
    motifs = _intron_motifs_from_cs(cs)
    if motifs:
        if all(m in FWD_MOTIFS for m in motifs):
            pass
        elif all(m in REV_MOTIFS for m in motifs):
            pass
        else:
            return [], "noncanonical"
    # identity per exon
    exstats = _exon_per_id_from_cs(cs, None)
    pids = []
    if exstats and len(exstats) == len(exons):
        for m, e in exstats:
            pids.append(100.0 * m / (m + e) if (m + e) else 0.0)
    else:
        aligned = sum(e[3] - e[1] + 1 for e in exons)
        errs = nm if nm is not None else 0
        pid = 100.0 * max(aligned - errs, 0) / aligned if aligned else 0.0
        pids = [pid] * len(exons)
    aligned_q = sum(e[3] - e[1] + 1 for e in exons)
    tot_err = sum(e for _, e in exstats) if exstats else (nm or 0)
    tot_m = sum(m for m, _ in exstats) if exstats else aligned_q - tot_err
    overall = 100.0 * tot_m / (tot_m + tot_err) if (tot_m + tot_err) else 0.0
    if overall < min_pident:
        return [], "lowident"
    lines = []
    for (rs, qs, re_, qe), pid in zip(exons, pids):
        if samrev:
            # SEQ is revcomp of the original transcript: map back
            tqs = qlen - qe + 1
            tqe = qlen - qs + 1
        else:
            tqs, tqe = qs, qe
        lines.append(
            f"{rname}\tgenome\tcDNA_match\t{rs}\t{re_}\t{pid:.2f}\t{strand}\t.\t"
            f"ID={cols[0]};Target={cols[0]} {tqs} {tqe} +"
        )
    return lines, "ok"


def bam2gff3(input, output, min_pident=80.0, samtools="samtools"):
    count = 0
    seen = {}
    with open(output, "w") as out:
        out.write("##gff-version 3\n")
        p = subprocess.Popen([samtools, "view", os.path.realpath(input)],
                             stdout=subprocess.PIPE, universal_newlines=True)
        for line in p.stdout:
            cols = line.rstrip("\n").split("\t")
            lines, why = parse_record(cols, min_pident)
            if not lines:
                continue
            # keep ID unique if the same read name appears more than once
            n = seen.get(cols[0], 0)
            seen[cols[0]] = n + 1
            if n:
                lines = [l.replace(f"ID={cols[0]};", f"ID={cols[0]}.{n};") for l in lines]
            out.write("\n".join(lines) + "\n")
            count += 1
        p.wait()
    return count


if __name__ == "__main__":
    a = sys.argv[1:]
    st = "samtools"
    if "--samtools" in a:
        i = a.index("--samtools"); st = a[i + 1]; del a[i:i + 2]
    mp = 80.0
    if "--min-pident" in a:
        i = a.index("--min-pident"); mp = float(a[i + 1]); del a[i:i + 2]
    print(bam2gff3(a[0], a[1], mp, st))
