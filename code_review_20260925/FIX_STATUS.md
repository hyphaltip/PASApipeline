# Fix status: PASA / funannotate training-set work

Last updated 2026-09-25 by the REVIEW session (Claude Opus 5.5). The full reasoning and evidence for each entry are in the decision log at `/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/do_pasa_rust_vs_perl/DECISIONS.md`; the D-numbers below point there. The review itself is in `../CODE_REVIEW_20260925.md`. **Nothing is committed or merged yet.** Every change is in an uncommitted worktree, waiting for review and user approval.

## Where the code changes are

| Repo | Worktree | Branch | Owner |
|---|---|---|---|
| funannotate | `~/projects/funannotate/funannotate-live-bam2gff3` | `fix/bam2gff3-cigar` (base `41a2fd7`) | REVIEW |
| funannotate | `~/projects/funannotate/funannotate-live` (main working tree) | uncommitted | SELECT (the other Claude session) |
| PASApipeline | `../PASApipeline-r2` | `fix/unspliced-orient-clustering` (base `23d67c0`) | REVIEW |
| PASApipeline | this tree (`rust_optimize`) | review files only | REVIEW |

Changelogs: `PASApipeline-r2/Changelog.txt` (Unreleased section) and `funannotate-live-bam2gff3/CHANGELOG.md` (Unreleased → Fixed).

## Status by fix

| ID | Fix | Where | State | Evidence | What it got wrong or blocked along the way |
|---|---|---|---|---|---|
| F1 | Duplicate GFF3 rows from `PASA_transcripts_and_assemblies_to_GFF3.dbi` | PASApipeline `4376a22` (already on `rust_optimize`) | **Done, released in rc.1** | 0 duplicate rows in rc.1 runs | The bug came from the fork's own commit `bce776a`, not upstream. Before the fix it corrupted the TransDecoder genome ORFs; the old A. fumigatus training set was 53% single-CDS against 36%. |
| R1 | `bam2gff3` / `bam2ExonsHints` CIGAR rewrite | funannotate `fix/bam2gff3-cigar` (commit `c134412`, pushed; PR #1210) | **PR open; merge after SELECT commits R3/R5** (D01, D16, D24, D32) | 24 unit tests. Real minimap2 check. A/B test in the rc.1 image: N. crassa valid spliced custom alignments 55 → 9,045, exact RefSeq chains +42%; A. nidulans +1.5% | v0 put the alignment strand on the predict/EVM evidence (caught by Fable review; fixed with `strand=`). v0 used a stricter identity formula (fixed). v1 left insertions right after an intron out of both exons, causing 137/274 PASA "Incontiguous" failures (fixed in v2). The running aligner arms use v1; the final reruns will use v2. |
| R2 | `--UNSPLICED_JOIN_SPLICED`: `?` alignments join the spliced gene covering them | PASApipeline `rust_optimize`, tag **v2.6.1-rc.2** | **Merged and tagged; passed the D39 acceptance rule** (D18, D26) | 14 unit tests; Perl syntax OK in the image | The first end-to-end runs **crashed** in `--ALT_SPLICE` (`orient='?'` hit a CHECK constraint), because mixed-orientation subclusters were not expected downstream. Fixed by storing the assigned orientation. The unit tests had not caught this. Not covered yet: `--gene_overlap` mode and the sharding script. |
| F4 | `--ONE_ALIGNMENT_PER_CDNA`: one alignment per transcript per cluster when several aligners are used | same branch | **Implemented; runs together with R2** (D18) | none yet on its own | The original condition in `Launch_PASA_pipeline.pl:821` (upstream since 2018) never ran with blat + custom alignments. Whether it is "inverted" is a judgment call (Fable argued it is speculation), so it is added as opt-in rather than changing upstream behavior. |
| gmap fallback | Chunked retry, skipping transcripts that crash gmap | same branch | **Implemented and tested; default-on (user-approved exception)** (D22, D27, D28) | 297/298 real contigs aligned; content identical when gmap does not fail | Found because the first gmap arm (29107331) died. The crash reproduces in 3 gmap versions, so pinning a version does not help. |
| gmap root cause | Patch the segfault in gmap itself | `/bigdata/stajichlab/jstajich/projects/funannotate/gmap_debug` | **In progress** (Fable subagent) (D31) | pending | none yet |
| R3 / R5 / F12 | Complete-ORF filter; one model per locus ranked by structure; transitive overlap clustering | funannotate main tree | **Implemented by SELECT; benchmarked** (D08, D09, D17, D23) | Exact-chain precision N. crassa 38.9 → 78.7%, A. nidulans 50.1 → 62.1%, Botrytis 75.5 → 89.4% | SELECT found that `gff2dict` stores `cds_transcript` in genome orientation for minus-strand genes, which would have failed every minus-strand gene in a completeness check. It now checks the protein. |
| R6 | Single-exon genes in Augustus/SNAP training | not started | **Needs a user decision** (D23) | 0% single-exon in every training set against 14-21% in RefSeq | The review first described the single-exon excess as a training-set problem. It is in the PASA models (EVM evidence). The Augustus/SNAP training sets have the opposite problem. |
| Identity gate | Detect divergent RNA-seq (N. crassa reads are from wild isolate HJDF) | proposal (SELECT gate code) | **Needs a user decision** (D19, D21) | N. crassa blat alignments are 74% at 95-98% identity; A. nidulans are 98.8% at ≥99.5% | The mapping-rate gate (95.8%) hid the divergence. The N. crassa benchmark must be read as a divergent-reads case. |
| Aligner choice | gmap versus fixed minimap2 versus blat | experiment `xs/` arms | **A. nidulans done; N. crassa running** (D22) | A. nidulans exact chains: fixed mm2 + blat 3,220; gmap only 3,185; blat only 3,167; fixed mm2 only 3,154 (all ~46% precision) | On clean data, aligner choice barely matters. The divergent N. crassa arms decide the question. |
| Selection × R1 predict arms | Do the training-set gains improve predictions? | SELECT predict arms | **Running** | pending | These arms use the image's **old** `bam2gff3`, so they measure selection changes only. The R1 predict arms will follow. |

## Process problems, and how they were handled

- **Decision-log numbering collisions:** two sessions wrote D16 at the same time, and a "last entry" check gave the wrong number. Rule now: take the highest number (`sort -n`), not the last line.
- **Runs cancelled or rerun:** the first R1 A/B jobs (29107189) were cancelled after the Fable review changed the identity formula. The first R2 jobs (29107369/70) crashed. The first gmap job (29107331) crashed. All were resubmitted with the fixes and recorded.
- **Frozen inputs:** each test run records the md5 of its patched `library.py` (and its extra binds) in `inputs.tsv`. Frozen files are never overwritten while runs use them (`r1_fix_library.py` = v1, `r1_fix_library_v2.py` = v2; `r2_pasa/`, `r3_pasa/`).
- **Test harness files** (`r2_test_train.py`, `r3_test_train.py`) are copies of the image's `train.py` with switches for the experiments. They belong to no branch and must not be merged.
