# Code Review: PASApipeline fork (rust_optimize) and funannotate PASA training path

**Date:** 2026-09-25
**Primary reviewer:** Claude Opus 5.5 (`claude-opus-5-5`), running in Claude Code
**Additional reviewers:** 5 Opus 5.5 subagents, 1 Claude Sonnet 5 subagent (`claude-sonnet-5`), 2 Claude Fable 5.1 subagents (`claude-fable-5-1`)
**Repositories and commits reviewed:**

| Repository | Branch | Commit |
|---|---|---|
| PASApipeline (this repo) | `rust_optimize` | `23d67c0` (baseline `master` = `381e96b`) |
| funannotate | `target_1.9/rust_EVM_trinity_PASA` | `41a2fd7` |

No repository code was changed during this review. The only new files are this report and the `code_review_20260925/` folder. That folder holds a proposed replacement converter and its test harness.

---

## 1. The request

The user asked for an expert programmer and bioinformatician review of this PASApipeline fork. The request had these parts:

- The fork replaces slow Perl steps with Rust and optimized C++ code. Check that the logic is still correct.
- Skip all `bench*` folders.
- The mappers (minimap2, blat, gmap) are important. The most important output is high-quality transcripts and spliced gene models for gene-predictor training.
- Consider how funannotate uses PASA (`funannotate/train.py` and `funannotate/update.py` in `~/projects/funannotate/funannotate-live`). These scripts use PASA to update gene models and to build a training set for Augustus and SNAP.
- Test runs show too many low-complexity (single-exon) gene models and poor mapping. This gives poor training sets and weaker evidence for EVM.
- Find hidden logic problems. Suggest improvements to scoring and to heuristics that build spliced gene models and the longest valid ORFs through spliced exons.

A follow-up request asked for three more things:
- Run `git blame` on the funannotate `bam2gff3` parser.
- Have Fable re-review that bug and its fix.
- Have other models debate whether other statistical and gene-prediction strategies are needed.

## 2. How the review was done

| Stage | Agent (model) | Scope | Method |
|---|---|---|---|
| 1 | Assembler reviewer (Opus 5.5) | `pasa_cpp/`, `pasa_rust/pasa-assembler`, `PASA_alignment_assembler.pm`, `assemble_clusters.dbi` | Built upstream (`cb98893`), `master` and HEAD `pasa` binaries and the Rust assembler. Diffed their outputs on 419 inputs (16 test files, 400 random, 3 large). |
| 1 | Alignment reviewer (Opus 5.5) | `Launch_PASA_pipeline.pl`, the aligner scripts, import, validation, clustering; funannotate `bam2gff3` | Read the code paths end to end. Ran `bam2gff3` on synthetic SAM records. |
| 1 | ORF / training reviewer (Opus 5.5) | `pasa_asmbls_to_training_set.dbi`, TransDecoder plugin; funannotate `train.py`, `predict.py`, `library.py` | Read the code. Checked real run outputs in `BENCHMARK/` and `benchmarking/results/`. Ran a synthetic `filterGenemark.pl` test. |
| 1 | Concurrency reviewer (Opus 5.5) | Threaded and sharded `.dbi` scripts, `DB_connect.pm`, slclust (C++ and Rust), cdbtools | Ran `cargo test`. Compared C++ `master`, C++ HEAD and Rust slclust on 3,200 random graphs. Compared cdbyank with cdbyank_rust. |
| 2 | Bug re-review (Fable 5.1) | funannotate `bam2gff3` | Checked the bug adversarially. Ran real minimap2 2.31 on a synthetic genome. Emulated PASA's splice-boundary check. Wrote and tested a fix. |
| 3 | Debate (Sonnet 5) | Statistical and gene-prediction strategy | Criticized the findings brief. Spot-checked 3 findings. |
| 3 | Debate (Fable 5.1) | Cause ranking, ablation design, ORF heuristics | Criticized the findings brief. Spot-checked findings. Measured the reference baseline. |
| All | Primary reviewer (Opus 5.5) | Synthesis | Checked each major claim again in the code or data before including it (see "Verified by" in each finding). |

**Evidence labels in this report:**
- **Executed**: a test or binary was run and the output checked.
- **Data**: checked against real pipeline output files.
- **Code**: the code path was read end to end.
- **Inference**: not tested. Treat these as hypotheses.

---

## 3. Executive summary

The single-exon excess has several causes that add together. The causes are ranked below by the evidence available. This ranking is not yet measured by an ablation (see section 7).

1. **Resolved: the fork's `master` wrote every PASA assembly GFF3 row twice (F1).** TransDecoder's genome mapper then placed CDS coordinates wrongly.
   - In the A. fumigatus runs, the share of single-CDS models in the training set went from 35.9% (v1.8.17) to 53.3% (v1.9.0-beta.5).
   - This is **fixed on HEAD** (`4376a22`, 2026-09-21). **Resolved (user, 2026-09-25):** the duplicated output came from an earlier run. Current runs use the fixed writer. F1 therefore explains only the old beta.5-era numbers. It does not explain the excess that remains today.
2. **Even the v1.8.17 level (35.9%) is 1.85 times the reference.** The Af293 RefSeq annotation has 19.4% single-CDS mRNAs (1,905 of 9,823, measured). The excess that remains comes from upstream PASA design and from funannotate's model selection:
   - Unspliced (`?`-orientation) alignments are kept apart from spliced ones (F3).
   - `getBestModel` ranks models by TPM before structure (F8).
   - TransDecoder keeps several ORFs per assembly, and start refinement changes their labels (F10).
   - No completeness filter is applied before training, and the single-exon filter is all-or-nothing (F7).
3. **funannotate's `bam2gff3` has real coordinate and strand bugs, present since 2018 (F2).**
   - PASA always runs its own blat (or gmap) alignment as well, so a lost minimap2 alignment is often recovered.
   - The most direct effect on single-exon excess is soft-clip inflation. A partial, unspliced alignment passes `MIN_PERCENT_ALIGNED`.
4. **The Rust assembler picks different assemblies from the C++ binary (F5).** The fork's `master` used Rust by default. HEAD prefers C++.
   - The optimized C++ in HEAD gives **byte-identical** output to upstream on all 419 inputs tested.
5. **Augustus accuracy reported by funannotate is inflated (F11).** `etraining` sees the test genes before the split. This does not cause bad models, but it hides them.

---

## 3b. Correction: where the single-exon excess actually lands (SELECT session, DECISIONS D23)
- The single-exon excess measured in this review (for example 55.5% in N. crassa) is in the PASA models: `pasa_assemblies.gff3` and `funannotate_train.pasa.gff3`.
- `selectTrainingModels` then drops all single-exon genes once 200 or more multi-exon models exist. In every benchmarked training set, old code and new, **0% of the Augustus/SNAP training genes are single-exon**. The RefSeq references have 19-21%.
- So the excess affects two things:
  - **EVM evidence:** PASA models get a default weight of 6.
  - **Model choice at each locus:** a single-exon PASA model can replace a spliced one.
- It does not put single-exon genes into Augustus/SNAP training. There the opposite problem exists: the predictors never see single-exon genes. R6, a cap at the BUSCO-estimated share instead of all-or-nothing, addresses this. It is logged as needing a user decision.
- SELECT's new selection (R3/R5) raises exact-match precision against RefSeq. See section 3a.

## 3a. Evidence from current (rc.1, F1-fixed) runs

This section adds data from the Fungi_BFD_runs session (a parallel Claude session working for the same user). The primary reviewer re-measured the items marked "re-measured".

**Runs:**
- Location: `/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/do_pasa_rust_vs_perl/results/<genome>.<rust|perl>/`.
- Images: `funannotate-1.9.0-rc.1.sif` and `rc.1-norust.sif`, rebuilt 2026-09-24 with the fixed GFF3 writer.
- Settings shared by all runs: the same Trinity-GG FASTA, `--pasa_db sqlite`, `--max_intronlen 3000`, `--jaccard_clip`, and **no `--stranded`**.

**F1 is absent.** Each `pasa_assemblies.gff3` has 0 duplicate rows (re-measured for N. crassa and A. nidulans).

**Training sets (`funannotate_train.pasa.gff3`, rust arm):**

| Genome | Reads mapped | Models | Complete ORF (peer) | Single-CDS models (re-measured) |
|---|---|---|---|---|
| *N. crassa* OR74A | 95.8% | 6,817 | 2,146 (31%) | 3,784 (55.5%) |
| *A. nidulans* FGSC A4 | 96.5% | 6,903 | ~5,185 (75%) | 1,664 (24.1%) |
| *M. guilliermondii* | 97.3% | 1,640 | ~1,232 (75%) | 1,475 (89.9%) |
| *C. siamense* CAD1 (reads from another species) | 0.5% | 806 | 159 (20%) | 344 (42.7%) |

- "Complete" means an ATG start, a stop codon, and a CDS length divisible by 3 (peer definition).
- *M. guilliermondii* is a Saccharomycotina yeast. A high single-exon share is expected for this group. No reference was measured.
- No local *N. crassa* reference was found, so its true single-exon share was not measured here. Its reads map well, yet only 31% of models are complete and 55.5% are single-CDS. This makes it the clearest test genome.

**Measured test of F3 (metric M2, re-measured on `pasa_assemblies.gff3`, rust arm).** An assembly "overlaps a spliced assembly" if it overlaps any multi-exon assembly on either strand. Strand was ignored because `?` assemblies have an arbitrary strand.

| Genome | Assemblies | Single-exon | …overlapping a spliced assembly | …inside a spliced assembly's span |
|---|---|---|---|---|
| *N. crassa* | 12,128 | 6,725 (55.5%) | 2,264 (33.7% of single) | 1,014 (15.1%) |
| *A. nidulans* | 9,968 | 2,979 (29.9%) | 776 (26.0%) | 336 (11.3%) |
| *M. guilliermondii* | 2,039 | 1,705 (83.6%) | 298 (17.5%) | 192 (11.3%) |

- **Interpretation:** fixing F3 (R2) can at most absorb the single-exon assemblies that overlap a spliced assembly. In *N. crassa* that is about one third. It is only 15% if absorption requires containment.
- The remaining two thirds (about 4,460 in *N. crassa*) are standalone single-exon loci. Only the training-set filters can deal with them: completeness (R3), intron support (R4), a single-exon cap (R6), minimum length (R11) and homology (R12).
- **Revised priority:** for current runs, R3 and R5 come first (per-model completeness, then structure-first ranking), then R2, then R1. This matches the peer session's suggested order.
- The script for this measurement is `code_review_20260925/m2_single_exon_overlap.py`.

**Where the missing start and stop codons come from (N. crassa rc.1, measured by the primary reviewer):**
- **The genome mapping is faithful.** Of the 6,817 chosen training models, TransDecoder labelled 2,146 `complete`, 1,869 `internal` (no start and no stop), 1,423 `3prime_partial` and 1,379 `5prime_partial`. The peer's sequence-based count, done on the genome, also finds exactly 2,146 complete. So the missing codons already exist in the PASA assembly sequence. They are not introduced when ORFs are mapped back to the genome.
- **Model selection explains only a small part.**
  - Across all assemblies, TransDecoder predicted 4,221 complete ORFs; 2,075 of them were not chosen.
  - 1,746 of those overlap a chosen complete model (redundant).
  - Only 329 overlap a chosen partial model and no chosen complete model. These are the loci where structure-first ranking (R5) can swap a partial model for a complete one. The best case is about 2,475 complete models (36%), up from 31%.
- **Most incomplete models (about 4,300 loci) have no complete ORF in any assembly.** At those loci, the PASA assembly itself is truncated. The possible causes are:
  1. The Trinity contig is a fragment.
  2. The aligners left transcript ends unplaced: soft clips, and short terminal exons that blat, gmap or minimap2 did not align. PASA builds assemblies only from aligned segments, and takes the assembly sequence from the genome, so unaligned ends are lost. `MIN_PERCENT_ALIGNED=90` still allows up to 10% of a transcript to be unaligned.
  3. A spliced alignment failed validation, so only a shorter piece remained.
- **These causes can be told apart** by checking, for each truncated assembly, whether a member transcript has a complete ORF. PASA already runs TransDecoder on the transcripts (`trinity.fasta.clean.transdecoder.gff3`) when `--TRANSDECODER` is set, and `pasa_assemblies_described.txt` links assemblies to transcripts. The comparison scripts delete both files (`run_train_compare.sh:45,53`), so this needs a rerun that keeps them (see R13).

**First R13 result: Cordyceps benchmark, measured with `code_review_20260925/r13_truncation_analysis.py`.**
- **Caveats about this run:** it is the older run. It had the F1 bug, which does not affect the transcript-level ORF calls used here. PASA also ran minimap2 itself, as well as importing it as custom alignments.
- **Across all assemblies:** 15,827 assemblies, of which 11,609 have a complete best ORF. Of the 4,218 incomplete ones:

  | Class | Count | Share of incomplete |
  |---|---|---|
  | A member transcript has a complete ORF (alignment or assembly truncation) | 577 | 13.7% |
  | No member transcript is complete (the transcripts are fragments) | 3,544 | 84.0% |
  | No member transcript has an ORF | 97 | 2.3% |
- **Training models only:** 2,288 are incomplete. 263 of them (11.5%) have a complete member transcript.
- **The truncation class is linked to clipping.**
  - When a complete member exists, the median ORF gain from using it is only 2 aa. So the assembly usually ends just short of the start or stop codon.
  - In that class, 21.5% have a member transcript that is ≥10% soft-clipped in the minimap2 BAM, against 7.3% in the fragment class. The median maximum clip is 3.0% against 0.0%.
- **Conclusion (one run, one genome, old configuration):**
  - Alignment-side trimming explains a minority of incomplete assemblies, roughly 1 in 7.
  - Most of the rest are Trinity fragments, and better alignment cannot recover those. The levers for them are the training filters (R3/R4/R6), or better transcripts: long reads, or re-running Trinity with different settings. The effect of Trinity settings was not tested.
  - The N. crassa rc.1 rerun with `KEEP_R13=1` should confirm this on current code.
- **F4 evidence:** 13,386 of 15,827 assemblies list the same transcript accession twice in `pasa_assemblies_described.txt`. That is, alignments of one transcript from two aligners were assembled together. In this run the two aligners were minimap2 inside PASA and the custom minimap2 import. The current funannotate removes minimap2 from PASA's aligner list, so the blat + custom pair is the remaining case.

**R13b/c/d: can incomplete N. crassa training models be rescued? (rc.1 rust arm, against the RefSeq annotation GCF_000182925.2)**
- **Reference baseline:** 2,312 of 10,812 RefSeq mRNAs (21.4%) have one CDS segment. The training set has 3,784 of 6,817 (55.5%).
- **Complete models:** 1,474 of 2,146 (69%) match a reference model exactly (start, stop and intron chain). 569 have different ends. 30 have the same ends but a different intron chain. 73 have no in-frame reference.
- **Incomplete models are fragments of real genes.** 4,474 of the 4,671 incomplete models (96%) overlap a reference CDS in the same frame. A model covers a median of 44% of its reference CDS.
- **A short genome extension would not rescue them** (`r13_genome_extension.py`, `r13_reference_compare.py`).
  - Only 2-4% of open ends lie within 9 codons of the reference start or stop.
  - The median gap is 146-189 codons when the missing part is in the same exon, and 365-435 codons when it crosses an intron.
  - The blind in-frame walk agrees:
    - ATG within 9 codons upstream in 13-14% of cases, against 9% in the shifted-frame baseline.
    - A stop within 9 codons downstream in only 4-5% of cases, against 26-29% in the baseline.
  - So the reading frame continues well beyond the model ends: the models end in the middle of the coding sequence.
- **Where the missing part is** (`r13_fragment_evidence.py`; 4,474 incomplete models with a reference):

  | Class | Models | Share | What could rescue it |
  |---|---|---|---|
  | Transcript alignments do not cover ≥90% of the gene (evidence gap) | 2,516 | 56.2% | Nothing on the alignment side. Needs more RNA-seq depth, long reads, or protein/BUSCO evidence |
  | Transcript pieces cover the gene, but no single transcript does | 1,613 | 36.1% | Joining the pieces (PASA assembly and clustering) |
  | One transcript covers ≥90% of the gene, but the model is incomplete | 345 | 7.7% | Alignment, validation or assembly fixes |

  - Median share of the reference CDS covered: the model 0.44; all PASA assemblies together 0.58; all minimap2 Trinity alignments together 0.84; the best single transcript 0.49.
  - **Evidence lost inside PASA:** in 1,377 models (30.8%), the Trinity alignments together cover ≥90% of the gene but the PASA assemblies together do not. Candidate causes: validation failures, which the F2 coordinate errors make more likely; clustering; and alignments dropped by `bam2gff3`. The R13 rerun keeps `failed_*_alignments.gff3` so that this can be tested directly.
  - In 616 models (13.8%), the PASA assemblies together cover the gene, but it is split across several assemblies or models.
- **Caveat:** the Trinity alignment GFF3 here was made by the unfixed `bam2gff3`. Its coordinates can shift by a few bases. That does not change coverage at this scale, but the set excludes the alignments the old parser dropped (F2 C3). So 0.84 is a lower bound.

**Important caveat (DECISIONS D19): the N. crassa RNA-seq comes from a divergent wild isolate, not from OR74A.**
- The SRA BioSample for SRR33994706 gives strain "Neurospora crassa HJDF", from red mold tofu in Fujian.
- Its transcripts are 95-98% identical to the OR74A genome (74% of blat alignments), and 19% are 90-95%.
- In A. nidulans, whose reads come from the reference strain, 98.8% of alignments are at least 99.5% identical. There, spliced alignments pass PASA validation at 75% (custom) and 89% (blat), against 0.8% and 34% in N. crassa.
- So the N. crassa numbers below describe the divergent-reads case. That case matters for BFD, because reads are chosen by species taxid. Same-strain data behave much better.
- The RNA-seq gate measures only the mapping rate (95.8% for N. crassa), so it does not detect this divergence.

**R13 rerun on current code (N. crassa rc.1, `KEEP_R13=1`, job 29106979): where the single-exon excess starts.**
- **Validation outcome by aligner** (`alignment.validations.output`, via `r13_lost_evidence.py`):

  | Aligner | Spliced alignments valid | Single-exon alignments valid |
  |---|---|---|
  | custom (minimap2 via unfixed `bam2gff3`) | **55 of 6,973 (0.8%)** | 8,432 of 13,713 |
  | blat | 6,193 of 18,320 (34%) | 5,910 of 10,119 |

  - About 70% of the valid alignments PASA uses are single-exon (14,342 against 6,248 spliced). The bias exists before assembly and before model selection.
  - Custom failures are mostly "Splice site validations failed" (6,725), as expected from F2.
  - blat failures are mostly "Incontiguous alignment + splice site failed" (7,100), so blat on Trinity contigs is also weak here.
- **R13 member test (training models):** 88% of the 4,609 incomplete models come from fragment transcripts. Alignment truncation of transcripts is a minor cause of missing start and stop codons.
- **Lost genes:** 947 genes are covered at least 90% by Trinity alignments but not by PASA assemblies.
  - 65% lose the coverage at validation.
  - 35% have valid alignments covering the gene and lose it in clustering or assembly.
- **Running A/B test of R1** (DECISIONS D12): the same rc.1 image, with only `library.py` replaced by the fixed version (the md5 of the image copy equals the branch base). The runs are on N. crassa and A. nidulans.

**F5 was not tested by the rust/norust comparison (checked by the primary reviewer).**
- The rc.1 image does not set `PASA_ASSEMBLER`, and it contains both `pasa` (C++) and `pasa_rust`. HEAD's `PASA_alignment_assembler.pm:84-86` then picks `pasa` first. So **both arms ran the C++ assembler**.
- The "rust" arm differed only in `slclust_rust`, `cdbyank_rust` and `faidx_rust`. These are preferred by `SingleLinkageClusterer.pm:77` and `CdbTools.pm:81`.
- The `pasa_assemblies.gff3` files differ line by line between the arms. But with IDs removed (chromosome, start, end, strand only), *N. crassa* and *A. nidulans* have **0 differing lines**. The differences are ID numbering only.
- This agrees with the slclust and cdbyank equivalence tests in section 8.

**Production scale (peer data; not re-measured):**
- Scale: 21,762 genomes checked. 8,010 trained Augustus on PASA models. 1,595 of these have fewer than 500 complete PASA models. There are 2,572 rerun candidates (`qc_training_sweep_20260925/`).
- Of the 1,857 flagged PASA-trained genomes, 1,295 belong to species whose reads map at least 50% to the genome. So most of the damage comes from how the training set is built, not from RNA-seq mismatch.
- No production Nextflow module passes `--stranded`. So F3 and F12 apply to all production data.
- Downstream effect: weak Augustus and SNAP models outvote GeneMark in EVM. EVM output fell to 33-66% of merged loci, and BUSCO fell by up to 50 points (for example, *Drepanopeziza* 97.1 → 46.8 and *Ophiostoma* 98.8 → 37).

**Gates already added in funannotate-live (uncommitted). The fixes proposed here must stay compatible with them:**
- `library.count_complete_orf_models()` and `pasa_training_gate()`.
- `predict --min_pasa_complete_models` (default 500) counts complete ORFs before `selectTrainingModels`. Below that number, Augustus and SNAP train from BUSCO.
- `train --min_rnaseq_map_rate` (default 10%) exits before Trinity and PASA when reads map poorly.
- R3 is the per-model version of this gate. The genome-level gate can stay as a backstop.

---

## 4. Findings: PASApipeline fork

### F1. RESOLVED (was MAJOR): Fork `master` printed every assembly GFF3/GTF row twice. This corrupted the TransDecoder genome ORFs. (Fixed on HEAD. The user confirmed the affected output came from an earlier run.)
- **Where:**
  - Writer: `master:scripts/PASA_transcripts_and_assemblies_to_GFF3.dbi` calls `print_GFF_row` at lines 187 and 290. The bug came from commit `bce776a` (2026-06-29, "fix N+1 query").
  - Consumer: `pasa-plugins/transdecoder/util/cdna_alignment_orf_to_genome_orf.pl:175-176` does not remove duplicate exons:
    ```perl
    if (my $struct = $cdna_alignments{$asmbl}) { push (@{$struct->{coords}}, [$lend, $rend]); }
    ```
- **Effect:**
  - Every assembly appears to have its exons twice. The transcript-to-genome coordinate walk then places CDS ends on the wrong copy.
  - Example (asmbl_1098): the transcript CDS is 193-1548 (452 aa, complete). It maps to a 171-nt genomic CDS.
  - Every unspliced assembly appears to have 2 or more exons. Line 275 therefore drops all of their minus-strand ORFs: 4,850 of them.
- **Measured (Data):**

  | Run | Duplicate rows | Single-CDS training models | Mean CDS |
  |---|---|---|---|
  | Cordyceps `BENCHMARK/trinity_rust/work_rust` | 47,990 / 95,980 | n/a | n/a |
  | A. fumigatus v1.8.17 | 0 | 35.9% | 1,118 nt |
  | A. fumigatus v1.9.0-beta.5 | 88,494 | 53.3% | 805 nt |
- **Fix:**
  - The writer is fixed in HEAD `4376a22`.
  - Also add a `%seen` exon filter at line 175 of the consumer, as a guard.
  - Add a check in `runPASAtrain` that fails when `pasa_assemblies.gff3` contains duplicate rows.
- **Verified by:** the primary reviewer. `git show master:` confirms the two calls. `sort | uniq -d` on the Cordyceps file gives 47,990.

### F3. MAJOR (upstream design): single-exon `?` alignments never join spliced loci.
- **Where:**
  - `scripts/assign_clusters_by_stringent_alignment_overlap.dbi:83`, and the same line in `assign_clusters_by_gene_intergene_overlap.dbi:94`:
    ```perl
    $asmbl_id .= ";$spliced_orient";
    ```
  - `scripts/subcluster_builder.dbi:179`:
    ```perl
    unless ($i_spliced_orientation eq $j_spliced_orientation) { next; }
    ```
- **Effect:**
  - A single-exon alignment has spliced orientation `?`. Under `--stringent_alignment_overlap`, which funannotate always passes, it is clustered apart from an overlapping spliced gene. It is then assembled alone and becomes its own subcluster and "gene".
  - Under default clustering (`reassign_clusters_via_valid_align_coords.dbi`), `?` alignments share the cluster. The assembler can then absorb them into a spliced assembly (`PASA_alignment_assembler.pm:187` tolerates `?`). This is a nuance of the Fable debate claim: removing `--stringent_alignment_overlap` helps, but `subcluster_builder` still keeps any `?` assembly that was not absorbed separate.
- **Fix, in order of preference:**
  1. Let a `?` alignment join the `+` or `-` partition of a spliced alignment that contains it, both in clustering and in `subcluster_builder`.
  2. Use `--transcribed_is_aligned_orient` for stranded data. funannotate already does this, except when paired-end and single-end reads are mixed (`train.py:1055` then forces `--stranded no`).
  3. Use `--INVALIDATE_SINGLE_EXON_ESTS`. It keeps full-length (`is_fli`) transcripts, so it fits with `--TRANSDECODER`.
- **Verified by:** the primary reviewer (Code, both locations). This is upstream PASA behavior, not a fork regression. Its share of the excess is Inference.

### F4. MAJOR (upstream, likely inverted condition): multiple alignments of one transcript from different aligners are never reduced to one.
- **Where:** `Launch_PASA_pipeline.pl:821`
  ```perl
  if ($NUM_TOP_ALIGNMENTS > 1 || scalar(@PRIMARY_ALIGNERS) == 0) {
  ```
- **Effect:**
  - `custom` is pushed onto `@PRIMARY_ALIGNERS` at line 496. So the `== 0` branch can only run with no aligners and no custom import (a Cufflinks-only run).
  - With blat + custom minimap2 (the funannotate default is `--aligners minimap2 blat`, `train.py:897`), each transcript can have 2 valid alignments in one cluster. `assemble_clusters.dbi` does not filter per cDNA.
  - Compatible copies merge harmlessly.
  - An unspliced blat copy of a spliced transcript becomes a separate `?` assembly. This interacts with F3.
- **Debate:**
  - Fable points out that the line has been unchanged upstream since 2018. It calls "inverted" speculation until someone measures how often blat and minimap2 disagree.
  - The comment under the condition says that alignments of the same cDNA from different aligners should not both be assembled. That comment supports the "inverted" reading.
- **Fix:** use `scalar(@PRIMARY_ALIGNERS) > 1`. Prefer the spliced alignment over the unspliced one before comparing scores.
- **Verified by:** the primary reviewer (Code). Effect size: not measured.

### F5. MAJOR (fork): the Rust assembler does not match the C++ assembler. `master` used Rust by default.
- **Where:** `pasa_rust/pasa-assembler/src/assembler.rs:157-162`. It emits every trace in a score bin that has any unconsumed alignment. The C++ code (`get_max_missing_Lobj`) repeatedly takes the trace with the most unconsumed alignments.
- **Measured (Executed):**
  - The set of assemblies differs on `rand.txt` (85 vs 87) and on 10 of 400 random inputs.
  - After emulating the Perl selection step, 8 of 400 still differ. Rust keeps extra, partly redundant assemblies and never keeps fewer.
  - Output order is not deterministic: 5 runs gave 5 different md5 sums, because of `HashSet` at lines 458-464.
- **Deployment:**
  - `master`'s `PASA_alignment_assembler.pm` looked for `pasa_rust` first, and `master`'s Makefile always built it.
  - HEAD prefers C++ `pasa`. It uses Rust only when `$PASA_ASSEMBLER` asks for it or `pasa` is missing.
  - The rc.1 production image contains both binaries and does not set `PASA_ASSEMBLER`, so it runs C++. The rc.1 rust/norust comparison therefore did not test F5 (section 3a). F5 matters only for old `master` builds, or if someone sets `PASA_ASSEMBLER=pasa_rust`.
- **Fix:** port `get_max_missing_Lobj` and the C++ tie order to Rust, or remove the fallback. Sort the output before returning it.

### F6. MINOR items (PASA fork)
| ID | Where | Issue | Evidence |
|---|---|---|---|
| F6a | `sharding/shard_subcluster_builder.sh` | If run before `alignment_assembly_loading.ok` exists, it can write an empty output and a `.ok` checkpoint. Launch then skips subclustering without an error. | Inference |
| F6b | `find_alternate_internal_exons.dbi:229-244` | `eval {RunMod}` only prints DB errors. Under concurrent SQLite writes, lost `alt_splice_link` rows go unreported. | Code |
| F6c | `assign_clusters_by_stringent_alignment_overlap.dbi` | With SQLite, each thread holds one write transaction for a whole scaffold (`AutoCommit=0`). This removes the parallelism. A long scaffold can exceed the busy timeout; this failure is loud. | Code + Inference |
| F6d | `pasa_cpp` Makefile / `pasa` | OpenMP is not limited per call, so `-T` workers oversubscribe the cores. Output does not change; performance only. | Executed (output only) |
| F6e | `pasa_asmbls_to_training_set.dbi:178-179` | `-m` is passed to TransDecoder.Predict, which rejects it. It should go to LongOrfs. funannotate does not pass `-m`. | Executed |
| F6f | `PerlLib/fasta.ph` | The new fallback to `fasta36` still passes `-p`. The meaning of `-p` in FASTA36 was not checked. This could change `cDNA_annotation_comparer` identities used by `funannotate update`. | Unverified |
| F6g | `pasa_cpp/tests`, `pasa_rust/*/tests` | No test compares output against the original binary. The Rust tests only check that output is not empty. The divergence in F5 passes all of them. | Code |
| F6h | `import_spliced_alignments.dbi:257, 305, 222` | Three issues: `$matches[0]` is used in a loop; per_id is computed with the wrong formula (later overwritten); and the cluster update is not limited to one aligner (`prog`). No current effect. | Code |
| F6i | `shard_subcluster_builder.sh:44` | `dirname "$0"` breaks when the script is submitted directly with sbatch. This is the same failure class as the BASH_SOURCE rule. | Code |
| F6j | `DB_connect.pm` | Permanent SQLite WAL on GPFS is safe only while all processes run on one node. Document this. | Inference |

---

## 5. Findings: funannotate (consumer of PASA)

### F2. MAJOR: `bam2gff3` gives wrong minimap2 alignment coordinates, Target ranges and strand.
- **Where:** `funannotate/library.py:1838-1950` (the cs-string walk is at 1873-1911).
  ```python
  if x == ":":   matches += n; position += n; querypos += n
  elif x == "-": gaps += 1          # deletion: genome position NOT advanced
  elif x == "+": gaps += 1; querypos += len(...)
  # "*" substitution: no branch -> position and querypos NOT advanced
  query = [1] ... query.append(len(cols[9]))   # soft clips counted as aligned
  ```
- **History (`git blame`):**
  - The walk was written in 2018 by the funannotate author: `a9ed8f2` (2018-02-18), `bc661b4` (2018-03-09), `b77692b` (2018-06-14, "fix bam2gff parsing for crick strand alignments").
  - Commit `ce48682` (2023-04-11) only reformatted the code (quote style and line wrapping). Its diff shows no logic change.
  - Even before June 2018, `*` was counted as a mismatch but never advanced `position`. So the coordinate bug has existed since the code was first written.
- **Fable re-review verdicts (Executed, minimap2 2.31 on a synthetic genome with GT-AG and GC-AG introns, both gene strands, both transcript orientations):**
  - **C1 CONFIRMED.** Every boundary after a substitution or deletion shifts.
    - Example: a transcript with 3 substitutions, a 2-bp insertion and a 2-bp deletion. The output gave exons (201,379)(680,898)(1149,1306)(1707,1945). The truth is (201,380)(681,900)(1151,1310)(1711,1950).
    - An emulated PASA 3-bp splice-boundary check rejects it at 5 boundaries.
  - **C2 CONFIRMED, and worse than first reported.**
    - A trailing clip gives `percent_aligned` = 100 when the true value is 96.4.
    - A leading clip shifts every Target coordinate, so the boundary check fails.
  - **C3 CONFIRMED.** The intron motif is compared with the SAM flag, but the cs motif is always written in reference-forward orientation (`ts:A` is relative to the read).
    - Flag-16 alignments of plus-strand genes and flag-0 alignments of minus-strand genes are rejected. So are GC-AG introns.
    - Only the last intron decides.
    - In the real test, a clean flag-16 transcript without errors was dropped.
  - **C4 PARTIAL.**
    - The Target flip `L - q + 1` for flag 16 is correct for PASA.
    - Column 7 must be the SAM strand, not `ts:A`. PASA works out the spliced orientation itself.
    - `mismatches = NM - gaps` counts indel events where NM counts bases. This has only a small effect at the 80% cutoff.
    - The secondary and supplementary filter and ID uniqueness are correct.
- **Why it went unnoticed (Fable; the first part is verified in code):**
  - funannotate removes minimap2 from PASA's `--ALIGNERS` (`train.py:549`). PASA therefore always runs its own blat or gmap on every transcript and assembles from the combined set. A transcript rejected from the custom set usually still enters through blat or gmap.
  - The loss shows up mainly as less support for long reads with errors and for one orientation class (Inference).
- **Effect on single-exon excess:**
  - The direct route is C2. A partly aligned transcript whose other exons were soft-clipped imports as a valid, "100% aligned" single-exon alignment. blat or gmap would have failed it at 90% aligned.
  - C1 and C3 mostly remove spliced support. Because of blat and gmap, they matter less than first thought.
- **Fix status (2026-09-25):** implemented test-first in the git worktree `~/projects/funannotate/funannotate-live-bam2gff3`, branch `fix/bam2gff3-cigar`. It is **not committed or merged**.
  - One shared parser, `parse_minimap2_splice_record()`, now feeds both `bam2gff3` and `bam2ExonsHints`.
  - For the Augustus hints, the strand of a spliced alignment comes from the intron motif.
  - Tests: `tests/test_bam2gff3.py` has 18 tests. 11 failed on the old code; all 18 pass on the new code. The existing test files also pass.
  - Real minimap2 check: the old code emits 2 of 4 transcripts and shifts the exons of the transcript with errors. The new code emits all 4 with exact exons, and the clipped transcript gets Target 31-830 instead of 1-830.
- **Original proposal:** Fable's drop-in CIGAR-based rewrite is in `code_review_20260925/bam2gff3_fixed.py`. Its test harness is `harness.py`, `make_test.py` and `synth.sam`. The rewrite:
  - Walks the CIGAR string: M/=/X, D, N, I, S and H.
  - Uses the true query length.
  - Flips the Target for flag 16, and uses the SAM strand in column 7.
  - Accepts consensus motifs in either orientation, and requires every intron to be consensus.
  - Computes per-exon identity from cs.
  - Drops flags 0x904.
  - After the fix, all real and synthetic cases give exact coordinates and pass the emulated PASA checks.
  - **Not done:** an end-to-end PASA run with the fixed converter.
- **Also:** `bam2ExonsHints` (about line 1960) has the same walk. Its output feeds Augustus hints. Fable flagged it; it was not reviewed in depth.

### F7. MAJOR: the training-model filters are weak.
- **Where:** `library.py:11096-11106` (`selectTrainingModels`), confirmed by the primary reviewer.
  ```python
  if keeperCheck and k not in keeperList: ...
  if multiCDScheck and len(v["CDS"][0]) < 2: ...
  if not v["protein"][0]: ...
  ```
- **Effect:**
  - When `keeperCheck` falls back to False, 5′-partial, 3′-partial and internal ORFs all pass. In v1.8.17 that is 1,538 of 5,438 models (28%) (Data).
  - Single-exon genes are removed only when at least 200 multi-exon models exist. Otherwise every single-exon gene is admitted.
  - DIAMOND redundancy removal (80% identity and 80% coverage) already exists.
- **Fix:** R1-R5 in section 7.

### F8. MAJOR: `getBestModel` ranks models by TPM first, and its overlap test is one-sided and ignores strand.
- **Where:** `train.py:763-830`. The overlap test is at line 789; the sort key at line 816 is `(TPM, nCDS, cdsLen, overlap)`.
- **Effect:**
  - kallisto TPM is normalized by length. A short 3′ fragment with the same read count as the full-length model gets a higher TPM, so it wins the locus. (Fable debate; this follows from how TPM is defined.)
  - Overlap is measured as the fraction of x that is covered. So a short model with high TPM inside a long model with low TPM keeps both. Measured overlaps between genes in the final sets:

    | Run | Same strand | Opposite strand |
    |---|---|---|
    | v1.8.17 | 433 | 314 |
    | Cordyceps | 283 | 630 |
  - The ignore list for incomplete models (`train.py:667`, `#PROT` lines) is dead code: 0 such lines exist in the checked outputs.
  - `getPASAtranscripts2genes` is never called.
- **Fix:** cluster the models by locus, using a symmetric overlap that respects strand. Pick one model per cluster. Rank by completeness first, then intron support, then CDS length, and use TPM last.

### F9. MAJOR: `filterGenemark.pl` misreads minus-strand genes in funannotate's GTF.
- **Where:** `aux_scripts/filterGenemark.pl:299, 321, 449`.
- **Effect:**
  - The script expects rows in ascending coordinate order. `dict2gtf` writes minus-strand genes in transcription order. As a result:
    - A minus-strand 2-exon gene with a supported intron was marked bad.
    - A complete minus-strand single-exon gene was also marked bad.
  - Line 449 passes every complete single-exon gene without any evidence.
  - The keeper set therefore leans toward plus-strand genes and single-exon genes (Executed, synthetic GTF).
- **Fix:** sort CDS rows by coordinate before this step. Or replace the script with a Python check that every intron matches `hints.ALL.gff` exactly (R3).

### F10. MINOR: TransDecoder settings allow several ORFs per assembly and label refined starts "complete".
- **Where:** `pasa_asmbls_to_training_set.dbi:182-186` does not pass `--single_best_only`. The other settings are in the TransDecoder plugin files:
  - `select_best_ORFs_per_transcript.pl:195` removes ORFs only when they overlap by 10%.
  - `start_codon_refinement.pl:147-148` changes `5prime_partial` to `complete` after moving the start up to 15% downstream.
- **Data (v1.8.17):**
  - 5,548 of 11,904 assemblies have more than one ORF.
  - 6,653 ORFs lie on the minus strand of their assembly.
  - 521 starts were revised.
- **Effect:** short ORFs confined to one exon of a spliced assembly can compete at the locus (see F8). "Complete" models with an N-terminal truncation teach Augustus the wrong start context (Inference).
- **Fix:** add an option for `--single_best_only`. Use `--no_refine_starts` for training sets, or drop models whose start was revised.

### F11. MINOR (hides the problem): Augustus test-set leakage.
- **Where:** `library.py:11505-11517` runs `etraining species TrainSet` on the full set before `randomSplit`. It then reports accuracy on `TrainSet.test`. Line 11574 repeats this pattern.
- **Effect:** the reported sensitivity and specificity are inflated. Poor training sets therefore look acceptable.
- **Fix:** run `etraining` on `TrainSet.train`.
- **Verified by:** the primary reviewer, the Sonnet debate agent and the Fable debate agent.

### F12. MINOR: the selection overlap removal is not transitive, and one training fallback discards strandedness.
- **Non-transitive overlap removal:** `library.py:11185`. If A overlaps both B and C, but B and C do not overlap, A and B can both pass. Fix: build overlap clusters.
- **Strandedness reset:** `train.py:1055` sets `--stranded no` when paired-end and single-end reads are combined. PASA then loses `--transcribed_is_aligned_orient`, so single-exon alignments stay `?` (this makes F3 worse).

---

## 6. Multi-model debate

The two debate agents received the same findings brief (`code_review_20260925/findings_brief.md`). Each was told to argue against it where the evidence justified that.

### 6.1 Points of agreement (all three models)
- The duplicate-row bug (F1), the `bam2gff3` bugs (F2), `?` isolation (F3), weak training filters (F7), TPM-first ranking (F8) and etraining leakage (F11) are real. Each was spot-checked in the code by at least two models.
- Fix the upstream root causes before adding any cap on the single-exon fraction.

### 6.2 Sonnet 5 position
- **Estimate the single-exon target from BUSCO, not from PASA.**
  - Use the exon-count distribution of BUSCO single-copy complete genes in the same assembly. funannotate already runs BUSCO.
  - Using high-confidence PASA models to set the target is circular.
  - Intron density varies 3-10 times across fungi. The example given was Saccharomycetales at about 5% intron-containing genes versus Pezizomycotina at 40-60% or more. Sonnet gave these numbers without a citation.
- **Down-sampling changes Augustus's intron model.**
  - A cap applied only by down-sampling changes the intron submodel. Use the cap as a QC gate, not as the main correction.
  - This is Sonnet's inference.
- **Split by locus and report confidence intervals.**
  - `randomSplit.pl` splits by gene record. Paralogs and near-identical loci can fall on both sides of the split. Split by locus or cluster instead.
  - Report binomial confidence intervals. With n=100 test genes, an observed 80% sensitivity has a 95% CI of about 71-87%.
- **Use other evidence-selection strategies.**
  - BRAKER3/TSEBRA-style selection by intron and CDS support fraction.
  - GeneMark-ETP.
  - Training sets from miniprot protein hints.
  - funannotate's existing BUSCO training path, which is a lower-variance entry point.
- **Fit EVM weights on held-out data.** Fit the weights from per-source accuracy on held-out genes. Do not use a fixed PASA weight of 6.
- **Filter training ORFs.**
  - Remove ORFs that overlap RepeatMasker repeats (TE-derived ORFs).
  - Check codon usage and GC content.
- **Tune minimap2 for fungi.**
  - Set `-G` to a fungal intron length (for example 2-5 kb) instead of the 200 kb default.
  - Consider `--splice-flank=no`.
  - Compare minimap2 alone with gmap alone on a well-annotated chromosome.
- **Other single-exon sources:**
  - Intron-retention isoforms.
  - Antisense transcripts and genomic DNA contamination.
  - Read-through or fused transcripts in gene-dense genomes.
- **Data gaps Sonnet stated:** test-set n, paralog leakage rate, TE-overlap fraction.

### 6.3 Fable 5.1 position
- **Baseline measured:** Af293 RefSeq has 19.4% single-CDS mRNAs. The primary reviewer re-measured this and got 1,905 of 9,823.
  - F1 explains about half of the beta.5 excess.
  - A 1.85-fold excess remains even at v1.8.17.
- **Disagreed with F3's fix:** removing `--stringent_alignment_overlap` alone will not fix `?` isolation, because `subcluster_builder.dbi` also requires the same orientation. The primary reviewer accepts this in part (see the F3 nuance).
- **Disagreed with F4's framing:** "likely inverted" is speculation. Measure how often blat and minimap2 disagree first.
- **Rated F11 (leakage) as not a cause, but still needed** for an honest ablation metric.
- **Cause ranking (Fable's estimate, not measured):** F1 > F2 > F3 > F8 > F10 > F4/F5 > others.
  - F1 has since been confirmed resolved in current runs. For current runs the ranking therefore starts at F2/F3. For mostly unstranded data, the primary reviewer puts F3 first (see section 9).
  - The Fable bug re-review later lowered the weight of F2, because blat and gmap recover the lost alignments.
- **ORF heuristics:**
  - Established: `--single_best_only`; complete ORFs only; a minimum protein length; `--retain_blastp_hits` / `--retain_pfam_hits`.
  - Use the NMD rule (a stop more than 50-55 nt upstream of the last exon junction) only as a penalty. It is weaker in fungi than in mammals.
  - Suggested, untested:
    - For a spliced assembly, require the ORF to span at least one intron.
    - Penalize an ORF confined to one exon.
    - Use this score: CDS length × (1 + introns spanned) × fraction of the transcript covered by the ORF, minus a penalty for introns in the 3′UTR.
    - For a single-exon assembly, require the ORF to cover at least 50% of the transcript.
    - Keep one ORF per assembly.

### 6.4 Adjudication by the primary reviewer
- Both debate agents are right that no fix should be accepted on argument alone. The ablation in 7.2 is needed.
- The BUSCO-based estimate of the single-exon target (Sonnet) is the better estimator. It is specific to the genome and independent of PASA. Use the Af293 reference only to validate the method.
- A locus-level train/test split and binomial CIs (Sonnet) should be part of the fix for F11.
- A higher EVM weight for PASA should not be trusted until the training set is fixed. Fitting weights on held-out data is established EVM practice (Haas et al. 2008, *Genome Biology*). The fitting was not tested here.
- The claim of "roughly half of spliced alignments lost" for F2 is an **expected value from the orientation rule. It was not measured on real data.** Measure it with one pass over the real BAM (see 7.2, metric M1).

---

## 7. Recommendations

### 7.1 Code changes, in priority order
| # | Change | Where | Basis |
|---|---|---|---|
| R0 | Optional guard only: F1 is fixed and current runs are clean. Add a `%seen` exon filter in the ORF mapper, or a duplicate-row check, so the bug cannot return unnoticed. | `cdna_alignment_orf_to_genome_orf.pl:175`; `train.py` `runPASAtrain` | F1, resolved |
| R1 | Replace `bam2gff3` with the CIGAR-based version. Also review `bam2ExonsHints`. | `library.py:1838`, `code_review_20260925/bam2gff3_fixed.py` | F2, executed |
| R2 | Let `?` alignments join a spliced partition that contains them, in clustering and in subcluster building. | `assign_clusters_by_*_overlap.dbi:83/94`, `subcluster_builder.dbi:179` | F3, code |
| R3 | Train only on complete ORFs. Drop models with an in-frame stop and models whose start was revised (or use `--no_refine_starts`). | `selectTrainingModels` (~11102); `pasa_asmbls_to_training_set.dbi` | Established: Augustus and BRAKER train on complete genes |
| R4 | Require every intron to be supported by RNA-seq junctions. Do this in Python and retire `filterGenemark.pl`. | `predict.py` / `library.py` | Established: BRAKER hint support. F9, executed |
| R5 | One model per locus, using strand-aware symmetric overlap clusters. Rank by completeness, then intron support, then CDS length, then TPM. | `getBestModel`, `selectTrainingModels` | Established: Augustus and SNAP need non-overlapping genes. F8, F12 |
| R6 | Replace the all-or-nothing single-exon filter with a cap at the BUSCO-estimated single-exon fraction. Use it as a QC gate. | `selectTrainingModels` | Suggestion (Sonnet); threshold must come from data |
| R7 | Fix the etraining leakage. Split by locus, and report binomial CIs. | `trainAugustus` (~11505, ~11574) | F11; suggestion (Sonnet) |
| R8 | Add a `--single_best_only` option. Consider the ORF-spans-intron rule for spliced assemblies. | `pasa_asmbls_to_training_set.dbi` | Established option; the rule is a suggestion (Fable) |
| R9 | Change `Launch_PASA_pipeline.pl:821` to `@PRIMARY_ALIGNERS > 1`, after measuring aligner disagreement. | Launch | F4 |
| R10 | Make Rust assembler selection match C++, or remove the fallback. Add golden-output tests against upstream `cb98893`. | `pasa_rust`, tests | F5, F6g |
| R11 | Minimum CDS length for training (for example 300 aa, as in PASA's own `extract_reference_orfs.pl`). | `selectTrainingModels` | Suggestion. No data supports a specific value. |
| R13 | **Implemented as opt-in `KEEP_R13=1` in `run_train_compare.sh`, plus `code_review_20260925/r13_truncation_analysis.py`. First result in section 3a.** Find why assemblies are truncated. Keep `trinity.fasta.clean.transdecoder.gff3`, `pasa_assemblies_described.txt`, `valid_*`/`failed_*_alignments.gff3` and `alignment.validations.output`. For each partial assembly ORF, test whether a member transcript has a complete ORF (yes = alignment or assembly truncation; no = the Trinity contig is a fragment). Also measure soft-clip length per alignment in the minimap2 BAM. | PASA work dir; `run_train_compare.sh` | Measured need: about 4,300 N. crassa loci have no complete ORF in any assembly |
| R12 | Homology or BUSCO support with start and stop agreement. Remove ORFs that overlap repeats. | `selectTrainingModels` or TransDecoder `--retain_blastp_hits` | Suggestion (Sonnet and Fable) |

### 7.2 Ablation plan to measure what each fix contributes
This combines the Fable design with Sonnet's additions.

- **Genomes:**
  - A. fumigatus Af293: RefSeq reference, 19.4% single-CDS.
  - One intron-rich genome with a curated reference, for example *N. crassa* OR74A.
  - Use the same Trinity assembly and BAM in every arm.
- **Arms:** each arm adds one change to the previous arm. Also run each change alone against A0 to detect interactions.
  - A0 = HEAD as-is
  - A1 = R1 (bam2gff3)
  - A2 = R2 (`?` joins spliced loci)
  - A3 = R5 (structure-first ranking)
  - A4 = R3 + R8 (complete ORFs, single best ORF, no start refinement)
  - A5 = R9
  - A6 = C++ vs Rust assembler
- **Metrics:**
  - M1: at the alignment stage, the share of spliced alignments kept, and the count of flag/motif disagreements in the BAM.
  - M2: at the assembly stage, the share of single-exon assemblies, and the share of those contained in a spliced assembly on the same strand.
  - M3: in the training set, the single-CDS share against the BUSCO estimate and against the reference; the median CDS length; and exact intron-chain matches to the reference (gffcompare `=`).
  - M4: for predictions, measure these against the reference:
    - gffcompare exon-level and intron-chain sensitivity and precision
    - false single-exon loci
    - BUSCO completeness
    - Augustus held-out accuracy after the R7 fix
- **Statistics:**
  - Run 5 seeds for `randomSplit` and for the Rust arm.
  - Use 1,000 bootstrap resamples of reference genes for 95% CIs on the gffcompare metrics.
- **Meaningful effect (Fable's proposal; the thresholds are judgment calls, not derived from data):**
  - The training-set single-CDS share falls within 5 percentage points of the target.
  - Intron-chain sensitivity improves by at least 3 points, with CIs that do not overlap.
  - BUSCO completeness does not drop.

---

## 8. Checked and found correct
- **C++ assembler:** HEAD, `master` and upstream `cb98893` gave byte-identical output on 419 inputs, with 1, 8 and 16 threads. `canMerge`, `mergeAlignments`, the sparse compatibility rows, the bitset Lobject, the early exit in the overlap scan and the sort-comparator fix keep behavior the same (Executed).
- **slclust:** C++ `master`, C++ HEAD and Rust gave identical partitions on 3,200 random graphs and on one graph with 300k nodes (Executed).
- **cdbyank_rust vs cdbyank:** identical output, except that Rust drops a trailing blank line (Executed). All Rust unit tests pass.
- **Threaded clustering:** partitions are disjoint. Each thread has its own DB connection. Thread failures make the step fail. Cluster membership matches `master`; only the ID numbering changes (Code).
- **Alignment validation:** `validate_alignments_in_db.dbi` applies the same thresholds to all aligners, including custom. The splice-boundary window and the consensus offsets are correct on both strands (Code).
- **HEAD `Ath1_cdnas.pm` fix (`3a170e7`):** restores upstream orientation handling. `master` forced the aligned orientation over the one computed from splice sites. Side effect: in stringent mode, HEAD puts all single-exon alignments in `?` partitions, so HEAD may show *more* single-exon assemblies than `master` until R2 is done (Code + Inference).
- **Fork fixes that change output on purpose:** the GFF3 duplicate fix (F1); `classify_alt_splice_isoforms.dbi` now keeps all subclusters instead of only the last; the removal of dead duplicate subs in `find_alternate_internal_exons.dbi`; the failed-alignment BED step, which never ran before, now runs.
- **The outputs funannotate needs** (`pasa_assemblies.gff3`, `assemblies.fasta`, `pasa_assemblies_described.txt`, `valid_*_alignments.gff3`) are still produced.

## 9. Open questions and data gaps
1. **Answered (user, 2026-09-25):** the duplicated rows came from an earlier run, and the duplication is fixed. Current single-exon excess must come from F2-F12, not F1. The A. fumigatus v1.8.17 set (35.9% single-CDS against the 19.4% reference, with no duplicate rows) is the relevant baseline for the current code.
2. **Answered (user, 2026-09-25): most RNA-seq libraries are unstranded; some are stranded.** For the unstranded majority, this has four effects:
   - **F3 is at full strength.** Without `--transcribed_is_aligned_orient`, no alignment gets an orientation from the library. Every single-exon alignment stays `?`, and in stringent mode none of them can join a spliced locus. So R2 (let `?` join a spliced partition that contains it) is the most important structural fix for these datasets. The strandedness flag cannot replace it.
   - **F2 (C3) applies to about half of the spliced Trinity contigs.** An unstranded Trinity contig can come out in either orientation. The flag/motif rule rejects spliced alignments in one orientation class. This 50% figure is expected from the rule; it is not measured. M1 would measure it.
   - **Unspliced assemblies have no strand information**, except from ORF direction. TransDecoder scans both strands of them. An antisense or non-coding single-exon assembly can therefore still produce a "complete" ORF. This is one more reason for R3, R4, R11 and R12 (complete ORFs, intron support, minimum length, homology), and for limiting single-exon models to the BUSCO-estimated fraction (R6).
   - **For the stranded minority**, pass `--stranded` so that PASA gets `--transcribed_is_aligned_orient`. Do not combine paired-end and single-end reads, because `train.py:1055` then resets strandedness to `no`.
3. On the real BAM, how many spliced minimap2 alignments show a flag/motif disagreement (M1)? This sets the real size of F2.
4. How often do blat and minimap2 alignments of the same transcript disagree on the intron chain (F4)?
5. What is the BUSCO-estimated single-exon fraction for each target genome (R6)?
6. What test-set size does `trainAugustus` use in practice? The CI width depends on it.
7. What does `-p` mean in FASTA36 (F6f)?

## 10. Files
- `CODE_REVIEW_20260925.md`: this report.
- `code_review_20260925/findings_brief.md`: the brief given to the debate agents.
- `code_review_20260925/bam2gff3_fixed.py`: the proposed CIGAR-based converter, by Fable 5.1. Not yet run through a full PASA pipeline.
- `code_review_20260925/harness.py`, `make_test.py`, `synth.sam`, `truth.json`: the test harness. It needs minimap2 and samtools, for example from `funannotate-live/.pixi/envs/default/bin`, and is run with `/usr/bin/python3.12`.

## 11. Decision log and collaboration
From 2026-09-25, all decisions and code changes are recorded in a shared, append-only log: `/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/do_pasa_rust_vs_perl/DECISIONS.md`.
- For each decision, the log gives the owner session, the status (proposed / accepted / done), the measured evidence with file paths, the alternatives considered, and the exact code location (repo, worktree, branch, function).
- The log also lists the worktrees, file ownership and merge order that the two Claude sessions (REVIEW and SELECT) use so they do not overwrite each other.

## 12. Fix status and changelogs
- **Overall status** of every fix: what is done, tested, pending, and what each got wrong along the way. See `code_review_20260925/FIX_STATUS.md`.
- **Changelogs:**
  - `../PASApipeline-r2/Changelog.txt`, "Unreleased" section: R2, F4 and the gmap crash fallback.
  - `~/projects/funannotate/funannotate-live-bam2gff3/CHANGELOG.md`, "Unreleased → Fixed": the R1 `bam2gff3` rewrite.

## 13. Assessment record (2026-09-26)
The assessment of PASA-trained vs BUSCO-trained predictors, the value of RNA-seq evidence, and the performance of the fixes (accuracy and runtime) is recorded in `docs/assessment_2026-09/README.md`. The folder also holds the result tables it cites (`data/`), snapshots of both methods drafts, and the decision log.
