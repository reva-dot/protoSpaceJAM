# protoSpaceJAM Presentation Notes

## What This Tool Does

`protoSpaceJAM` is a donor-design and guide-selection workflow for CRISPR knock-in experiments.

At a high level, it:

- takes transcript- or coordinate-based user input
- resolves an intended insertion/edit position
- finds candidate gRNAs near that position
- designs HDR donors with homology arms and payloads
- annotates donor features in CSV and GenBank outputs
- performs guide-avoidance recoding to reduce donor recutting

The workflow supports:

- N-terminal tagging
- C-terminal tagging
- direct coordinate-based insertion
- insertion payloads
- SNP-style payload replacement
- dsDNA donors
- ssODN donors
- transcript-aware design from Ensembl

## Core Design Logic

The tool fundamentally separates:

- the intended insertion site
- the gRNA cut site
- the donor sequence

This matters because:

- the donor should be anchored to the insertion site, not the cut site
- guides can be chosen upstream or downstream of the insertion
- `Cut2Ins_dist` is used to rank/filter guides by how close the cut is to the intended insertion

## Major Improvements Made

### 1. CHOPCHOP Guide Coordinate Import Was Fixed

The CHOPCHOP import path had multiple issues that were corrected:

- guide `start/end` are now kept in genomic low-to-high order
- guide strand/orientation is handled separately from genomic coordinates
- PAM interpretation is strand-aware
- cut-site calculation is strand-aware
- reverse-strand guides are now interpreted correctly

This was important because previously:

- some reverse guides had flipped coordinate logic
- PAM placement could be inferred from the wrong side
- cut positions could be off

### 2. Preferred Region Resolution Was Extended and Clarified

The tool now better supports transcript-aware region targeting such as:

- `CDS`
- `CDS_exonN`
- `5UTR`
- `3UTR`
- exon-specific coding regions

Important semantics clarified during development:

- `CDS` refers to the transcript’s translated region, not “everything after a start codon in every exon”
- exon 1 is not guaranteed to contain the CDS start
- some transcripts begin coding in exon 2 or later
- `Preferred_region=CDS` with `Preferred_anchor=start` means the first translated base of the transcript CDS
- `Preferred_region=CDS` with `Preferred_anchor=end` means the last coding base before the stop codon

### 3. Sequence-Window CHOPCHOP Submission Was Added/Fixed

Originally, CHOPCHOP mostly behaved like gene-mode targeting.

The tool was extended so that for region/coordinate-driven workflows it can:

- resolve the desired genomic region
- fetch the actual genomic sequence window from Ensembl
- submit that sequence to CHOPCHOP as `fastaInput`
- recover intronic and nearby noncoding guides from that local genomic window

This was especially important for cases where:

- the relevant coding exon is tiny
- the best cut is in a nearby intron
- CHOPCHOP gene mode does not return the needed intronic guides

Key behavior:

- `--chopchop_target_type WHOLE` is now used as the wrapper-level signal for sequence-window mode
- when a coordinate/resolved site is present, the code forces sequence submission instead of gene-name-only behavior
- the local region can be expanded with `--chopchop_window_padding`

### 4. CHOPCHOP Sequence-Mode Coordinate Remapping Was Fixed

When CHOPCHOP is run on pasted sequence, it often returns local sequence-relative coordinates.

The importer was fixed so it now:

- maps local `seq:<pos>` results back to genomic coordinates
- restores correct `chr/start/end/cut_pos`
- recalculates `Cut2Ins_dist`
- supports downstream region annotation and ranking

This was a crucial fix because without remapping:

- guides looked like `chr=seq`
- cut-region labels were blank or wrong
- distance-to-insert calculations were unusable

### 5. Guide Context Annotation Was Added

The tool now annotates where each guide cuts relative to transcript structure.

Current output fields include cleaned-up labels such as:

- `cut_annotation` in the lean result table
- `cut_region`
- `cut_region_detailed`
- `cut_exon_number`
- `cut_intron_number` in the audit/guide tables

These let users see whether a guide cuts in:

- CDS
- exon
- intron
- 5'UTR
- 3'UTR

And, when applicable:

- which exon
- which intron

This makes it much easier to explain and validate guide choice biologically.

### 6. Guide Ordering Sequence Column Was Added

The original guide sequence fields used DNA alphabet (`T`), which is useful for genomic interpretation but not ideal for ordering.

A separate ordering column was added:

- `guide_ordering_seq`

This provides an RNA-style sequence (`U` instead of `T`) that can be copied more directly for guide ordering.

### 7. Resolved Region Audit Columns Were Added

To make transcript-aware targeting transparent, resolved-region audit columns were added so outputs can show:

- what region was requested
- what exon-specific region was actually used
- the resolved coordinate
- resolved region boundaries

This helps users verify that requests like:

- `CDS,start`
- `CDS,end`
- `CDS, exon 2, start`

resolve the way they intended.

### 8. GenBank Transcript Feature Projection Was Improved

The GenBank output was updated to project actual transcript regions onto donors instead of relying on vague exon-like chunks.

Features now more explicitly reflect:

- exons
- CDS_exonN
- 5UTR_exonN
- 3UTR_exonN

This improved the interpretability of donor maps, especially near exon boundaries and UTR transitions.

### 9. Split-Codon CDS Translation Annotation Was Fixed

One major issue was that CDS translation in GenBank was wrong for exons where codons are split across exon junctions.

Examples like:

- CTNNB1 exon 15
- GZMB exon 3

revealed that the old export logic treated exon-local CDS slices as if they started a fresh ORF.

This was fixed by:

- using the precomputed Ensembl per-position codon phase information
- carrying donor phase tracks through donor feature construction
- building translated CDS features from the coding track plus correct `codon_start`
- avoiding naive local ORF reconstruction where possible

Result:

- GenBank CDS translation now better matches UCSC-style transcript phase behavior
- split codons at exon boundaries are handled much more accurately

## How the Tool Uses Transcript Annotation

The tool uses transcript-aware annotation derived from Ensembl.

The relevant upstream precompute script is:

- [extract_gene_models_info.py](/c:/Users/ok/Desktop/REVA/Masters/Spring%202026/Goodman%20lab/code/protoSpaceJAM/protoSpaceJAM/precompute/scripts/extract_gene_models_info.py)

That script precomputes:

- transcript structure
- exon boundaries
- CDS boundaries
- UTR boundaries
- intron/junction context
- per-genomic-position codon phase for each transcript

The codon-phase information is especially important because it allows correct treatment of cases where:

- the first codon of an exon is completed by bases from the previous exon
- the last codon of an exon is completed by the next exon after splicing

This was the basis for fixing translation annotation behavior in donor GenBanks.

## Important Biology Semantics Clarified During Development

### CDS Start Does Not Have To Be In Exon 1

The tool now explicitly supports the fact that:

- exon 1 may be entirely 5'UTR
- CDS can begin in exon 2 or later

So users should not assume:

- `Preferred_exon=1` means “start of the CDS”

Instead:

- use `Preferred_region=CDS` and `Preferred_anchor=start` for the true transcript CDS start

### Exon 3 Internal Start vs Canonical Start

For CTNNB1, the canonical transcript CDS starts in exon 2.

This means:

- an ATG in exon 3 is an internal ATG, not the canonical transcript start

So if a user targets:

- after the exon 2 start codon
- after the first ATG in exon 3

those are biologically different designs.

### C-Terminus Coordinate Semantics Were Fixed

There was an off-by-one problem in C-terminal resolution because:

- coordinate mode inserts *after* the chosen base
- resolving to the stop-codon base itself places insertion after the stop’s first base

This was corrected so that:

- C-terminal `CDS,end,0` resolves to the base immediately upstream of the stop codon in the way the donor logic expects

## Donor Behavior and Recoding

### What the Donor Workflow Does

For each design, the donor workflow:

- extracts homology arms around the intended insertion site
- inserts the payload between arms
- performs guide-avoidance recoding
- optionally trims dsDNA donors if required by synthesis constraints

### Important Recoding Caveat

The tool can recode donor sequence, but recoding is primarily driven by:

- recut avoidance
- PAM disruption
- CFD reduction
- chimeric cut-site suppression

It is **not** a general synthesis optimizer.

That means:

- problematic sequence far out in a homology arm may remain untouched
- recoding is focused near cut/guide-related regions

### Payload Preservation Is Still Not a True Feature Yet

One important limitation remains:

- there is still no dedicated `preserve_payload` mode

This matters because:

- payload sequence can be biologically sensitive
- even synonymous or local sequence changes may be undesirable

Current workaround:

- users can disable recoding entirely with `--recoding_off`
- or manually reuse homology arms with the original payload sequence

But a real payload-lock feature remains future work.

### Recoding Mutation Reporting Was Fixed

The recoding summary had a bug where case-only differences were being reported as fake mutations, such as:

- `c>C`

This was fixed by:

- ignoring case-only differences
- comparing against the final donor sequence rather than an intermediate donor state

## Homology Arm Features and Synthesis Constraints

### Symmetric Arms

The original donor design mode assumes symmetric dsDNA homology arms, for example:

- `--HA_len 800`

### Asymmetric Arm Support

Support was added for:

- `--left_HA_len`
- `--right_HA_len`

This was intended for cases like:

- LHA can be long and synthesis-friendly
- RHA becomes hard to synthesize beyond a shorter length

However, asymmetric-arm handling exposed several bugs and became an unstable area.

Key asymmetric-arm issues found and fixed:

- wrong slicing using left-arm length when reconstructing the right arm
- broken donor geometry during later recoding phases
- bad GenBank projections when arms were asymmetric

Even with fixes, this path was more fragile than the original symmetric-arm workflow.

### Donor Trimming

The tool supports dsDNA trimming for synthesisability.

Important parameters:

- `--HA_len`
- `--MinArmLenPostTrim`

Behavior:

- start with the requested arm length
- trim if needed
- do not trim below the configured minimum

Example:

- `--HA_len 800 --MinArmLenPostTrim 580`

means:

- start from 800 bp arms
- allow trimming
- do not let arms drop below 580 bp

### Practical Synthesis Strategy Discussed

For synthesis-constrained donors, the recommended practical strategy was:

- use the longest synthesizeable homology arms
- let trimming help where needed
- if total donor length must be increased for vector packaging reasons, consider non-homology filler elsewhere in the construct rather than forcing a bad arm to be longer

## Caching and Performance Improvements

### Ensembl Annotation Cache

The tool already used Ensembl-derived transcript annotation bundles and could cache:

- transcript structure
- region catalogs
- codon-phase maps

### New Persistent Region-Level Sequence Cache

Repeated runs were still slow because genomic sequence windows were being fetched repeatedly from Ensembl REST whenever whole-chromosome pickle files were absent.

This was improved by adding persistent region-level sequence caching under:

- `--ensembl_cache_dir/sequence_cache/...`

What this cache stores:

- genome version
- chromosome
- start/end
- strand
- sequence

What this improves:

- repeated HA extraction
- CHOPCHOP sequence-window mode
- repeated debugging around the same locus

This avoids needing bulky whole-genome/chromosome caches while still speeding repeated design runs.

### Fetch-Path Cleanup

The sequence-fetch logic was also cleaned up so that:

- duplicated CHOPCHOP sequence-fetch logic was consolidated
- cache validation was improved
- returned sequence lengths are checked against the requested genomic interval

## CHOPCHOP Integration Notes

The tool now supports two broad guide retrieval modes:

- gene/transcript-oriented behavior
- local sequence-window behavior

When using local sequence-window behavior:

- the tool resolves a desired region/coordinate
- fetches the local genomic sequence
- submits that sequence to CHOPCHOP
- remaps returned guide positions back to genomic coordinates

This made it possible to recover nearby intronic guides that would otherwise be absent in CHOPCHOP gene mode.

## Output Files and Their Meaning

Important outputs include:

- `guides_from_chopchop.csv`
- `result.csv`
- `result_audit.csv`
- `homology_arms.csv`
- `recoding_mutations.csv`
- donor GenBank files
- payloadless `_gRNAonly_noPayload.gb` files

### `guides_from_chopchop.csv`

Contains:

- guide coordinates
- strand
- cut position
- cut-to-insert distance
- region annotation
- guide ordering sequence

### `result.csv`

Contains the lean user-facing donor/guide summary with:

- insertion position
- guide info
- distances
- combined cut-location label (`cut_annotation`)
- final recut CFD
- donor sequences
- synthesis flags

### `result_audit.csv`

Contains the full audit/detail version of the final design table, including:

- target-region and resolved-region metadata
- intermediate recoding CFD columns
- donor names
- strand/junction diagnostics

### `recoding_mutations.csv`

Summarizes actual sequence changes introduced by donor recoding.

### Donor GenBank Files

Used for visualizing:

- left HA
- payload
- right HA
- transcript-region overlaps
- guide-related features
- recoded regions
- translated CDS features

### Payloadless `_gRNAonly_noPayload.gb`

Useful for visualizing the wild-type arm sequence with the guide context, without the inserted payload complicating the map.

## Important Troubleshooting Lessons Learned

### Guide Distance in Donor Map vs Result Table

The CSV `Cut2Ins_dist` is based on genomic cut vs genomic insert positions.

This can differ visually from donor GenBank maps because:

- the payload splits the donor sequence
- guide remnants or guide-marked regions can appear far from the insert junction in donor coordinates even when the original cut was close

This made some donor `gRNA+PAM` annotations visually misleading.

### Donor Translation vs Exon Slice Translation

A standalone exon-local CDS slice does **not** have to be divisible by 3.

Why:

- codons can be split across exon boundaries

So a GenBank viewer that treats a projected `CDS_exonN` slice as an independent CDS can show misleading frames/stops.

This is why phase-aware translated CDS export was needed.

## Commands Frequently Used

### Example GZMB Annotation Test

```powershell
python -m protoSpaceJAM.protoSpaceJAM --path2csv input\GZMB.csv --annotation_source ensembl --ensembl_cache_dir output\ensembl_cache --num_gRNA_per_design 20 --max_cut2ins_dist 50 --payload_type insertion --payload_file payloads\bfp.txt --outdir output\GZMB --Donor_type dsDNA
```

### Example CTNNB1 CHOPCHOP Sequence-Window Run

```powershell
python -m protoSpaceJAM.protoSpaceJAM --path2csv input\CTNNB1_input.csv --annotation_source ensembl --ensembl_cache_dir output\ensembl_cache --num_gRNA_per_design 20 --max_cut2ins_dist 50 --payload_type insertion --payload_file payloads\bcatenin_mcherry_stop_3utr_polyA.txt --outdir output\bcatenin_mcherry --Donor_type dsDNA --chopchop_target_type WHOLE --clean_genbank_dir --MinArmLenPostTrim 0 --HA_len 450 --chopchop_window_padding 250
```

## Limitations and Deferred Work

### 1. True Payload Preservation Mode

Still not implemented.

Desired future feature:

- preserve payload exactly
- allow HA recoding/trimming around it
- prevent unintended payload sequence drift

### 2. More GenBank Fusion/CDS Polish

Translation annotation is much better now, but still needs broader testing across more genes.

Potential future improvements:

- more explicit fusion CDS labeling
- clearer donor-vs-transcript translation tracks
- better visualization for cut/junction remnants in full donor GenBanks

### 3. More Technical Debt Cleanup in `hdr.py`

`hdr.py` still contains a lot of accumulated workflow logic and historical complexity.

This includes:

- large multi-phase recoding logic
- feature-coordinate derivation
- donor trimming logic
- guide/chimeric-site handling

It works, but it remains a high-complexity file that would benefit from refactoring.

### 4. Broader Validation Across Genes

The split-codon translation fix worked for the tested GZMB case, but it should still be validated on:

- more human genes
- more reverse-strand transcripts
- more mixed UTR/CDS exon cases

## Branch / Recent Development State

Recent development was done on:

- `feat/chopchop-ensembl-ha`

Recent commits included:

- CHOPCHOP import and donor-annotation workflow refinements
- CDS annotation improvements
- sequence caching improvements

## Presentation Framing Suggestions

If presenting this tool, a useful story is:

1. The tool is a transcript-aware CRISPR donor design system.
2. It integrates guide discovery, donor design, recoding, and feature annotation.
3. A major focus of the recent work was making the biology representation more faithful:
   - correct guide orientation
   - correct cut-site placement
   - transcript-aware region targeting
   - intronic guide recovery
   - phase-aware CDS translation at exon boundaries
4. Another focus was improving usability:
   - region audit columns
   - cut-context labels
   - RNA guide ordering sequences
   - persistent Ensembl sequence caching

## Short Summary

In this development cycle, `protoSpaceJAM` was substantially improved in four major ways:

- guide import and CHOPCHOP integration became more accurate and transcript-aware
- donor/genbank annotation became much more biologically faithful
- user-facing CSV outputs became more interpretable
- repeated Ensembl-backed runs became more efficient through persistent sequence caching
