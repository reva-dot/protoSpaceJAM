![Python version](https://img.shields.io/badge/python-3.9%20%7C%203.10-blue)
![license](https://img.shields.io/badge/license-BSD--3-brightgreen)

# protoSpaceJAM

`protoSpaceJAM` designs CRISPR knock-in donors from transcript- or coordinate-level inputs.

The current default workflow is:
- use CHOPCHOP to generate guides
- use Ensembl REST for transcript/exon/CDS/UTR annotation
- use protoSpaceJAM to build donors, recode them, and export design files

You do **not** need to download the old precomputed genome-wide guide database to use this workflow.

## What A New User Needs To Know

The simplest supported path today is:
1. install the package in a Python/Conda environment
2. provide either a CSV of targets or one direct CLI target
3. let protoSpaceJAM fetch transcript annotations from Ensembl and guides from CHOPCHOP
4. inspect `result.csv`, `guides_from_chopchop.csv`, `homology_arms.csv`, GenBank files, and recoding summaries

This means:
- internet access is required for CHOPCHOP and Ensembl-based runs
- precomputed local guide resources are optional, not required
- the default guide ranking preserves CHOPCHOP order within your cut-distance filter

## Installation

### 1. Clone the repository

```sh
git clone https://github.com/czbiohub/protoSpaceJAM
cd protoSpaceJAM
```

### 2. Create and activate an environment

```sh
conda create -y -n protospacejam python=3.9
conda activate protospacejam
```

### 3. Install the package

```sh
pip install .
```

### 4. Sanity-check the CLI

```sh
python -m protoSpaceJAM.protoSpaceJAM --help
```

If that prints the CLI help text, the installation is working.

## First Run

### Option A: CSV input for one or many targets

A fill-in template is included at:
- `input/multi_target_template.csv`

Run it like this:

```sh
python -m protoSpaceJAM.protoSpaceJAM --path2csv input/multi_target_template.csv --outdir output/example_run --annotation_source ensembl --ensembl_cache_dir output/ensembl_cache --num_gRNA_per_design 20 --max_cut2ins_dist 50 --payload_type insertion --payload_file payloads/bfp.txt --Donor_type dsDNA
```

### Option B: One target directly from the command line

```sh
python -m protoSpaceJAM.protoSpaceJAM --input_mode direct --direct_gene_name PDCD1 --direct_ensembl_id ENST00000334409.10 --direct_target_terminus N --outdir output/pdcd1_example --annotation_source ensembl --ensembl_cache_dir output/ensembl_cache --num_gRNA_per_design 20 --max_cut2ins_dist 50 --payload_type insertion --payload_file payloads/bfp.txt --Donor_type dsDNA
```

## Defaults That Matter

You usually do **not** need to pass these explicitly anymore:
- `--guide_source chopchop`
- `--specificity_backend default`

Those are already the defaults.

Current default behavior:
- guide source: `chopchop`
- annotation source: `auto` unless you set `--annotation_source ensembl`
- specificity backend: `default`
- CHOPCHOP ordering is preserved after filtering
- `max_cut2ins_dist` is used as a cutoff, not a reranking score

## Input Modes

### CSV mode

Use:
- `--path2csv <file>`

This is the recommended mode for:
- multiple genes
- multiple regions per gene
- reproducible batch runs

### Direct mode

Use:
- `--input_mode direct`

This is useful for:
- one-off tests
- debugging a single transcript
- quickly checking one N-term or C-term design

## CSV Input Template

A fill-in template is available at:
- `input/multi_target_template.csv`

The input CSV supports three targeting modes per row:
- transcript terminus targeting with `Target_terminus`
- explicit coordinate targeting with `Chromosome` and `Coordinate`
- preferred-region targeting with `Preferred_region` and related columns

### Required column

- `Ensembl_ID`: transcript ID used for annotation, guide retrieval, and donor design

### Optional columns

- `Entry`: custom row label written into outputs
- `Gene_Name`: gene symbol used for CHOPCHOP web gene-based queries
- `Target_terminus`: `N`, `C`, or `ALL`
- `Chromosome`: chromosome name for explicit coordinate targeting
- `Coordinate`: genomic coordinate for explicit coordinate targeting
- `Preferred_region`: one of `exon`, `CDS`, `5UTR`, `3UTR`, or `transcript`
- `Preferred_exon`: exon number for region requests like `exon` or `CDS` within a specific exon
- `Preferred_anchor`: `start`, `center`, or `end`
- `Preferred_offset`: integer offset in transcript orientation after the anchor is chosen

### Precedence rules

Per row, protoSpaceJAM resolves targeting in this order:
1. `Chromosome` + `Coordinate`
2. preferred-region targeting
3. `Target_terminus`

So:
- if `Chromosome` and `Coordinate` are filled, they take precedence
- if you use preferred-region targeting, `Target_terminus` becomes mostly descriptive
- do not mix explicit coordinates with preferred-region fields on the same row

### CDS and exon targeting

This is the most important rule for knock-in rescue designs:
- `Preferred_region=CDS` means coding sequence only
- `Preferred_region=CDS` + `Preferred_exon=1` means the CDS portion of exon 1, not the whole exon

This is useful when you want guides whose **cut position** lands in the coding segment of exon 1 rather than in UTR sequence.

### Default interpretation

Standard terminus targeting is CDS-anchored by default:
- `Target_terminus=N` targets the CDS start / start-codon boundary
- `Target_terminus=C` targets the CDS end / stop-codon boundary

### Example rows

```csv
Entry,Ensembl_ID,Gene_Name,Target_terminus,Chromosome,Coordinate,Preferred_region,Preferred_exon,Preferred_anchor,Preferred_offset
1,ENST00000334409.10,PDCD1,ALL,,,,,,
2,ENST00000334409.10,PDCD1,,,,CDS,,start,0
3,ENST00000334409.10,PDCD1,,,,CDS,1,start,0
4,ENST00000334409.10,PDCD1,,,,exon,1,center,0
5,ENST00000614167.2,B2M,C,,,,,,
6,ENST00000611116.1,TRAC,,,,3UTR,,start,0
7,ENST00000334409.10,PDCD1,,chr2,241849884,,,,
```

## Ensembl Annotation And Region Resolution

When Ensembl-based annotation is used, protoSpaceJAM fetches transcript structure from Ensembl REST and builds a normalized region map locally.

The raw annotation bundle is cached as:
- `GRCh38_<ENST>.annotation.json`

The normalized region catalog is cached as:
- `GRCh38_<ENST>.regions.json`

The region file includes:
- transcript span
- full `CDS`
- full `5UTR`
- full `3UTR`
- transcript-ordered exons: `exon1`, `exon2`, ...
- CDS-within-exon regions: `CDS_exon1`, `CDS_exon2`, ...

This file is the easiest way to verify what protoSpaceJAM actually used for a request like:
- `Preferred_region=CDS` + `Preferred_exon=1`

### Cache reuse

To reuse previously fetched Ensembl annotations, keep using the same cache directory:

```sh
--annotation_source ensembl --ensembl_cache_dir output/ensembl_cache
```

There is no separate "use cache" flag.
If the relevant files already exist in `--ensembl_cache_dir`, protoSpaceJAM reuses them.

## Outputs

Each run writes a main output folder containing the design tables and per-design GenBank files.

### Main tables

- `result.csv`
  - one row per final donor design
  - includes guide coordinates, cut/edit distance, recut CFD summaries, donor names, and final donor sequence

- `guides_from_chopchop.csv`
  - guide candidates that survived filtering
  - includes:
    - `chopchop_rank`
    - `cut_pos`
    - `target_region_label`
    - `target_region_start`
    - `target_region_end`
    - `MM0`, `MM1`, `MM2`, `MM3`
    - `Cut2Ins_dist`

- `homology_arms.csv`
  - donor broken into major sequence pieces
  - includes:
    - `left_HA`
    - `payload`
    - `right_HA`
    - `donor_final`

- `recoding_mutations.csv`
  - recoding audit table
  - includes:
    - `donor_before_recoding`
    - `donor_after_recoding`
    - `donor_final`
    - `mutation_count`
    - `recoding_mutations`

### GenBank files

GenBank files are written to:
- `genbank_files/`

These contain the final donor sequence and donor-feature annotations.

## CFD / Recut Interpretation

protoSpaceJAM computes several donor recut metrics:
- `cfd_before_recoding`
- `cfd_after_recoding`
- `cfd_after_windowScan_and_recoding`
- `max_recut_cfd`

In plain language:
- `before` = predicted recut risk before donor recoding
- `after_recoding` = predicted recut risk after the main recoding step
- `after_windowScan` = worst remaining high-risk local window after scan/refinement
- `max_recut_cfd` = final recut-risk summary reported for the donor

## Existing Donor CFD Utility

A standalone utility is included for scoring an existing donor or homology arm sequence against a guide:
- `protoSpaceJAM/util/score_existing_donor_cfd.py`

Example:

```sh
python protoSpaceJAM/util/score_existing_donor_cfd.py --guide-seq GAGTCTCTCCTCTTCTTTGA --sequence-file my_donor.fa --top-n 25 --output-csv output/existing_donor_cfd_hits.csv
```

This scans the supplied sequence on both strands and reports the highest-scoring CFD windows.

## Advanced / Legacy Local Workflows

The old precomputed genome-wide guide-resource workflow is still part of the repository, but it is no longer required for standard CHOPCHOP + Ensembl usage.

If you specifically want to use precomputed local guide resources, follow the legacy/precompute tooling in:
- `protoSpaceJAM/precompute`

For most new users, you can ignore that directory.

## Troubleshooting

### My guides are not in the region I expected
Check:
- `guides_from_chopchop.csv`
- `target_region_label`
- `target_region_start`
- `target_region_end`
- `cut_pos`
- `GRCh38_<ENST>.regions.json`

### My row is not using the preferred region I asked for
Check whether the CSV row also has:
- `Chromosome`
- `Coordinate`

Explicit coordinates take precedence over region targeting.

### I want the start of the 3' UTR
Use:
- `Preferred_region=3UTR`
- `Preferred_anchor=start`

Do not rely on `Target_terminus=N/C` to steer a preferred-region request.

### The Ensembl fetch step is slow or flaky
Use a persistent cache directory:
- `--ensembl_cache_dir output/ensembl_cache`

That lets future runs reuse transcript annotation without re-fetching it.

## License

Distributed under the terms of the BSD-3 license, `protoSpaceJAM` is free and open source software.
