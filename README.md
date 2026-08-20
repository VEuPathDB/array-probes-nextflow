# array-probes-nextflow

A Nextflow pipeline that maps microarray probe sequences onto a genome to build (or refresh) gene-to-probe mappings and platform definition files.

## Overview

Microarray platforms are defined by probe sequences and the genes (or genomic loci) those probes represent. As genome assemblies and annotations are updated, the gene targets of existing probes can change, so probe-to-gene mappings must be periodically regenerated. This pipeline aligns a platform's probe FASTA against a current genome assembly and annotation, derives an up-to-date probe-to-gene mapping, and — depending on platform type — produces vendor-specific platform files (Affymetrix CDF or NimbleGen NDF). It is used within VEuPathDB's ReFlow-orchestrated genomic data workflows to keep microarray platform definitions current with the latest genome builds.

The pipeline supports two probe types, controlled by `params.platformType`:

- **`expression`** — spliced (mRNA-based) probes. Probes are mapped to the genome with GSNAP (spliced-aware alignment against a GMAP index built from the genome and, optionally, known splice sites from a GTF), intersected with gene annotation to build a gene-to-probe mapping, and used to regenerate vendor platform files.
- Any other value (e.g. genomic/unspliced probes) — probes are mapped to the genome with Bowtie 2 and only a sorted, indexed BAM alignment is produced.

## Requirements

- [Nextflow](https://www.nextflow.io/) (DSL2)
- A container engine: [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/singularity/)/[Apptainer](https://apptainer.org/)

All processes run from published containers (GMAP/GSNAP, Bowtie 2, samtools, bedtools, and `bioperl/bioperl`) — no local Dockerfile is built for this repo.

Execution profiles are provided under `conf/`:
- `conf/docker.config` — run with Docker (default, included via `nextflow.config`)
- `conf/singularity.config` — run with Singularity
- `conf/lsf.config` — run with Singularity on an LSF cluster

## Usage

The pipeline has a single, unnamed entry point (no `-entry` flag is needed); behavior is switched by `params.platformType`.

Spliced/expression probes, generating a gene-to-probe mapping and a CDF file:

```bash
nextflow run VEuPathDB/array-probes-nextflow -r main -resume \
  --platformType expression \
  --genomeFastaFile /path/to/genome.fasta \
  --gtfFile         /path/to/genome.gtf \
  --probesFastaFile /path/to/probes.fsa \
  --vendorMappingFile /path/to/vendor.cdf \
  --makeCdfFile true \
  --outputDir /path/to/output \
  -C conf/docker.config
```

Genomic/unspliced probes, producing only a sorted BAM of probe alignments:

```bash
nextflow run VEuPathDB/array-probes-nextflow -r main -resume \
  --platformType genomic \
  --genomeFastaFile /path/to/genome.fasta \
  --probesFastaFile /path/to/probes.fsa \
  --outputDir /path/to/output \
  -C conf/docker.config
```

## Key Parameters

| Parameter               | Default                | Description |
| ------------------------ | ----------------------- | ------------ |
| `platformType`            | `expression`             | `expression` for spliced/mRNA probes (GSNAP + splice-aware mapping); any other value uses Bowtie 2 for straight genomic mapping. |
| `genomeFastaFile`         | `data/genome.fasta`      | Genome FASTA the probes are mapped against. |
| `probesFastaFile`         | `data/probes.fsa`         | FASTA of probe sequences to map (required). |
| `gtfFile`                 | `data/genome.gtf`         | Gene annotation used for splice sites (GSNAP) and for intersecting probe alignments with genes. |
| `wantSplicedAlignments`   | `true`                   | Whether GSNAP should use known splice sites from `gtfFile` when mapping expression probes. |
| `outputDir`               | `output`                 | Directory results are published to. |
| `outputFileName`          | `probes.bam`              | Name of the sorted, indexed BAM of probe alignments. |
| `outputMappingFileName`   | `geneProbeMapping.tab`    | Name of the generated probe-to-gene mapping file (expression probes only). |
| `makeCdfFile`             | `false`                  | Whether to regenerate an Affymetrix CDF platform file from the new gene-to-probe mapping. |
| `makeNdfFile`             | `true`                   | Whether to regenerate a NimbleGen NDF platform file from the new gene-to-probe mapping. |
| `vendorMappingFile`       | —                        | Path to the existing vendor platform file (CDF/NDF) used as a template when regenerating platform files. |
| `arrayRows` / `arrayColumns` | `712` / `712`         | Physical array dimensions used when constructing a new CDF header. |

## Output

- `outputDir/<outputFileName>` (and `.bai` index) — sorted, indexed BAM of probe alignments to the genome.
- For `platformType = expression`:
  - `outputDir/<outputMappingFileName>` — tab-delimited gene-to-probe mapping derived from intersecting probe alignments with the GTF annotation.
  - A regenerated Affymetrix CDF file (if `makeCdfFile` is true) and/or a regenerated NimbleGen NDF file (if `makeNdfFile` is true), reflecting the updated probe-to-gene assignments.
