# WES → ClinVar Variant Report Pipeline

A script-based pipeline for **whole exome sequencing (WES)** data: align paired-end
Illumina reads, call germline variants, annotate against **ClinVar**, and generate
a clinical-style report — variants bucketed into **Pathogenic / Likely Pathogenic**,
**Uncertain Significance (VUS)**, **Benign / Likely Benign**, **Other**
(drug-response / risk-factor), and **Not in ClinVar** — each with the associated
condition, ClinVar review confidence, and a plain-language explanation.

Works with **any** reference genome + WES sample, not tied to a specific dataset —
just point the scripts at your own files (see [Directory structure](#directory-structure) below).

> **Research/educational use only.** Not validated as a clinical diagnostic
> pipeline. Any Pathogenic/Likely Pathogenic finding must be confirmed by an
> accredited clinical laboratory and interpreted by a qualified clinician or
> genetic counselor (e.g. under ACMG/AMP guidelines) before being acted on.

---

## Table of contents
- [Workflow overview](#workflow-overview)
- [Supported inputs](#supported-inputs)
- [Directory structure](#directory-structure)
- [Tools & databases](#tools--databases)
- [Installation](#installation)
- [Usage](#usage)
- [Commands & parameters explained](#commands--parameters-explained)
- [Sample output report](#sample-output-report)
- [Troubleshooting](#troubleshooting)
- [Disclaimer](#disclaimer)

---

## Workflow overview

```
 Raw FASTQ (R1 + R2, Illumina paired-end WES)
        │
        │  fastp — adapter trim + quality filter
        ▼
 Trimmed FASTQ
        │  BWA-MEM — align to reference exome/genome
        │  samtools sort/index
        ▼
 Sorted BAM
        │  GATK MarkDuplicates — flag PCR/optical duplicates
        ▼
 Deduplicated BAM
        │  GATK HaplotypeCaller — germline SNV/indel calling
        │  (optionally restricted to capture-kit target BED)
        ▼
 Raw VCF
        │  GATK VariantFiltration — hard-filter low-confidence calls
        │  bcftools norm — left-align & split multiallelic sites
        ▼
 Filtered, normalized VCF
        │  SnpSift annotate — join against ClinVar by position
        ▼
 ClinVar-annotated VCF
        │  02_generate_report.py — bucket + explain
        ▼
 Markdown / HTML / TSV clinical report
```

## Supported inputs

| Input | Supported formats | Notes |
|---|---|---|
| Sequencing platform | Illumina paired-end WES (HiSeq/NovaSeq/NextSeq) | Not designed for long-read (PacBio/ONT) or single-end data |
| Reads | `.fastq`, `.fastq.gz` | Must be a matched R1/R2 pair from the same sample/library |
| Reference genome | `.fa` / `.fasta`, uncompressed | GRCh38 (recommended) or GRCh37/hg19 — must match the ClinVar build you download |
| Capture kit target regions (optional) | `.bed` | e.g. Agilent SureSelect, IDT xGen Exome, Twist Exome — restricts calling/reporting to on-target exons |
| Variant database | ClinVar VCF (`vcf_GRCh38` or `vcf_GRCh37`) | Auto-downloaded and chromosome-renamed by `00b_get_clinvar.sh` |

## Directory structure

Any layout works as long as you pass correct paths to `-r/-1/-2/-c/-s`, but a
typical project looks like:

```
my_project/
├── scripts/
│   ├── 00_install_tools.sh
│   ├── 00b_get_clinvar.sh
│   ├── 01_pipeline.sh
│   └── 02_generate_report.py
├── reference/
│   ├── genome.fa                  # your reference FASTA (any build)
│   └── genome.fa.{fai,bwt,...}    # built automatically on first run if missing
├── capture_kit/
│   └── targets.bed                # optional: your exome kit's target BED
├── clinvar/
│   └── clinvar_exome.vcf.gz       # produced by 00b_get_clinvar.sh
├── raw_reads/
│   ├── <SAMPLE>_R1.fastq.gz
│   └── <SAMPLE>_R2.fastq.gz
└── results/
    └── <SAMPLE>/
        ├── qc/       (fastp reports)
        ├── align/    (BAMs, flagstat, dup metrics)
        ├── vcf/      (raw, filtered, normalized, ClinVar-annotated VCFs)
        └── report/   (final .md / .html / .tsv report)
```

Swap in **your own** `genome.fa`, capture BED, and `<SAMPLE>_R1/R2.fastq.gz` —
nothing in the scripts is hardcoded to a specific sample or chromosome anymore.

## Tools & databases

| Tool | Role | Version used | Link |
|---|---|---|---|
| [fastp](https://github.com/OpenGene/fastp) | Read QC + adapter trimming | 0.23.4 | github.com/OpenGene/fastp |
| [BWA](https://github.com/lh3/bwa) | Short-read alignment (BWA-MEM) | 0.7.17 | github.com/lh3/bwa |
| [SAMtools](http://www.htslib.org/) | BAM sort/index/stats | 1.19 | htslib.org |
| [GATK4](https://github.com/broadinstitute/gatk) | Duplicate marking, HaplotypeCaller, hard-filtering | 4.5.0.0 | gatk.broadinstitute.org |
| [BCFtools](http://www.htslib.org/) | VCF normalization, filtering | 1.19 | htslib.org |
| [SnpSift](https://pcingola.github.io/SnpEff/) (SnpEff suite) | ClinVar annotation | 5.2 | pcingola.github.io/SnpEff |
| [BEDtools](https://bedtools.readthedocs.io/) | Target-region handling | 2.31.1 | bedtools.readthedocs.io |
| [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) | Pathogenicity/clinical significance database | rolling (NCBI weekly release) | ncbi.nlm.nih.gov/clinvar · [VCF FTP](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/) |
| Reference genome | GRCh38 (or GRCh37) assembly | build-dependent | [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/) · [UCSC](https://hgdownload.soe.ucsc.edu/downloads.html) · [Ensembl](https://www.ensembl.org/info/data/ftp/index.html) |

ClinVar updates weekly — re-run `00b_get_clinvar.sh` periodically if you want
current classifications; a variant's `CLNSIG` label (e.g. VUS → Likely
Pathogenic) can and does change as new evidence is submitted.

## Installation

```bash
chmod +x scripts/*.sh scripts/*.py
./scripts/00_install_tools.sh wes-pipeline   # creates a conda env named 'wes-pipeline'
conda activate wes-pipeline                  # (or: mamba activate wes-pipeline)
```

This installs all seven tools above with pinned versions in one isolated
environment. Verify:

```bash
fastp --version
bwa 2>&1 | head -3
samtools --version | head -1
bcftools --version | head -1
gatk --version
SnpSift 2>&1 | head -5
bedtools --version
```

## Usage

**1. Fetch ClinVar, matched to your reference's chromosome naming**
```bash
./scripts/00b_get_clinvar.sh reference/genome.fa clinvar/ GRCh38
```

**2. Run the pipeline on your sample**
```bash
./scripts/01_pipeline.sh \
  -r reference/genome.fa \
  -1 raw_reads/<SAMPLE>_R1.fastq.gz \
  -2 raw_reads/<SAMPLE>_R2.fastq.gz \
  -c clinvar/clinvar_exome.vcf.gz \
  -s <SAMPLE> \
  -t 8 \
  -b capture_kit/targets.bed \
  -o results
```
(`-b` and `-t` are optional — omit `-b` to call variants across the whole
reference rather than restricting to a capture kit's target regions.)

**3. Generate the report**
```bash
python scripts/02_generate_report.py \
  results/<SAMPLE>/vcf/<SAMPLE>.clinvar_annotated.vcf \
  <SAMPLE> \
  results/<SAMPLE>/report
```

## Commands & parameters explained

**Trimming (`fastp`)**
| Flag | Meaning |
|---|---|
| `--detect_adapter_for_pe` | Auto-detects Illumina adapter sequence from paired-end read overlap, no need to specify it manually |
| `-q 20` | Minimum Phred base quality to count as "high quality" (Q20 = 1% error rate) |
| `-l 50` | Discard trimmed reads shorter than 50bp |
| `-5 -3` | Sliding-window quality trimming from both the 5' and 3' ends |

**Alignment (`bwa mem`)**
| Flag | Meaning |
|---|---|
| `-t` | Threads |
| `-R` | Read group string (`SM` = sample name — must be set correctly, GATK uses it downstream) |

**Variant calling (`gatk HaplotypeCaller`)**
| Flag | Meaning |
|---|---|
| `--sample-ploidy 2` | Diploid calling (standard for autosomal/human germline) |
| `-L <bed>` | Restrict calling to capture-kit target intervals |
| `-ip 20` | Pad 20bp around each target interval, catches variants right at exon boundaries |

**Hard filtering (`gatk VariantFiltration`)**
| Filter | Meaning | Typical rationale |
|---|---|---|
| `QD < 2.0` | Variant confidence normalized by depth | Below 2.0 = weak signal relative to coverage, common artifact threshold |
| `FS > 60.0` | Fisher Strand bias (Phred-scaled) | High = reads supporting the variant are lopsided by strand, suggests artifact |
| `MQ < 40.0` | Root-mean-square mapping quality | Low = reads map ambiguously (repeats, paralogs) |
| `DP < 10` | Read depth at the site | Below 10x, genotype calls are statistically unreliable |

These are GATK's standard germline hard-filter thresholds — reasonable
defaults for WES, but **exome capture edges naturally have lower/variable
depth**, so loosen `DP`/`QD` if you're seeing too many real variants dropped
near target boundaries (check with the filter-breakdown command in
[Troubleshooting](#troubleshooting)).

**Normalization (`bcftools norm -m -both`)**
Splits multiallelic sites into one record per alt allele and left-aligns
indels — required so variant positions line up exactly with ClinVar's
representation; without this step, real matches get silently missed.

## Sample output report

Given an annotated VCF with 5 example variants, `02_generate_report.py` produces:

```markdown
# Clinical Variant Report — DEMO_SAMPLE

Total variants processed: 5
- Pathogenic / Likely Pathogenic: 2
- Uncertain Significance (VUS): 1
- Other (drug response / risk factor / association): 0
- Benign / Likely Benign: 1
- Not in ClinVar: 1

## Pathogenic / Likely Pathogenic (2)

### KCNE1 — 21:5033832 C>T
- HGVS: NC_000021.9:g.5033832C>T
- ClinVar significance: Pathogenic
- Associated condition(s): Long QT syndrome
- Review status: criteria provided, multiple submitters, no conflicts
- Explanation: ClinVar submitters agree this variant causes or strongly
  predisposes to the associated disease, based on clinical, functional,
  and/or segregation evidence.

### SOD1 — 21:16340889 G>A
- HGVS: NC_000021.9:g.16340889G>A
- ClinVar significance: Likely_pathogenic
- Associated condition(s): Amyotrophic lateral sclerosis
- Review status: criteria provided, single submitter
- Explanation: Evidence favors disease causation (>90% confidence per
  ACMG/AMP guidelines) but does not yet meet the full threshold for
  'Pathogenic'.

## Uncertain Significance (VUS) (1)

### RUNX1 — 21:34196343 A>G
- ClinVar significance: Uncertain_significance
- Associated condition(s): Platelet disorder
- Explanation: Available evidence is insufficient or conflicting to
  classify this variant as disease-causing or benign (a 'VUS').

## Benign / Likely Benign (1)

### COL6A1 — 21:43053191 T>C
- ClinVar significance: Benign
- Explanation: This variant is common and/or has strong evidence showing
  it does not cause the associated disease.

## Not in ClinVar (1)

### 21:27543384 G>T
- Explanation: No ClinVar classification available for this variant.
```

The `.html` version renders the same content as color-coded cards
(red = Pathogenic, orange = VUS, purple = Other, green = Benign, grey = Not
in ClinVar), and the `.tsv` gives every variant as one flat row for
filtering in Excel or pandas.

## Troubleshooting

**0 variants annotated, but variants were called**
Chromosome naming mismatch (`chr21` vs `21`) between reference and ClinVar.
`00b_get_clinvar.sh` auto-detects and renames — confirm afterward with:
```bash
grep "^>" reference/genome.fa | head -1
zcat clinvar/clinvar_exome.vcf.gz | grep -v "^#" | cut -f1 | sort -u | head
```
These should use the same naming convention.

**0 variants called at all**
Check alignment first — near-0% mapped in `flagstat.txt` means the reads
don't match the reference provided (wrong build, wrong sample, wrong species).

**`samtools sort` gets `Killed` mid-run**
Out-of-memory kill. `01_pipeline.sh` already caps `samtools sort` at
`-m 512M` per thread — if it still happens, lower `-t` (threads) and check
`free -h` before rerunning.

**Very few variants survive filtering**
```bash
zcat results/<SAMPLE>/vcf/<SAMPLE>.filtered.vcf.gz | grep -v "^#" | cut -f7 | sort | uniq -c
```
Shows which hard filter (`lowQD`/`highFS`/`lowMQ`/`lowDP`) is removing the
most variants, so you know which threshold to relax for your data's actual
coverage/quality profile.

## Disclaimer

ClinVar classifications reflect the database snapshot at the time you ran
`00b_get_clinvar.sh` and change as new clinical evidence is submitted. This
pipeline reports ClinVar's existing classifications only — it does not
perform independent ACMG/AMP variant classification. Findings should be
confirmed by an accredited clinical laboratory before any clinical action
is taken.
