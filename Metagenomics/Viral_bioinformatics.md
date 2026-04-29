# Viral Metagenomics Analysis Skill

## Overview
Identify and characterize viral sequences from metagenomic data. Supports both viral species identification (clinical/diagnostic) and novel virus discovery (virome profiling).

## Triggers
Use when user asks for: viral metagenomics, virome, virsorter, checkv, find viruses, bacteriophage, phage discovery, viral contigs, vOTUs, viral agents, viral species identification

## Available Tools
megahit, bowtie2, samtools, fastp, kaiju, virsorter2, checkv, diamond, drep, spades (metaspades)

## Input Requirements

| Input | Required | Description |
|-------|----------|-------------|
| Metagenomic reads | **YES** | Paired-end FASTQ |
| Host organism | **YES** | Needed for host removal — ask user if not obvious |

---

## Choosing the Right Workflow

| Goal | Workflow | When to use |
|------|----------|-------------|
| **Identify viral species** | Host removal → Assembly → Kaiju classification | Clinical samples, "what viruses are present", species-level ID |
| **Virome discovery** | Host removal → Assembly → VirSorter2 → CheckV | Novel phage discovery, vOTU profiling, virome characterization |

Both workflows share steps 1-3 (QC, host removal, assembly), then diverge.

---

## Shared Steps (Both Workflows)

### 1. QC & Trimming

```bash
mkdir -p trimmed

fastp -i reads/sample_R1.fastq.gz -I reads/sample_R2.fastq.gz \
    -o trimmed/sample_R1.fastq -O trimmed/sample_R2.fastq \
    --detect_adapter_for_pe --length_required 30 \
    --cut_front --cut_tail --cut_mean_quality 10 \
    -w 4 --html trimmed/fastp_report.html
```

### 2. Host Removal (MANDATORY)

**CRITICAL: Always remove host reads before assembly.** Without this step, >90% of reads may be host genome, burying viral signal and producing host contigs that waste compute and confuse classification.

**Step A: Get host reference genome**

| Host | Reference | Source |
|------|-----------|--------|
| Human | GRCh38 | Ensembl or NCBI |
| Mouse | GRCm39 | Ensembl |
| Dolphin | turTru1 | `ftp.ensembl.org/pub/release-113/fasta/tursiops_truncatus/dna/` |
| Chicken | GRCg7b | Ensembl |
| Other | Search Ensembl/NCBI for species | Use top-level DNA FASTA |

```bash
mkdir -p reference

# Download host genome (example: dolphin)
wget -O reference/host_genome.fa.gz [URL_TO_HOST_GENOME]
gunzip reference/host_genome.fa.gz
```

**Step B: Build index and remove host reads**
```bash
mkdir -p host_removal

# Build bowtie2 index
bowtie2-build reference/host_genome.fa reference/host_index --threads 4

# Map reads to host — keep UNMAPPED reads
bowtie2 -x reference/host_index \
    -1 trimmed/sample_R1.fastq \
    -2 trimmed/sample_R2.fastq \
    -S host_removal/mapped.sam \
    --un-conc host_removal/unmapped_%.fastq \
    --threads 4

# unmapped_1.fastq and unmapped_2.fastq are now host-free
echo "Host removal stats:"
samtools flagstat host_removal/mapped.sam | head -5
```

The `--un-conc` flag outputs read pairs where NEITHER mate mapped to host. These are your non-host reads for assembly.

**If host genome is not available**, skip this step but note that results will contain host contamination and viral sensitivity will be lower.

### 3. Assembly

**Assemble HOST-REMOVED reads, not raw reads:**

```bash
mkdir -p assembly

megahit \
    -1 host_removal/unmapped_1.fastq \
    -2 host_removal/unmapped_2.fastq \
    -o assembly/megahit_out \
    --min-contig-len 1000 \
    --num-cpu-threads 4

# Check assembly stats
echo "Contigs produced:"
grep -c ">" assembly/megahit_out/final.contigs.fa
```

**Use `--min-contig-len 1000`** — viral genomes need longer contigs for reliable classification.

| Assembler | When to use |
|-----------|-------------|
| MEGAHIT | Default — fast, memory-efficient |
| metaSPAdes | If MEGAHIT produces few contigs — better for complex samples |

---

## Workflow A: Viral Species Identification

Use when the task asks to identify what viruses are present, classify viral agents, or report species-level taxonomy.

### 4A. Kaiju Classification (Viral DB)

**Kaiju classifies contigs against a viral protein database** — more sensitive than nucleotide-based classifiers for divergent viruses.

```bash
mkdir -p classification

# Download Kaiju viral database (if not available)
wget https://kaiju-idx.s3.eu-central-1.amazonaws.com/2024/kaiju_db_viruses_2024-08-15.tgz
tar xzf kaiju_db_viruses_2024-08-15.tgz
rm kaiju_db_viruses_2024-08-15.tgz

# Classify assembled contigs
kaiju -t nodes.dmp \
    -f kaiju_db_viruses.fmi \
    -i assembly/megahit_out/final.contigs.fa \
    -o classification/kaiju_output.tsv \
    -z 4 -v

# Add taxonomic names
kaiju-addTaxonNames \
    -t nodes.dmp \
    -n names.dmp \
    -i classification/kaiju_output.tsv \
    -o classification/kaiju_names.tsv \
    -r superkingdom,phylum,class,order,family,genus,species

# Generate Krona-format file (viewable if KronaTools installed)
kaiju2krona -t nodes.dmp \
    -n names.dmp \
    -i classification/kaiju_output.tsv \
    -o classification/kaiju.krona \
    -u
```

### 5A. Build Output Table

```python
import pandas as pd
from collections import Counter

# Parse Kaiju output
classified = []
unclassified = 0

with open('classification/kaiju_names.tsv') as f:
    for line in f:
        parts = line.strip().split('\t')
        status = parts[0]  # C=classified, U=unclassified
        if status == 'U':
            unclassified += 1
        elif status == 'C':
            # Last field contains taxonomy string
            tax = parts[-1] if len(parts) > 3 else ''
            # Extract species (last non-empty level)
            levels = [l.strip() for l in tax.split(';') if l.strip()]
            species = levels[-1] if levels else 'Unknown'
            domain = levels[0] if levels else 'Unknown'
            classified.append({'domain': domain, 'species': species})

# Count contigs per species
df = pd.DataFrame(classified)
if len(df) > 0:
    counts = df.groupby(['domain', 'species']).size().reset_index(name='contig_count')
else:
    counts = pd.DataFrame(columns=['domain', 'species', 'contig_count'])

# Add unclassified
if unclassified > 0:
    unclass_row = pd.DataFrame([{
        'contig_count': unclassified,
        'domain': '',
        'species': 'Unclassified'
    }])
    counts = pd.concat([counts, unclass_row], ignore_index=True)

# Reorder columns to match expected format
counts = counts[['contig_count', 'domain', 'species']]
counts = counts.sort_values('contig_count', ascending=False)

counts.to_csv('results/taxonomy.csv', index=False)
print(counts.to_string(index=False))
```

**Expected output format:**
```
contig_count,domain,species
44,Viruses,Bottlenose dolphin adenovirus 1
2,Viruses,Equid gammaherpesvirus 5
68,,Unclassified
```

---

## Workflow B: Virome Discovery (vOTU Profiling)

Use when the task asks for virome characterization, phage discovery, or novel virus identification.

### 4B. Viral Identification (VirSorter2)
```bash
mkdir -p virsorter2_output

virsorter run \
    -w virsorter2_output/sample \
    -i assembly/megahit_out/final.contigs.fa \
    --min-length 1000 \
    --min-score 0.5 \
    -j 4 \
    all

cp virsorter2_output/sample/final-viral-combined.fa virsorter2_output/sample_viral.fa
```

| Setting | --min-score | Use case |
|---------|-------------|----------|
| Default | 0.5 | Balanced |
| Stringent | 0.9 | High confidence |
| Sensitive | 0.3 | Maximum discovery |

### 5B. Quality Assessment (CheckV)
```bash
mkdir -p checkv_output

checkv end_to_end virsorter2_output/sample_viral.fa checkv_output/sample -t 4
```

| Quality Tier | Completeness | Description |
|-------------|--------------|-------------|
| Complete | 100% | Circular or DTR |
| High | >90% | Near-complete |
| Medium | 50-90% | Partial |
| Low | <50% | Fragment |

### 6B. Dereplication (Multi-sample)

Only if multiple samples — cluster at 95% ANI using nucmer/ANI comparison:
```bash
cat checkv_output/*/viruses.fna > all_viruses.fna

# Use dRep for dereplication at 95% ANI
# First, create per-genome FASTAs in a directory
mkdir -p derep_input
awk '/^>/{f="derep_input/"++n".fa"} {print > f}' all_viruses.fna

dRep dereplicate derep_output \
    -g derep_input/*.fa \
    -sa 0.95 \
    -comp 0 -con 1000 \
    --ignoreGenomeQuality \
    -p 4

cat derep_output/dereplicated_genomes/*.fa > viral_representatives.fna
```

### 7B. Taxonomic Classification of vOTUs

```bash
kaiju -t nodes.dmp \
    -f kaiju_db_viruses.fmi \
    -i viral_representatives.fna \
    -o taxonomy/viral_classification.tsv \
    -z 4 -v

kaiju-addTaxonNames \
    -t nodes.dmp -n names.dmp \
    -i taxonomy/viral_classification.tsv \
    -o taxonomy/viral_classification_names.tsv
```

### 8B. Abundance Quantification
```bash
mkdir -p abundance

bowtie2-build viral_representatives.fna viral_index

for r1 in reads/*_R1.fastq.gz; do
    prefix=$(basename $r1 _R1.fastq.gz)
    r2=reads/${prefix}_R2.fastq.gz

    bowtie2 -x viral_index \
        -1 $r1 -2 $r2 \
        -p 4 --no-unal \
        | samtools view -bS -F 4 - \
        | samtools sort -o abundance/${prefix}.sorted.bam -@ 4

    samtools index abundance/${prefix}.sorted.bam
    samtools coverage abundance/${prefix}.sorted.bam > abundance/${prefix}_coverage.tsv
done
```

---

## Multi-Sample Batch Processing

```bash
SAMPLES=("sample1" "sample2" "sample3")
THREADS=4

for sample in "${SAMPLES[@]}"; do
    echo "=== Processing $sample ==="

    # Trim
    fastp -i reads/${sample}_R1.fastq.gz -I reads/${sample}_R2.fastq.gz \
        -o trimmed/${sample}_R1.fastq -O trimmed/${sample}_R2.fastq \
        --detect_adapter_for_pe --length_required 30 -w $THREADS

    # Host removal
    bowtie2 -x reference/host_index \
        -1 trimmed/${sample}_R1.fastq -2 trimmed/${sample}_R2.fastq \
        --un-conc host_removal/${sample}_unmapped_%.fastq \
        -S /dev/null --threads $THREADS

    # Assembly
    megahit -1 host_removal/${sample}_unmapped_1.fastq \
            -2 host_removal/${sample}_unmapped_2.fastq \
            -o assembly/${sample}/ --min-contig-len 1000 -t $THREADS

    # Classify
    kaiju -t nodes.dmp -f kaiju_db_viruses.fmi \
        -i assembly/${sample}/final.contigs.fa \
        -o classification/${sample}_kaiju.tsv -z $THREADS
done
```

---

## Required Outputs

### Workflow A (Species Identification)
| Output | Description |
|--------|-------------|
| `results/taxonomy.csv` | Contig counts per domain/species + unclassified |
| `classification/kaiju_names.tsv` | Full Kaiju classification |
| `classification/kaiju.krona` | Krona visualization |

### Workflow B (Virome Discovery)
| Output | Description |
|--------|-------------|
| `virsorter2_output/*_viral.fa` | Identified viral contigs |
| `checkv_output/*/quality_summary.tsv` | Quality assessment |
| `viral_representatives.fna` | Dereplicated vOTUs |
| `taxonomy/viral_classification_names.tsv` | vOTU taxonomy |
| `abundance/*_coverage.tsv` | Per-sample abundance |

---

## Database Requirements

| Tool | Database | Size | Download |
|------|----------|------|----------|
| Kaiju | kaiju_db_viruses | ~10GB | `kaiju-idx.s3.eu-central-1.amazonaws.com` |
| VirSorter2 | Built-in | ~12GB | `virsorter setup -d vs2_db/` |
| CheckV | checkv-db | ~2GB | `checkv download_database /path/` |
| DIAMOND | viral_refseq | Variable | NCBI |

---

## Quality Thresholds

| Metric | Good | Investigate |
|--------|------|-------------|
| Host removal rate | 80-99% mapped to host | <50% (wrong host reference?) |
| Contigs assembled | >50 | <10 (low viral load or poor assembly) |
| Classified contigs | >10% | <5% (novel viruses or wrong DB) |
| Unclassified contigs | Common (30-80%) | Many novel viruses expected |

---

## Edge Cases & Troubleshooting

| Situation | Solution |
|-----------|----------|
| Don't know the host organism | Ask user; if clinical sample, check sample metadata |
| Host genome not available | Skip host removal, but note lower sensitivity |
| Very few contigs after assembly | Try metaSPAdes instead of MEGAHIT; lower `--min-contig-len` to 500 |
| Kaiju DB not available | Download from kaiju-idx.s3.eu-central-1.amazonaws.com |
| All contigs unclassified | Normal for novel viruses — try VirSorter2 workflow for functional characterization |
| Task asks "identify viral agents" | Use Workflow A (Kaiju classification) |
| Task asks "virome" or "phage discovery" | Use Workflow B (VirSorter2 + CheckV) |
| Human clinical sample | Host removal against GRCh38 is mandatory |
| Environmental sample (water/soil) | Host removal may not be needed if no obvious host |
| VirSorter2 finds no viral contigs | Try `--min-score 0.3`; or contigs may not be viral |
| Task asks for contig counts | Report classified + unclassified totals |

## Critical Rules

1. **ALWAYS REMOVE HOST READS** before assembly (unless environmental sample with no host)
2. **ASSEMBLE HOST-REMOVED READS** — never assemble raw reads for viral metagenomics
3. **USE VIRAL-SPECIFIC DB** — Kaiju viral DB, not general Kraken2 standard
4. **USE `--min-contig-len 1000`** — viral genomes need longer contigs
5. **REPORT UNCLASSIFIED CONTIGS** — many viral contigs won't classify (novel viruses)
6. **RUN CheckV** for virome discovery — essential for quality assessment
