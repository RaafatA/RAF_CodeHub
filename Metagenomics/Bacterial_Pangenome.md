# Bacterial Pan-Genome Analysis Skill

## Overview
Analyze pan-genome across multiple bacterial isolates to identify core, accessory, and unique genes.

## Triggers
Use when user asks for: pan-genome, pangenome, core genome, accessory genome, compare isolates, roary, gene presence absence, multiple isolates comparison

## Required Tools
megahit, prokka, roary

Optional: fasttree (for phylogeny), matplotlib/R (for plots)

## Input Requirements

| Input | Required | Description |
|-------|----------|-------------|
| Multiple isolates | **YES** | ≥3 isolates (paired-end reads each) |
| Same species | **YES** | Pan-genome requires related organisms |

**Minimum:** 3 isolates for meaningful analysis

## Workflow

### 1. Assembly (per isolate)
```bash
mkdir -p assemblies

for r1 in reads/*_1.fastq.gz; do
    prefix=$(basename $r1 _1.fastq.gz)
    r2=reads/${prefix}_2.fastq.gz
    
    megahit -1 $r1 -2 $r2 \
        -o megahit_${prefix} \
        --out-prefix ${prefix} \
        --min-contig-len 500
    
    mv megahit_${prefix}/${prefix}.contigs.fa assemblies/
    rm -r megahit_${prefix}
done
```

### 2. Annotation (per isolate)
```bash
mkdir -p annotation

for assembly in assemblies/*.contigs.fa; do
    prefix=$(basename $assembly .contigs.fa)
    
    prokka --outdir prokka_${prefix} \
        --prefix ${prefix} \
        --kingdom Bacteria \
        --fast \
        --cpus 4 \
        $assembly
    
    cp prokka_${prefix}/${prefix}.gff annotation/
    rm -r prokka_${prefix}
done
```

### 3. Pan-Genome Analysis (Roary)
```bash
roary -f roary_output \
    -e -n \
    -p 4 \
    annotation/*.gff
```

**Key flags:**
- `-e`: Create core gene alignment
- `-n`: Fast core gene alignment with MAFFT
- `-p`: Threads

## Pan-Genome Categories

| Category | Definition | Typical % |
|----------|------------|-----------|
| Core | 100% of isolates | 10-30% |
| Soft-core | 95-99% of isolates | 5-10% |
| Shell | 15-95% of isolates | 20-40% |
| Cloud | <15% of isolates | 30-60% |

## Required Outputs

| Output | Path | Description |
|--------|------|-------------|
| Gene matrix | `roary_output/gene_presence_absence.csv` | Main result |
| Statistics | `roary_output/summary_statistics.txt` | Counts |
| Core alignment | `roary_output/core_gene_alignment.aln` | For phylogeny |

## Optional Visualizations

**If user requests plots:**

### Frequency Histogram
```python
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('roary_output/gene_presence_absence.csv', low_memory=False)
n_isolates = len(df.columns) - 14  # First 14 cols are metadata

# Count presence per gene
presence = df.iloc[:, 14:].notna().sum(axis=1)

plt.figure(figsize=(10, 6))
plt.hist(presence, bins=n_isolates, edgecolor='black')
plt.xlabel('Number of Isolates')
plt.ylabel('Number of Genes')
plt.title('Gene Frequency Distribution')
plt.savefig('gene_frequency.png', dpi=150)
```

### Pie Chart
```python
# Categorize genes
core = (presence == n_isolates).sum()
soft_core = ((presence >= 0.95*n_isolates) & (presence < n_isolates)).sum()
shell = ((presence >= 0.15*n_isolates) & (presence < 0.95*n_isolates)).sum()
cloud = (presence < 0.15*n_isolates).sum()

plt.figure(figsize=(8, 8))
plt.pie([core, soft_core, shell, cloud], 
        labels=['Core', 'Soft-core', 'Shell', 'Cloud'],
        autopct='%1.1f%%')
plt.title('Pan-genome Composition')
plt.savefig('pangenome_pie.png', dpi=150)
```

### Phylogenetic Tree (from core alignment)
```bash
# If roary ran with -e flag
FastTree -nt roary_output/core_gene_alignment.aln > core_tree.nwk
```

## Troubleshooting

| Issue | Cause | Solution |
|-------|-------|----------|
| Too few core genes | Different species mixed | Verify all same species |
| Roary fails | GFF format issues | Check Prokka completed |
| No alignment | Missing -e flag | Rerun with `-e -n` |

## Critical Rules

1. **MINIMUM 3 isolates** - Pan-genome needs multiple genomes
2. **SAME SPECIES** - Don't mix unrelated bacteria
3. **USE -e flag** if phylogenetic tree needed
4. **PROKKA annotation required** - Roary needs GFF3 format

## Notes

- Roary uses 95% BLASTP identity by default (change with `-i`)
- Core alignment useful for phylogenetic tree
- For >100 isolates, consider using `-cd 100` for faster clustering
- Memory: ~1GB per isolate
