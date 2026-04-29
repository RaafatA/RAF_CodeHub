# Phylogenetic Tree Construction Skill

## Overview
Build phylogenetic trees from DNA or protein sequences: alignment, trimming, tree inference, and visualization.

## Triggers
Use when user asks for: phylogenetic tree, phylogeny, evolutionary tree, build tree, mafft, iqtree, raxml, fasttree, sequence alignment, newick, tree visualization

## Required Tools
mafft, trimal, iqtree (or raxml-ng, fasttree)

Python: ete3 (visualization)

## Input Requirements

| Input | Required | Description |
|-------|----------|-------------|
| Sequences | **YES** | FASTA file (DNA or protein) |
| Minimum | **≥4 sequences** | Tree needs multiple taxa |
| Same type | **YES** | All DNA or all protein |

## Workflow

### 1. Input Validation
```bash
# Check sequence count and type
seqmagick info input.fasta
```

**Verify:**
- ≥4 sequences
- No duplicate IDs
- All same type (DNA or protein)

### 2. Multiple Sequence Alignment (MAFFT)

| Dataset Size | Command |
|--------------|---------|
| <200 seqs | `mafft --auto input.fasta > aligned.fasta` |
| 200-1000 seqs | `mafft --maxiterate 1000 --localpair input.fasta > aligned.fasta` |
| >1000 seqs | `mafft --retree 1 --maxiterate 0 input.fasta > aligned.fasta` |

```bash
# Default (auto-select algorithm)
mafft --auto input.fasta > aligned.fasta
```

### 3. Alignment Trimming (trimAl)

| Mode | Command | Use case |
|------|---------|----------|
| Gappyout | `trimal -gappyout` | Default, removes gappy columns |
| Automated | `trimal -automated1` | Heuristic selection |
| Strict | `trimal -nogaps` | No gaps allowed |

```bash
trimal -in aligned.fasta -out trimmed.fasta -gappyout
```

### 4. Tree Inference

**IQ-TREE (recommended):**
```bash
iqtree -s trimmed.fasta \
    -m MFP \
    -bb 1000 \
    -alrt 1000 \
    -nt 4 \
    --prefix output
```

| Flag | Purpose |
|------|---------|
| `-m MFP` | Auto model selection |
| `-bb 1000` | Ultrafast bootstrap |
| `-alrt 1000` | SH-aLRT test |
| `-nt AUTO` | Auto threads |

**RAxML-NG (alternative):**
```bash
raxml-ng --all \
    --msa trimmed.fasta \
    --model GTR+G \
    --bs-trees 100 \
    --threads 4 \
    --prefix output
```

**FastTree (quick/large datasets):**
```bash
# DNA
FastTree -gtr -nt < trimmed.fasta > fasttree.nwk

# Protein
FastTree < trimmed.fasta > fasttree.nwk
```

### 5. Tree Rooting (ETE3)

**Midpoint rooting (default):**
```python
from ete3 import Tree

tree = Tree("output.treefile")
tree.set_outgroup(tree.get_midpoint_outgroup())
tree.write(outfile="output_rooted.nwk")
```

**Outgroup rooting (if specified):**
```python
from ete3 import Tree

tree = Tree("output.treefile")
outgroup = tree.search_nodes(name="outgroup_taxon")[0]
tree.set_outgroup(outgroup)
tree.write(outfile="output_rooted.nwk")
```

### 6. Visualization (ETE3)
```python
from ete3 import Tree, TreeStyle, NodeStyle

tree = Tree("output_rooted.nwk")

ts = TreeStyle()
ts.show_leaf_name = True
ts.show_branch_support = True

# Color nodes by bootstrap support
for node in tree.traverse():
    if not node.is_leaf() and hasattr(node, 'support'):
        nstyle = NodeStyle()
        if node.support >= 95:
            nstyle["fgcolor"] = "darkgreen"
        elif node.support >= 70:
            nstyle["fgcolor"] = "orange"
        else:
            nstyle["fgcolor"] = "red"
        nstyle["size"] = 6
        node.set_style(nstyle)

tree.render("tree.png", w=800, units="px", tree_style=ts)
tree.render("tree.pdf", tree_style=ts)
```

### 7. Tree Statistics
```python
from ete3 import Tree

tree = Tree("output_rooted.nwk")

n_taxa = len(tree.get_leaves())
total_length = sum(n.dist for n in tree.traverse())
supports = [n.support for n in tree.traverse() if not n.is_leaf() and hasattr(n, 'support')]

high_support = sum(1 for s in supports if s >= 95)
low_support = sum(1 for s in supports if s < 70)

print(f"Taxa: {n_taxa}")
print(f"Total branch length: {total_length:.4f}")
print(f"High support nodes (≥95%): {high_support}")
print(f"Low support nodes (<70%): {low_support}")
```

## Bootstrap Interpretation

| Support | Interpretation |
|---------|----------------|
| ≥95% | Strong support |
| 70-94% | Moderate support |
| <70% | Weak support |

## Tool Selection

| Tool | Speed | Accuracy | Use Case |
|------|-------|----------|----------|
| IQ-TREE | Medium | High | Default, publication |
| RAxML-NG | Slow | High | Traditional ML |
| FastTree | Fast | Medium | Quick preview, >1000 seqs |

## Required Outputs

| Output | Path | Description |
|--------|------|-------------|
| Alignment | `aligned.fasta` | MSA |
| Trimmed | `trimmed.fasta` | Cleaned alignment |
| Tree | `output.treefile` | Best ML tree |
| Rooted | `output_rooted.nwk` | Rooted tree |
| Image | `tree.png` | Visualization |
| PDF | `tree.pdf` | Vector format |
| Stats | `tree_statistics.txt` | Summary |

## Troubleshooting

| Issue | Cause | Solution |
|-------|-------|----------|
| MAFFT slow | Large dataset | Use `--retree 1` |
| All gaps after trim | Alignment too divergent | Use `-automated1` or skip trim |
| Low bootstrap | Short sequences or few informative sites | Need more data |
| ETE3 render fails | Display issues | Use `tree.render()` not `tree.show()` |

## Critical Rules

1. **MINIMUM 4 sequences** - Trees need ≥4 taxa
2. **SAME sequence type** - Don't mix DNA and protein
3. **UNIQUE IDs** - Duplicate names cause errors
4. **ALWAYS TRIM** - Remove poorly aligned regions
5. **ROOT the tree** - Midpoint if no outgroup
6. **REPORT bootstrap** - Essential for interpretation

## Model Selection

IQ-TREE `-m MFP` auto-selects, but common models:

| Data | Model |
|------|-------|
| DNA | GTR+G, HKY+G |
| Protein | LG+G, WAG+G, JTT+G |
| Codon | GY+G |

## Notes

- For protein-coding DNA, consider codon-aware alignment
- IQ-TREE output: `.treefile` (best tree), `.contree` (consensus)
- FastTree is ~100x faster but less accurate
- Large alignments (>10k sites): consider `-fast` flag in IQ-TREE
- ETE3 requires X server or use `render()` for headless
