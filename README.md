# Yeast mRNA GC Content Analysis

## Research Question

This analysis examines the distribution of GC content (proportion of guanine and cytosine nucleotides) across yeast mRNA sequences from *Saccharomyces cerevisiae*. Understanding GC content variation can reveal insights into codon usage biases, gene expression patterns, and evolutionary pressures on different gene classes.

## What the Script Does

The `analyze_mrna_gc.py` script:
1. Streams through a gzipped FASTA file of yeast mRNA sequences without loading the entire file into memory
2. Extracts the accession ID from each sequence header
3. Computes the GC content (fraction of G and C nucleotides) for each sequence
4. Outputs a TSV file with accession, sequence length, and GC content (sorted by GC content descending)
5. Generates a histogram plot showing the distribution of GC content with KDE density overlay, mean/median lines, and summary statistics

## Setup and Running

### Create and activate the conda environment:

```bash
conda env create -f environment.yml
conda activate bch709_hw1
```

### Run the analysis:

```bash
python analyze_mrna_gc.py
```

## Expected Outputs

- **results/mrna_metrics.tsv**
  - Tab-separated values file
  - Columns: accession, length, gc_content
  - gc_content formatted to 4 decimal places
  - Sorted by gc_content in descending order (highest to lowest)

- **results/gc_content_distribution.png**
  - Histogram plot (1600×900 px at 200 dpi)
  - Shows distribution of GC content across all mRNAs
  - Overlay: KDE density curve (blue)
  - Vertical dashed lines for mean (red) and median (green)
  - Statistics box: n, mean, median, standard deviation

## Interpretation

### Why Yeast mRNAs Show GC-Content Variation:

**1. Codon Usage Bias and Expression Level**
Yeast genes can use different codons depending on how often they are expressed. Highly expressed genes often use preferred codons that contain more G and C bases. This creates distinct subpopulations of mRNAs with different overall GC fractions, as different genes are optimized for different expression levels. The codon preferences reflect abundant tRNAs in the cell and provide selective advantages for translation efficiency.

**2. Genomic Regional Variation and Gene Function**
Different parts of the genome can experience different mutation patterns over time. These different genomic regions, variation in gene function, and untranslated regions (UTRs) can lead to changes in GC content between mRNAs. Functional gene classes may have different GC content preferences depending on their role in the cell, and regional mutation/repair biases create apparent subpopulations in the observed GC-content distribution.

