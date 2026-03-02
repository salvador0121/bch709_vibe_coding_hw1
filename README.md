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

**1. Codon Usage and Translational Selection**
Highly expressed yeast genes tend to use optimal codons that are abundant in their cellular tRNA pools. Since codon bias often correlates with GC content (certain amino acids have GC-rich or AT-rich codon families), this creates distinct subpopulations: highly expressed, codon-optimized genes may cluster at higher or lower GC content depending on the amino acid composition and the organism's tRNA pool, leading to observed heterogeneity in the overall mRNA population.

**2. Functional Gene Compartmentalization and Isochores**
Different functional gene classes (e.g., ribosomal proteins vs. metabolic enzymes, constitutively expressed vs. stress-responsive) may occupy distinct genomic "isochores" or have different mutational/repair biases. Additionally, untranslated regions (UTRs) can have systematically different GC content than coding sequences. These functional and structural differences can produce apparent GC-content subpopulations, where genes fall into clusters reflecting their biological role and genomic context rather than a continuous uniform distribution.

