"""
Write a complete Python script named analyze_mrna_gc.py for Homework 1.

The goal is to analyze GC content in yeast mRNA sequences.

The input file is:
data/mrna.fa.gz

It is a gzipped FASTA file. Each sequence looks like this:

>BC001547 /gb=BC001547 /gi=12654078 /ug=Sc.3456 /len=1254
ATGTCTGCTCCAGCTAGCAGTGAAACTTTATTCAGAAACTGCTTAG...

Instructions:

1) Read the gzipped FASTA file using streaming (do not load the entire file into memory).
2) For each FASTA record:
   - The accession is the first word after ">" in the header (for example, BC001547).
   - Combine all sequence lines.
   - Count only the letters A, C, G, T, and N when calculating sequence length.
   - Compute GC content as (number of G + number of C) divided by length.
   - If length is zero, set GC content to 0.

3) Create a directory named "results" if it does not exist.

4) Write a file named:
   results/mrna_metrics.tsv

   The file must:
   - Have the header:
     accession<TAB>length<TAB>gc_content
   - Format gc_content to exactly 4 decimal places.
   - Be sorted from highest GC content to lowest.

5) Create a plot saved as:
   results/gc_content_distribution.png

   The plot must:
   - Use matplotlib only.
   - Be exactly 1600×900 pixels at 200 dpi.
   - Show a histogram of GC content (x-axis from 0 to 1).
   - Normalize the histogram to density.
   - Overlay a smooth density curve (KDE).
   - Show vertical dashed lines for mean and median.
   - Display a small text box showing:
     n, mean, median, and standard deviation.

6) If the file data/mrna.fa.gz does not exist, print a clear error message and exit.

7) Print a summary to the terminal after running (n, mean, median, sd, and output file names).

At the very top of the script, include this entire prompt inside a triple-quoted comment block so it can be submitted.

At the bottom of the script, include a short interpretation section explaining two possible biological reasons why yeast mRNA sequences might show different GC-content groups.

Use only:
- Python standard library
- numpy
- matplotlib
- scipy (optional for KDE)

Return only the full Python script.
"""

import os
import sys
import gzip
import numpy as np
import matplotlib.pyplot as plt

# try to import KDE from scipy if available
try:
    from scipy.stats import gaussian_kde
    _has_scipy = True
except ImportError:
    _has_scipy = False


def compute_metrics(fasta_path):
    """
    Stream-parse a gzipped FASTA file and compute (accession, length, gc_content)
    for each record.
    """
    metrics = []
    if not os.path.exists(fasta_path):
        sys.stderr.write(f"ERROR: input file '{fasta_path}' not found\n")
        sys.exit(1)

    with gzip.open(fasta_path, "rt") as fh:
        accession = None
        seq_chunks = []
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                # process previous record
                if accession is not None:
                    seq = "".join(seq_chunks).upper()
                    valid = [c for c in seq if c in "ACGTN"]
                    length = len(valid)
                    gc = (valid.count("G") + valid.count("C")) / length if length > 0 else 0.0
                    metrics.append((accession, length, gc))
                header = line[1:].strip()
                accession = header.split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        # last record
        if accession is not None:
            seq = "".join(seq_chunks).upper()
            valid = [c for c in seq if c in "ACGTN"]
            length = len(valid)
            gc = (valid.count("G") + valid.count("C")) / length if length > 0 else 0.0
            metrics.append((accession, length, gc))

    return metrics


def write_tsv(metrics, out_path):
    """
    Write metrics list to TSV sorted by gc_content descending.
    """
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    sorted_metrics = sorted(metrics, key=lambda x: x[2], reverse=True)
    with open(out_path, "w") as out:
        out.write("accession\tlength\tgc_content\n")
        for acc, length, gc in sorted_metrics:
            out.write(f"{acc}\t{length}\t{gc:.4f}\n")


def make_kde(xs):
    """
    return a callable KDE estimating density over xs.
    If scipy is available use gaussian_kde; otherwise,
    implement simple Gaussian KDE using numpy.
    """
    if _has_scipy:
        kde = gaussian_kde(xs)
        return lambda x: kde(x)
    else:
        # simple KDE with Silverman's rule for bandwidth
        n = len(xs)
        std = np.std(xs, ddof=1)
        h = 1.06 * std * n ** (-0.2) if std > 0 else 0.1
        def kde_func(x):
            X = np.atleast_1d(x)
            dens = np.zeros_like(X, dtype=float)
            coef = 1 / (n * h * np.sqrt(2 * np.pi))
            for xi in xs:
                dens += np.exp(-0.5 * ((X - xi) / h) ** 2)
            return coef * dens
        return kde_func


def plot_distribution(gc_values, out_png):
    """
    Create histogram+KDE plot of GC content values.
    """
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    n = len(gc_values)
    mean = np.mean(gc_values)
    median = np.median(gc_values)
    sd = np.std(gc_values, ddof=1)

    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=200)
    # histogram
    ax.hist(gc_values, bins=50, range=(0, 1), density=True,
            color="lightgray", edgecolor="black")
    # KDE
    kde = make_kde(gc_values)
    xs = np.linspace(0, 1, 1000)
    ax.plot(xs, kde(xs), color="blue", lw=1.5)

    # vertical lines
    ax.axvline(mean, color="red", linestyle="--", label="mean")
    ax.axvline(median, color="green", linestyle="--", label="median")

    # caption box
    text = f"n={n}\nmean={mean:.4f}\nmedian={median:.4f}\nsd={sd:.4f}"
    ax.text(0.95, 0.95, text, transform=ax.transAxes,
            ha="right", va="top", bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"))
    ax.set_xlabel("GC content")
    ax.set_ylabel("Density")
    ax.set_xlim(0, 1)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)
    return n, mean, median, sd


if __name__ == "__main__":
    fasta = "data/mrna.fa.gz"
    metrics = compute_metrics(fasta)
    if not metrics:
        sys.stderr.write("no records found in FASTA\n")
        sys.exit(1)
    tsv_path = "results/mrna_metrics.tsv"
    png_path = "results/gc_content_distribution.png"
    write_tsv(metrics, tsv_path)
    gc_list = np.array([m[2] for m in metrics])
    n, mean, median, sd = plot_distribution(gc_list, png_path)
    # print summary
    sys.stdout.write(
        f"n={n} mean={mean:.4f} median={median:.4f} sd={sd:.4f}\n"
        f"metrics: {tsv_path}\nplot: {png_path}\n"
    )

# Interpretation: Two biological reasons for the GC content distribution pattern in yeast mRNA
#
# 1) Codon usage bias and expression level:
#    Yeast genes can use different codons depending how often they are expressed.
#    Highly expressed genes often use preferred codons that contain more G and C bases.
#    This creates distinct subpopulations of mRNAs with different overall GC fractions,
#    as different genes are optimized for different expression levels.
#
# 2) Genomic regional variation and gene function:
#    Different parts of the genome can experience different mutation patterns over time.
#    These different parts, variation in gene function, and untranslated regions (UTRs)
#    can lead to changes in GC content between mRNAs. Functional gene classes may occupy
#    distinct genomic regions with their own mutation/repair biases, producing a
#    heterogeneous GC distribution.
