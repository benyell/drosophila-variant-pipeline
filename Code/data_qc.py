import os
import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import numpy as np

# Set paths
REF_PATH = "data/ref/genome.fa"
R1_PATH = "data/raw/sample_R1.fastq"
R2_PATH = "data/raw/sample_R2.fastq"
OUTPUT_DIR = "results/plots"

os.makedirs(OUTPUT_DIR, exist_ok=True)

def analyze_reference(path):
    print(f"Analyzing Reference: {path}")
    record = next(SeqIO.parse(path, "fasta"))
    seq_len = len(record.seq)
    gc = gc_fraction(record.seq) * 100
    
    # Count bases
    counts = {base: record.seq.count(base) for base in ["A", "C", "G", "T", "N"]}
    
    return seq_len, gc, counts

def analyze_fastq(path, limit=10000):
    print(f"Analyzing FASTQ: {path}")
    gc_contents = []
    qualities = []
    
    for i, record in enumerate(SeqIO.parse(path, "fastq")):
        if i >= limit: break
        gc_contents.append(gc_fraction(record.seq) * 100)
        qualities.append(record.letter_annotations["phred_quality"])
    
    return gc_contents, np.array(qualities)

# 1. Process Reference
ref_len, ref_gc, ref_counts = analyze_reference(REF_PATH)

# 2. Process Reads
r1_gc, r1_qual = analyze_fastq(R1_PATH)
r2_gc, r2_qual = analyze_fastq(R2_PATH)

# --- VISUALIZATION ---

fig, axs = plt.subplots(2, 2, figsize=(15, 12))

# Plot 1: Reference Base Distribution
axs[0, 0].bar(ref_counts.keys(), ref_counts.values(), color=['#ff9999','#66b3ff','#99ff99','#ffcc99','#c2c2f0'])
axs[0, 0].set_title(f"Reference (Chr4) Base Composition\nTotal Length: {ref_len:,} bp")
axs[0, 0].set_ylabel("Count")

# Plot 2: GC Content Comparison
axs[0, 1].hist(r1_gc, bins=50, alpha=0.5, label='Read 1 (Forward)', color='blue')
axs[0, 1].hist(r2_gc, bins=50, alpha=0.5, label='Read 2 (Reverse)', color='green')
axs[0, 1].axvline(ref_gc, color='red', linestyle='dashed', linewidth=2, label=f'Ref GC ({ref_gc:.1f}%)')
axs[0, 1].set_title("GC Content Distribution")
axs[0, 1].set_xlabel("% GC")
axs[0, 1].legend()

# Plot 3: Quality Score Profile (R1)
avg_qual_r1 = np.mean(r1_qual, axis=0)
axs[1, 0].plot(avg_qual_r1, color='blue')
axs[1, 0].axhline(30, color='red', linestyle='--', label='Q30 (99.9% Accuracy)')
axs[1, 0].set_title("Mean Quality Scores per Position (R1)")
axs[1, 0].set_xlabel("Position in Read (bp)")
axs[1, 0].set_ylabel("Phred Score")
axs[1, 0].set_ylim(0, 42)
axs[1, 0].legend()

# Plot 4: Quality Score Profile (R2)
avg_qual_r2 = np.mean(r2_qual, axis=0)
axs[1, 1].plot(avg_qual_r2, color='green')
axs[1, 1].axhline(30, color='red', linestyle='--', label='Q30')
axs[1, 1].set_title("Mean Quality Scores per Position (R2)")
axs[1, 1].set_xlabel("Position in Read (bp)")
axs[1, 1].set_ylim(0, 42)

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/data_quality_report3.png")
print(f"\nAnalysis complete! Plot saved to {OUTPUT_DIR}/data_quality_report3.png")