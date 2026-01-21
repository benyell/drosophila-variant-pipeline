import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pysam
import io
import os

# --- PATHS ---
VCF_PATH = "results/vcf/sample_filtered.vcf"
BAM_PATH = "results/vcf/sample.bam"
OUTPUT_DIR = "results/insights"
os.makedirs(OUTPUT_DIR, exist_ok=True)

def parse_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(io.StringIO(''.join(lines)), sep='\t')

print("Starting Comprehensive Analysis...")

# 1. Load Data
df = parse_vcf(VCF_PATH)
bam = pysam.AlignmentFile(BAM_PATH, "rb")

# 2. Extract Biological Insights: Ti/Tv Ratio
# Transitions (A<->G, C<->T) vs Transversions (others)
transitions = ['AG', 'GA', 'CT', 'TC']
transversions = ['AC', 'CA', 'AT', 'TA', 'GC', 'CG', 'GT', 'TG']

def get_titv(ref, alt):
    pair = f"{ref}{alt}"
    if pair in transitions: return "Transition"
    if pair in transversions: return "Transversion"
    return "Other"

df['mutation_type'] = df.apply(lambda x: get_titv(x['REF'], x['ALT']), axis=1)
titv_counts = df['mutation_type'].value_counts()
titv_ratio = titv_counts.get('Transition', 0) / max(1, titv_counts.get('Transversion', 0))

# 3. Extract Alignment Insights: Depth of Coverage
depths = [col.nsegments for col in bam.pileup()]
avg_depth = sum(depths) / len(depths) if depths else 0

# --- REPORT GENERATION ---

print(f"\n--- RESEARCH SUMMARY ---")
print(f"Total Variants Identified: {len(df)}")
print(f"Mean Read Depth: {avg_depth:.2f}x")
print(f"Ti/Tv Ratio: {titv_ratio:.2f}")

# Plotting
fig, axes = plt.subplots(2, 2, figsize=(16, 12))

# Plot 1: Variant Types (SNP vs Complex)
df['TYPE'] = df['INFO'].str.extract(r'TYPE=([^;]+)')
sns.countplot(data=df, x='TYPE', ax=axes[0,0], palette='viridis')
axes[0,0].set_title("Distribution of Variant Types")

# Plot 2: Ti/Tv Distribution
titv_counts.plot(kind='pie', autopct='%1.1f%%', ax=axes[0,1], colors=['#ff9999','#66b3ff'])
axes[0,1].set_title(f"Ti/Tv Ratio: {titv_ratio:.2f}")

# Plot 3: Quality Scores (QUAL)
sns.histplot(df['QUAL'], kde=True, ax=axes[1,0], color='purple')
axes[1,0].set_title("Variant Confidence (QUAL) Distribution")

# Plot 4: Position Density Across Chromosome
sns.scatterplot(data=df, x='POS', y='QUAL', alpha=0.5, ax=axes[1,1], color='orange')
axes[1,1].set_title("Variant Quality by Genomic Position")

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/genomic_insights_report3.png")
print(f"Full Report Saved to: {OUTPUT_DIR}/genomic_insights_report3.png")