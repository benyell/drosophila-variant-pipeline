import pandas as pd
import io

def generate_final_report(vcf_path, output_csv):
    print(f"Extracting high-quality biological insights from {vcf_path}...")
    
    with open(vcf_path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    
    df = pd.read_csv(io.StringIO(''.join(lines)), sep='\t')
    
    # Select key columns for the researcher
    summary = df[['#CHROM', 'POS', 'REF', 'ALT', 'QUAL']].copy()
    
    # Extract Depth (DP) and Type from the INFO field
    summary['DEPTH'] = df['INFO'].str.extract(r'DP=(\d+)').astype(int)
    summary['TYPE'] = df['INFO'].str.extract(r'TYPE=([^;]+)')
    
    # Filter for high-confidence only (in case not already filtered)
    final_hits = summary[summary['QUAL'] > 30].sort_values(by='QUAL', ascending=False)
    
    final_hits.to_csv(output_csv, index=False)
    print(f"Final Research Table saved to {output_csv}")
    print("\n--- TOP BIOLOGICAL TARGETS ---")
    print(final_hits.head(10))

if __name__ == "__main__":
    generate_final_report("results/vcf/sample_filtered.vcf", "final_variant_summary.csv")