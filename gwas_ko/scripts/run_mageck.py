import os
import subprocess
import pandas as pd
import sys

# ===== CONFIG =====
base_dir = '/data/JE_misc/GW_screen/'
library_file = os.path.join(base_dir, 'brunello_library.txt')
sample_info_file = os.path.join(base_dir, 'sample_info.csv')
# ==================

def generate_library_file(input_file, output_file):
    df = pd.read_csv(input_file, sep='\t')
    library_df = df[['sgRNA Target Sequence', 'Target Gene Symbol']].copy()
    library_df.columns = ['sgRNA Sequence', 'Target Gene']
    library_df.insert(0, 'sgRNA ID', [f'sgRNA_{i+1}' for i in range(len(library_df))])
    library_df.to_csv(output_file, sep='\t', index=False, header=False)

def run_mageck(donor, condition):
    df = pd.read_csv(sample_info_file)
    
    # Filter to matching donor and condition
    group_df = df[(df['Donor'].str.replace(' ', '').str.lower() == donor.strip().lower().replace(' ', '')) &
                  (df['Condition'].str.replace(' ', '').str.lower() == condition.strip().lower().replace(' ', ''))]
    
    if group_df.empty:
        print(f"No matching samples for {donor} {condition}")
        return

    # Group by replicate index — here we assume 1 vs 4, 31 vs 34, etc.
    # Group Q1 and Q4 separately
    q1 = group_df[group_df['Quartile'] == 'Q1'].sort_values(by='sample').reset_index(drop=True)
    q4 = group_df[group_df['Quartile'] == 'Q4'].sort_values(by='sample').reset_index(drop=True)

    if len(q1) != len(q4):
        print(f"Mismatched number of Q1 and Q4 samples for {donor} {condition}")
        return

    for i in range(len(q1)):
        rep_name = f"rep{i+1}"
        q1_sample = str(int(q1.loc[i, 'sample']))
        q4_sample = str(int(q4.loc[i, 'sample']))

        # Locate matching FASTQ files
        group_dir = os.path.join(base_dir, donor.strip(), condition.strip(), rep_name)
        os.makedirs(group_dir, exist_ok=True)

        fastq_base_dir = os.path.join(base_dir, donor.strip(), condition.strip())

        q1_fastq = next((f for f in os.listdir(fastq_base_dir) if f.startswith(f"{q1_sample}_S")), None)
        q4_fastq = next((f for f in os.listdir(fastq_base_dir) if f.startswith(f"{q4_sample}_S")), None)

        if not q1_fastq or not q4_fastq:
            print(f"Missing FASTQ for {donor} {condition} {rep_name}: {q1_sample} or {q4_sample}")
            continue

        q1_path = os.path.join(fastq_base_dir, q1_fastq)
        q4_path = os.path.join(fastq_base_dir, q4_fastq)

        output_prefix = f"mageck_{donor}_{condition}_{rep_name}"
        print(f"[{donor} {condition} {rep_name}] Running MAGeCK count...")

        count_cmd = [
            'mageck', 'count',
            '-l', library_file,
            '--fastq', q1_path, q4_path,
            '--sample-label', 'control,treatment',
            '-n', output_prefix
        ]
        subprocess.run(count_cmd, cwd=group_dir, check=True)

        print(f"[{donor} {condition} {rep_name}] Running MAGeCK test...")
        test_cmd = [
            'mageck', 'test',
            '-k', f"{output_prefix}.count.txt",
            '-t', 'treatment',
            '-c', 'control',
            '-n', f"{output_prefix}_test"
        ]
        subprocess.run(test_cmd, cwd=group_dir, check=True)

        print(f"[{donor} {condition} {rep_name}] Analysis complete.")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python run_mageck_single.py <Donor> <Condition>")
        sys.exit(1)

    donor_arg = sys.argv[1]
    condition_arg = sys.argv[2]

    
    run_mageck(donor_arg, condition_arg)
