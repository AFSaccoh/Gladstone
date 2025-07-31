
import pandas as pd 
import os 


def extract_snp_and_chromosome(folder_path, sep="\t"):
    all_entries = []

    for file in os.listdir(folder_path):
        if file.endswith(".clumps"):
            file_path = os.path.join(folder_path, file)
            try:
                df = pd.read_csv(file_path, sep=sep)
                

                if "ID" in df.columns and "#CHROM" in df.columns:
                    df = df[["ID", "#CHROM"]].copy()
                    df.rename(columns={"#CHROM": "chromosome"}, inplace=True)
                    all_entries.append(df)
                else:
                    print(f"Skipping {file}: Required columns not found.")
            except Exception as e:
                print(f"Error loading {file}: {e}")

    return pd.concat(all_entries, ignore_index=True)

folder_path="/Users/adama/Downloads/Tchandjieu_Lab/adama_snps_all"
pruned_df=extract_snp_and_chromosome(folder_path, sep="\t")


pruned_df['ID'] = pruned_df['ID'].astype(str)

# Large file
input_file = "/Users/adama/Downloads/Tchandjieu_Lab/NMF_python/raw_data/all_cardiac_traits_raw.csv"  # or .zip
output_file = "pre_bnmf_merged_output.csv"

chunksize = 100_000
first_write = True

for chunk in pd.read_csv(input_file, chunksize=chunksize):
    # Make sure required columns are there
    required_cols = ["CHROM", "GENPOS", "ALLELE1", "ALLELE0"]
    if not all(col in chunk.columns for col in required_cols):
        continue

    # Convert to string and create full_id_10
    chunk['full_id_10'] = (
        chunk["CHROM"].astype(str) + ":" +
        chunk["GENPOS"].astype(str) + ":" +
        chunk["ALLELE1"].astype(str) + ":" +
        chunk["ALLELE0"].astype(str)
    )

    # Merge with pruned table
    merged = pd.merge(chunk, pruned_df, left_on='full_id_10', right_on='ID', how='inner')

    # Save only if there's something to save
    if not merged.empty:
        merged.to_csv(output_file, mode='a', index=False, header=first_write)
        first_write = False

print("✅ Merging complete. Output saved to:", output_file)