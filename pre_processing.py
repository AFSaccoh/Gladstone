import os 
from glob import glob
import gzip
import numpy as np
import matplotlib as mpl
import pandas as pd 
import seaborn as sns 
from functools import reduce





def adama_process_summary_stats_before_pruning(
    files_by_trait,
    cols_keep,
    logp_thresh=7.3,
    pval_thresh=5e-8,  # unused for now
    output_folder="/Users/adama/Downloads/Tchandjieu_Lab/NMF_python/pre_pruned_files"
):
    os.makedirs(output_folder, exist_ok=True)

    df_dict = {}
    combined_table = []

    for trait, file_list in files_by_trait.items():
        filtered_trait_dfs = []

        for filepath in file_list:
            try:
                if not filepath.endswith('.regenie') and not filepath.endswith('.gz'):
                    print(f"Skipping invalid file: {filepath}")
                    continue

                df = pd.read_csv(filepath, delim_whitespace=True)
                df=df[cols_keep]
                df = df[df['ID'].str.startswith('rs')]

                df_filtered = df[df['LOG10P'] >= logp_thresh].copy()
                df_filtered['p_value'] = 10 ** (-df_filtered['LOG10P'])
                df_filtered['z_score']=df_filtered['BETA']/df['SE']
                if df_filtered.empty:
                    continue

                # Get label from filename
                base = os.path.basename(filepath).replace('.regenie.gz', '').replace('.regenie', '')
                label = '_'.join(base.split('_')[-3:])

                # Filter and format
                df_filtered['trait'] = label
                df_filtered = df_filtered[cols_keep + ['trait','p_value']]
                #df_filtered=df_filtered['ID', 'trait','p_value','CHROM', 'GENPOS']
                filtered_trait_dfs.append(df_filtered)

            except Exception as e:
                print(f"Error processing {filepath}: {e}")

        if filtered_trait_dfs:
            combined_df = pd.concat(filtered_trait_dfs, ignore_index=True)

            # Save individual trait file
            output_path = os.path.join(output_folder, f"{trait}_prepruned.csv")
            combined_df.to_csv(output_path, index=False)
            print(f"✅ Saved: {output_path}")

            df_dict[trait] = combined_df
            combined_table.append(combined_df)

    # Save combined file
    if combined_table:
        all_traits_df = pd.concat(combined_table, ignore_index=True)
        all_traits_df.to_csv(os.path.join(output_folder, "all_cardiac_traits_prepruned.csv"), index=False)
        print("✅ All traits combined and saved.")
    else:
        all_traits_df = pd.DataFrame()
        print("⚠️ No traits passed filtering.")

    return df_dict, all_traits_df



def all_combine_nofilter(
    files_by_trait,
    cols_keep,
    output_folder="/Users/adama/Downloads/Tchandjieu_Lab/NMF_python/raw_data"
):
    os.makedirs(output_folder, exist_ok=True)

    df_dict = {}
    combined_table = []

    for trait, file_list in files_by_trait.items():
        filtered_trait_dfs = []

        for filepath in file_list:
            try:
                if not filepath.endswith('.regenie') and not filepath.endswith('.gz'):
                    print(f"Skipping invalid file: {filepath}")
                    continue

                df = pd.read_csv(filepath, delim_whitespace=True)
                df=df[cols_keep]
                df = df[df['ID'].str.startswith('rs')]
            
                # Get label from filename
                base = os.path.basename(filepath).replace('.regenie.gz', '').replace('.regenie', '')
                label = '_'.join(base.split('_')[-3:])

                # Filter and format
                df_filtered['trait'] = label
                df_filtered = df_filtered[cols_keep + ['trait','p_value']]
                #df_filtered=df_filtered['ID', 'trait','p_value','CHROM', 'GENPOS']
                filtered_trait_dfs.append(df_filtered)
              

            except Exception as e:
                print(f"Error processing {filepath}: {e}")

        if filtered_trait_dfs:
            combined_df = pd.concat(filtered_trait_dfs, ignore_index=True)

            # Save individual trait file
            output_path = os.path.join(output_folder, f"{trait}_prepruned.csv")
            combined_df.to_csv(output_path, index=False)
            print(f"✅ Saved: {output_path}")

            df_dict[trait] = combined_df
            combined_table.append(combined_df)

    # Save combined file
    if combined_table:
        all_traits_df = pd.concat(combined_table, ignore_index=True)
        all_traits_df.to_csv(os.path.join(output_folder, "all_cardiac_traits_raw.csv"), index=False)
        print("✅ All traits combined and saved.")
    else:
        all_traits_df = pd.DataFrame()
        print("⚠️ No traits passed filtering.")

    return df_dict, all_traits_df




