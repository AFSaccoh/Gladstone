
import os 
import gc
from glob import glob
import gzip
import numpy as np
import matplotlib as mpl
import pandas as pd 
import seaborn as sns 
import matplotlib.pyplot as plt
from functools import reduce

#from checks import check_files



def load_paths(traits,main_folder):
    """ loads files from the main folder and returns a list of files for each trait 
    Args:
        traits (list): List of trait names to search for.
        main_folder (str): Path to the main folder containing trait directories.
    """
    files_by_trait = {}

    for trait in traits:

        trait_path = os.path.join(main_folder, trait)

        step_dirs = [
        d for d in os.listdir(trait_path)
        if d.startswith("step_2") and os.path.isdir(os.path.join(trait_path, d))
        ]
    
        if not step_dirs:
            print(f"No 'step_2_' folder in {trait_path}")
            files_by_trait[trait] = []
            
            continue

        step_dir = step_dirs[0]  # Use first match
        merged_path = os.path.join(trait_path, step_dir, "merged")

        if not os.path.exists(merged_path):
             print(f"No 'merged' folder in {step_dir}")
             files_by_trait[trait] = []
             continue

    # Get all files (adjust pattern as needed)
        files = glob(os.path.join(merged_path, "*")) #cant use suffix .regenie.gz as missing from aorta, just skip the first file from the aortic traits you can then filter later on 
        files_by_trait[trait] = files

# Optional: Flatten into one list - this  contains file paths for everything 
        all_files = sum(files_by_trait.values(), [])

    return all_files, files_by_trait




def process_file(dataframe, logp_thresh,output_file="zscore_matrix_wide.csv",):

    df_list = []

    for i in dataframe:
        try:
            if not i.endswith('.regenie') and not i.endswith('.gz'):
                print(f"Skipping invalid file {i}")
                continue

        # Load relevant columns
            df = pd.read_csv(i, delim_whitespace=True, usecols=['ID', 'BETA', 'SE', 'LOG10P'])
            df = df[df['ID'].str.startswith('rs')]  # Only keep rsIDs

        # Calculate z-score
            df['z_score'] = df['BETA'] / df['SE']

        # Extract label
            base = os.path.basename(i).replace('.regenie.gz', '').replace('.regenie', '')
            label = '_'.join(base.split('_')[-3:])
            z_col = f"z_score_{label}"

        # Filter SNPs
            filtered = df[df['LOG10P'] >= logp_thresh]

        # Keep only ID and z_score, rename column
            z_df = filtered[['ID', 'z_score']].rename(columns={'z_score': z_col})
            df_list.append(z_df)

            del df, filtered, z_df
            gc.collect()

        except Exception as e:
            print(f"Error processing {i}: {e}")


    zscore_matrix = reduce(lambda left, right: pd.merge(left, right, on='ID', how='outer'), df_list)

# Save full matrix
    zscore_matrix.to_csv(output_file, index=False)
    print(f"Saved matrix: shape={zscore_matrix.shape}")


def plot_trait_correlation_heatmaps(input_folder, output_folder, corr_threshold=0.85):
    os.makedirs(output_folder, exist_ok=True)

    for file in os.listdir(input_folder):
        if file.endswith(".csv"):
            file_path = os.path.join(input_folder, file)

            try:
                zscore_df = pd.read_csv(file_path)

                if 'ID' in zscore_df.columns:
                    zscore_df = zscore_df.drop(columns='ID')

                zscore_df = zscore_df.fillna(0)
                corr_matrix = zscore_df.corr(method='pearson')

                to_drop = set()
                for i in range(len(corr_matrix.columns)):
                    trait_i = corr_matrix.columns[i]
                    if trait_i in to_drop:
                        continue
                    for j in range(i + 1, len(corr_matrix.columns)):
                        trait_j = corr_matrix.columns[j]
                        if trait_j in to_drop:
                            continue
                        if corr_matrix.iloc[i, j] >= corr_threshold:
                            to_drop.add(trait_j)

                pruned_mat = corr_matrix.drop(columns=to_drop, errors='ignore')

                # Plot and save heatmap
                plt.figure(figsize=(12, 10))
                sns.heatmap(corr_matrix, annot=True, fmt=".2f", cmap="coolwarm", square=True, cbar_kws={"shrink": 0.7})
                plt.title(f"Genetic Correlation: {file}")
                plt.tight_layout()

                output_filename = os.path.splitext(file)[0] + "_heatmap.png"
                output_path = os.path.join(output_folder, output_filename)
                plt.savefig(output_path)
                plt.close()

                print(f"Saved heatmap to {output_path}")
                print(f"Pruned traits for {file}: {sorted(to_drop)}\n")

            except Exception as e:
                print(f"Failed to process {file}: {e}")
                
