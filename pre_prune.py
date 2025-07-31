import os 
import glob


from pre_processing import  adama_process_summary_stats_before_pruning
from utils import load_paths

#path for all gwas files 
main_folder = "/Users/adama/Gladstone Dropbox/Adama Saccoh/March_18_ukbrap_wynton"

#folder names for each group of traits 
traits = ["aortic_area", "atrial_volume","ventrical_volume" ,"wall_thickness"]

#load each individual path with function load paths 
files_use=load_paths(traits,main_folder)


# combine_all_traits_memory_efficient(files_use[1],
#     logp_thresh=7,
#     cols_keep=['ID', 'LOG10P'],
#     output_folder="/Users/adama/Downloads/Tchandjieu_Lab/NMF_python/pre_pruned_files"

# )

df_dict, all_traits_df=adama_process_summary_stats_before_pruning(files_use[1],
    logp_thresh=7,
    pval_thresh=5e-8,  # still unused, kept for future extension
    cols_keep=['ID', 'LOG10P','CHROM', 'GENPOS','ALLELE0','ALLELE1','A1FREQ','BETA','SE'],
    output_folder="/Users/adama/Downloads/Tchandjieu_Lab/Adama_bnmf_handover/files/gwas_sumstat"

)
