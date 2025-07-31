
import os
import pandas as pd
from utils import load_paths

main_folder = "/Users/adama/Gladstone Dropbox/Adama Saccoh/March_18_ukbrap_wynton"
traits = ["aortic_area", "atrial_volume","ventrical_volume" ,"wall_thickness"]

files_use=load_paths(traits,main_folder)[1]
print(files_use)


df_dict = {}  
cols_keep=['CHROM', 'GENPOS', 'ID', 'ALLELE0', 'ALLELE1', 'A1FREQ','LOG10P','z_score', 'p_value',
       'label', 'N']
zscore_tables_unfiltered=[]
zscore_tables_filtered=[]

pval_thresh=5e-8
logp_thresh=7.3

rename_dict={'CHROM':'CHR',
             'ALLELE0':'REF' ,
             'ALLELE1':'ALT'

}

for i in files_use["atrial_volume"]:

    try:
        if not i.endswith('.regenie') and not i.endswith('.gz'):
            print(f"skipping invalid file {i}")
            continue  # Skip to the next file

        # Example: read the file (adjust as needed)
        df = pd.read_csv(i, delim_whitespace=True)
        df_new=df.copy()
        df_new=df[df['ID'].str.startswith('rs')].copy()
        df_new=df_new[df_new['ID'].str.startswith('rs')]
        df_new['z_score']=df_new['BETA']/df['SE']
        df_new['p_value']= 10**(-df['LOG10P']) 



         # Extract label from filename (e.g., "DAo_max_area")
        base = os.path.basename(i)
        base = base.replace('.regenie.gz', '').replace('.regenie', '')
        label = '_'.join(base.split('_')[-3:])  # Adjust depending on the pattern
        df_new['label'] = label

        df_dict[label] = df_new 
        df_new=df_new[cols_keep]

        #df=df.rename(columns=rename_dict,inplace=True)
        
        filtered_logp=df_new[df_new['LOG10P']>= logp_thresh]
        filtered_pval=df_new[df_new['p_value']<= pval_thresh]
        
        #new z score columns
        z_col_name = f"z_score_{label}"
        z_df = df_new[['ID', 'z_score']].rename(columns={'z_score': z_col_name})
        zscore_tables_unfiltered.append(z_df)

        z_fil_df =filtered_logp[['ID', 'z_score']].rename(columns={'z_score': z_col_name})

        zscore_tables_filtered.append(z_fil_df)

        #print(df)
        
        # ... process df ...
    except Exception as e:
        print(f"Error processing {i}: {e}")


print(zscore_tables_filtered)