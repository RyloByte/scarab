import glob
from os.path import join as joinpath

import pandas as pd

# specify that all columns should be shown
pd.set_option('max_columns', None)

# Input file

# Files and dirs
scarabout_path = '/home/ryan/SCARAB_local/benchmarking_output/' \
                'single_binner_bench/MGE_6/SCARAB_single/0/' \
                '*/*/'
denovo_file_list = glob.glob(joinpath(scarabout_path, '*.denovo_clusters.tsv'))
trusted_file_list = glob.glob(joinpath(scarabout_path, '*.hdbscan_clusters.tsv'))
ocsvm_file_list = glob.glob(joinpath(scarabout_path, '*.ocsvm_clusters.tsv'))
inter_file_list = glob.glob(joinpath(scarabout_path, '*.inter_clusters.tsv'))

scarab_single_file = '/home/ryan/SCARAB_local/benchmarking_output/errstat_inputs/SCARAB.single.errstat.tsv'
unitem_single_file = '/home/ryan/SCARAB_local/benchmarking_output/errstat_inputs/UniteM.single.errstat.tsv'
vamb_multi_file = '/home/ryan/SCARAB_local/benchmarking_output/errstat_inputs/VAMB.multi.errstat.tsv'

# Load stats tables
scarab_single_df = pd.read_csv(scarab_single_file, header=0, sep='\t')
scarab_single_df['binner'] = ['_'.join(['SCARAB', str(x), str(y), str(z)])
                             for x, y, z in
                             zip(scarab_single_df['algorithm'],
                                 scarab_single_df['mode'],
                                 scarab_single_df['param_set']
                                 )
                             ]
scarab_mge_df = scarab_single_df.query("sample_type == 'MGE_6'")
scarab_mge_df.rename(columns={'>20Kb': 'over20Kb'}, inplace=True)
scarab_mg_df = scarab_mge_df.query("over20Kb == 'Yes' & "
                                 "MQ_bins == 'Yes' &"
                                 "level == 'strain_absolute'"
                                 )
print(scarab_mg_df.head())
flurp
unitem_single_df = pd.read_csv(unitem_single_file, header=0, sep='\t')
unitem_mge_df = unitem_single_df.query("sample_type == 'MGE_6'")

vamb_multi_df = pd.read_csv(vamb_multi_file, header=0, sep='\t')
vamb_multi_df['binner'] = 'VAMB'
vamb_mge_df = vamb_multi_df.query("sample_type == 'MGE_6'")

bin_cat_df = pd.concat([scarab_single_df, unitem_single_df,
                        vamb_multi_df
                        ])
mge_df = bin_cat_df.query("sample_type == 'MGE_6'")
print(mge_df.head())
