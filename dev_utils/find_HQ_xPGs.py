import os

import pandas as pd

# specify that all columns should be shown
pd.set_option('max_columns', None)

cwd_dir = "/home/ryan/Desktop/HQ_stuff"
'''
# parse all HQ trnascan output files
trnascan_file_list = glob.glob(
    os.path.join(cwd_dir, 'HQ_trnascan/scratch/st-shallam-1/mcglock/'
                          'SI60/SCARAB_output/SI060_*/*/*/trnascan/*.tsv'
                 ))
trna_df_list = []
for trna_file in trnascan_file_list:
    print(trna_file)
    sample_id = trna_file.split('/')[11]
    mode = trna_file.split('/')[12]
    set = trna_file.split('/')[13]
    sag_id = trna_file.split('/')[15].split('.', 1)[0]
    filter_len = trna_file.split('/')[15].split('.')[4]
    trna_df = pd.read_csv(trna_file, sep='\t', header=None,
                          skiprows=3
                          )
    trna_df.replace({' ': ''}, inplace=True)
    trna_df.columns = ['contig_id', 'trna_cnt', 'start',
                       'end', 'type', 'anti_codon',
                       'intron_start', 'intron_end',
                       'inf_score', 'notes'
                       ]
    trna_df['sample_id'] = sample_id
    trna_df['mode'] = mode
    trna_df['set'] = set
    trna_df['sag_id'] = sag_id
    trna_df['filter_len'] = filter_len
    trna_df_list.append(trna_df)

trna_concat_df = pd.concat(trna_df_list)
trna_concat_df.to_csv(os.path.join(cwd_dir, 'trnascan_all_output.tsv'),
                      sep='\t', index=False
                      )

trna_concat_df = pd.read_csv(os.path.join(cwd_dir, 'trnascan_all_output.tsv'),
                             sep='\t', header=0
                             )


trna_trim_df = trna_concat_df[['contig_id', 'sample_id', 'mode', 'set',
                              'sag_id', 'filter_len', 'type',
                              'trna_cnt', 'notes'
                              ]]
trna_filter_df = trna_trim_df.query("notes != 'pseudo'")  # TODO: need to check if this is ok
index_list = ['sample_id', 'mode', 'set', 'sag_id', 'filter_len']
cols = 'type'
vals = 'trna_cnt'
trna_sum_df = pd.pivot_table(values=vals, index=index_list,
                             columns=cols, data=trna_filter_df,
                             aggfunc=sum
                             ).reset_index().fillna(0)

new_row_list = []
for i, row in trna_sum_df.iterrows():
    trnas = row[5:]
    row_vals = list(row.values)
    uniq_ts = sum([1 if x != 0 else 0 for x in trnas])
    sum_ts = sum(trnas)
    row_vals.append(uniq_ts)
    row_vals.append(sum_ts)
    new_row_list.append(row_vals)
new_cols = list(trna_sum_df.columns)
new_cols.append('unique_trnas')
new_cols.append('total_trnas')
trna_cnt_df = pd.DataFrame(new_row_list, columns=new_cols)
trna_cnt_df.to_csv(os.path.join(cwd_dir, 'trnascan_cnt_output.tsv'),
                   sep='\t', index=False
                   )
'''
trna_cnt_df = pd.read_csv(os.path.join(cwd_dir, 'trnascan_cnt_output.tsv'),
                          sep='\t', header=0
                          )

trna_hq_df = trna_cnt_df.query("unique_trnas >= 18")
trna_hq_df.to_csv(os.path.join(cwd_dir, 'trnascan_HQ_output.tsv'),
                  sep='\t', index=False
                  )
'''
# parse all HQ barrnap output files
barrnap_file_list = glob.glob(
    os.path.join(cwd_dir, 'HQ_barrnap/scratch/st-shallam-1/mcglock/'
                          'SI60/SCARAB_output/SI060_*/*/*/barrnap/*.gff'
                 ))
bar_df_list = []
for bar_file in barrnap_file_list:
    print(bar_file)
    try:
        sample_id = bar_file.split('/')[11]
        mode = bar_file.split('/')[12]
        set = bar_file.split('/')[13]
        sag_id = bar_file.split('/')[15].split('.', 1)[0]
        filter_len = bar_file.split('/')[15].split('.')[4]
        bar_df = pd.read_csv(bar_file, sep='\t', header=None,
                             skiprows=1
                             )
        bar_df.columns = ['contig_id', 'software_version',
                          'seq_type', 'start', 'end',
                          'e-value', 'unknown1',
                          'unknown2', 'description'
                          ]
        bar_df['sample_id'] = sample_id
        bar_df['mode'] = mode
        bar_df['set'] = set
        bar_df['sag_id'] = sag_id
        bar_df['filter_len'] = filter_len
        bar_df['subunit'] = [x.split(';', 1)[0].split('=', 1)[1].split('_', 1)[0]
                             for x in bar_df['description']
                             ]
        bar_df['product'] = [x.split(';', 1)[1].split('=', 1)[1]
                             for x in bar_df['description']
                             ]
        bar_df['notes'] = [x.split(';note=', 1)[1] if 'note' in x
                           else '' for x in bar_df['description']
                           ]
        bar_df['completeness'] = ['full' if x is '' else 'partial'
                                  for x in bar_df['notes']
                                  ]
        bar_df['percent'] = [100 if x is '' else
                             int(x.split('only ', 1
                                         )[1].split(' percent', 1)[0]
                                 )
                             for x in bar_df['notes']
                             ]
        bar_df_list.append(bar_df)
    except:
        print('Somthing wrong with gff file, probably empty...')

bar_concat_df = pd.concat(bar_df_list)
bar_concat_df.to_csv(os.path.join(cwd_dir, 'barrnap_all_output.tsv'),
                     sep='\t', index=False
                     )
'''
bar_concat_df = pd.read_csv(os.path.join(cwd_dir, 'barrnap_all_output.tsv'),
                            sep='\t', header=0
                            )
bar_trim_df = bar_concat_df[['contig_id', 'sample_id', 'mode', 'set',
                             'sag_id', 'filter_len', 'subunit',
                             'completeness', 'percent'
                             ]]
index_list = ['sample_id', 'mode', 'set', 'sag_id', 'filter_len']
cols = 'subunit'
vals = 'percent'
bar_piv_df = pd.pivot_table(values=vals, index=index_list,
                            columns=cols, data=bar_trim_df,
                            aggfunc=sum
                            ).reset_index().fillna(0)
bar_piv_df.rename(columns={"16S": "rRNA_16S",
                           "23S": "rRNA_23S",
                           "5S": "rRNA_5S"
                           }, inplace=True)

bar_hq_df = bar_piv_df.query("rRNA_16S >= 100 &"
                             "rRNA_23S >= 100 &"
                             "rRNA_5S >= 100"
                             )
bar_hq_df.to_csv(os.path.join(cwd_dir, 'barrnap_HQ_output.tsv'),
                 sep='\t', index=False
                 )

# Merge them all together
bar_trna_df = pd.merge(bar_hq_df, trna_hq_df,
                       on=index_list,
                       how='inner'
                       )
bar_trna_df.to_csv(os.path.join(cwd_dir, 'barrnap_trnascan_HQ_output.tsv'),
                   sep='\t', index=False
                   )
print("TRNASCAN: ", len(trna_hq_df['sag_id'].unique()))
print("BARRNAP: ", len(bar_hq_df['sag_id'].unique()))
print("BOTH: ", len(bar_trna_df['sag_id'].unique()))

# Now Merge with the CheckM stats
checkm_file = os.path.join(cwd_dir, 'HQ_checkm.tsv')
checkm_df = pd.read_csv(checkm_file, sep='\t', header=0)

checkm_df['sample_id'] = [x + '_' + y for x, y in
                          zip(checkm_df['sample'],
                              checkm_df['depth']
                              )]
checkm_df.rename(columns={'SAG_ID': 'sag_id'}, inplace=True)
check_bar_trna_df = pd.merge(bar_trna_df, checkm_df,
                             on=index_list,
                             how='left'
                             )
check_bar_trna_df.to_csv(os.path.join(cwd_dir, 'barrnap_trnascan_checkm_HQ_output.tsv'),
                         sep='\t', index=False
                         )

# Lastly, add the assembly stats from seqkit
sk_file = os.path.join(cwd_dir, 'HQ_seqkit_stats.tsv')
sk_df = pd.read_csv(sk_file, sep='\t', header=0)

check_bar_trna_sk_df = pd.merge(check_bar_trna_df, sk_df,
                                on=index_list,
                                how='left'
                                )
check_bar_trna_sk_df.to_csv(os.path.join(cwd_dir, 'Final_HQ_output.tsv'),
                            sep='\t', index=False
                            )
