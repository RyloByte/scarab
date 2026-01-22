import glob
import os

import pandas as pd

# import plotly.graph_objects as go
pd.set_option('display.max_columns', None)


def build_cv_tables(working_dir, sample_list, config_list):
    cv_concat_list = []
    cv_all_list = []
    for sample in sample_list:
        cv_df = pd.read_csv(sample, sep='\t', header=0)
        if len(cv_df.columns) == 16:
            col_list = list(cv_df.columns)
            col_list[-1] = 'sample_id'
            cv_df.columns = col_list
            sample_id = sample.rsplit('/')[-2].split('_', 1)[1]
            cv_df['sample_id'] = [sample_id + '_' + str(x) for x in cv_df['sample_id']]
        else:
            sample_id = sample.rsplit('/')[-2].split('_', 1)[1]
            cv_df['sample_id'] = sample_id
        # get best performing configs
        for config in config_list:
            cv_algo = config[0]
            algo = config[1]
            config_df = cv_df.query("cv_algo == @cv_algo & algo == @algo")
            for level in config_df['level'].unique():
                level_df = config_df.query("level == @level")
                nc_max = level_df['nc_cnt'].max()
                mq_max = level_df['mq_cnt'].max()
                nc_df = level_df.query("nc_cnt == @nc_max")
                mq_df = level_df.query("mq_cnt == @mq_max")
                nc_mq_df = nc_df.query("mq_cnt == @mq_max")
                mq_nc_df = mq_df.query("nc_cnt == @nc_max")
                if nc_mq_df.shape[0] != 0:
                    cv_concat_list.append(nc_mq_df)
                else:
                    cv_concat_list.append(nc_df)

                if mq_nc_df.shape[0] != 0:
                    cv_concat_list.append(mq_nc_df)
                else:
                    cv_concat_list.append(mq_df)
        cv_all_list.append(cv_df)
    cv_concat_df = pd.concat(cv_concat_list)
    cv_best_file = os.path.join(working_dir, "CV_best_table.tsv")
    cv_concat_df.to_csv(cv_best_file, sep='\t', index=False)
    cv_all_df = pd.concat(cv_all_list)
    cv_all_file = os.path.join(working_dir, "CV_all_table.tsv")
    cv_all_df.to_csv(cv_all_file, sep='\t', index=False)

    return cv_concat_df, cv_all_df


def merge_cluster_table(working_dir, cv_concat_df, cv_all_df):
    cluster_file = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                'renyi_entropy/references/cluster_clean.tsv'
                                )
    cluster_df = pd.read_csv(cluster_file, sep='\t', header=0)
    cluster_sub_df = cluster_df[['sample_id', 'sample_type', 'best_cluster', 'best_match'
                                 ]].drop_duplicates()
    clust_cv_df = cv_concat_df.merge(cluster_sub_df[['sample_id', 'sample_type',
                                                     'best_cluster', 'best_match']],
                                     on='sample_id', how='left'
                                     )
    clust_cv_df.drop(['mq_avg_p', 'mq_avg_r', 'mq_avg_mcc',
                      'nc_avg_p', 'nc_avg_r', 'nc_avg_mcc'
                      ], axis=1, inplace=True
                     )
    clust_cv_df.drop_duplicates(inplace=True)
    clust_all_df = cv_all_df.merge(cluster_sub_df[['sample_id', 'sample_type',
                                                   'best_cluster', 'best_match']],
                                   on='sample_id', how='left'
                                   )
    clust_all_df.drop(['mq_avg_p', 'mq_avg_r', 'mq_avg_mcc',
                       'nc_avg_p', 'nc_avg_r', 'nc_avg_mcc'
                       ], axis=1, inplace=True
                      )
    clust_all_df.drop_duplicates(inplace=True)
    clust_cv_df['majority_rule'] = 'MR'
    clust_all_df['majority_rule'] = 'MR'
    clust_all_file = os.path.join(working_dir, "CV_clust_table.tsv")
    clust_all_df.to_csv(clust_all_file, sep='\t', index=False)

    return clust_cv_df, clust_all_df


def aggregate_best_params(working_dir, config_list, clust_all_df):
    group_list = ['sample_type', 'best_cluster', 'sample_id', 'best_match', 'majority_rule']
    nc_agg_list = []
    mq_agg_list = []
    for config in config_list:
        cv_algo = config[0]
        algo = config[1]
        config_df = clust_all_df.query("cv_algo == @cv_algo & algo == @algo")
        for level in config_df['level'].unique():
            level_df = config_df.query("level == @level")
            for group in group_list:
                group_vals = level_df[group].unique()
                group_params = []
                for val in group_vals:
                    type_df = level_df.loc[level_df[group] == val]
                    nc_max_df = type_df.groupby(['sample_id', 'cv_param1', 'cv_param2', 'cv_val1', 'cv_val2']
                                                )['nc_cnt'].max().reset_index()
                    mq_max_df = type_df.groupby(['sample_id', 'cv_param1', 'cv_param2', 'cv_val1', 'cv_val2']
                                                )['mq_cnt'].max().reset_index()
                    nc_sort_df = nc_max_df.sort_values(['nc_cnt'], ascending=[False])
                    mq_sort_df = mq_max_df.sort_values(['mq_cnt'], ascending=[False])

                    nc_dup_df = nc_sort_df.drop_duplicates(subset=['sample_id'])
                    mq_dup_df = mq_sort_df.drop_duplicates(subset=['sample_id'])

                    nc_expected = nc_dup_df['nc_cnt'].sum()
                    mq_expected = mq_dup_df['mq_cnt'].sum()

                    nc_sum_df = type_df.groupby(['cv_param1', 'cv_param2', 'cv_val1', 'cv_val2']
                                                )['nc_cnt'].sum().reset_index()
                    mq_sum_df = type_df.groupby(['cv_param1', 'cv_param2', 'cv_val1', 'cv_val2']
                                                )['mq_cnt'].sum().reset_index()
                    nc_sum_df['nc_expected'] = nc_expected
                    mq_sum_df['mq_expected'] = mq_expected
                    nc_sum_df['nc_exp_r'] = nc_sum_df['nc_cnt'] / nc_sum_df['nc_expected']
                    mq_sum_df['mq_exp_r'] = mq_sum_df['mq_cnt'] / mq_sum_df['mq_expected']
                    nc_val_max = nc_sum_df['nc_exp_r'].max()
                    mq_val_max = mq_sum_df['mq_exp_r'].max()
                    nc_params_df = nc_sum_df.query("nc_exp_r == @nc_val_max")
                    mq_params_df = mq_sum_df.query("mq_exp_r == @mq_val_max")
                    if nc_params_df.shape[0] > 1:
                        if cv_algo == 'hdbscan':
                            nc_sort_hdb_df = nc_params_df.sort_values(['cv_val1', 'cv_val2'],
                                                                      ascending=[True, True]
                                                                      )
                            nc_params_df = nc_sort_hdb_df.drop_duplicates(subset=['cv_val1', 'cv_val2'])
                        elif cv_algo == 'ocsvm':
                            nc_sort_oc_df = nc_params_df.sort_values(['cv_val1', 'cv_val2'],
                                                                     ascending=[False, True]
                                                                     )
                            nc_params_df = nc_sort_oc_df.drop_duplicates(subset=['cv_val1', 'cv_val2'])
                    if mq_params_df.shape[0] > 1:
                        if cv_algo == 'hdbscan':
                            mq_sort_hdb_df = mq_params_df.sort_values(['cv_val1', 'cv_val2'],
                                                                      ascending=[True, True]
                                                                      )
                            mq_params_df = mq_sort_hdb_df.drop_duplicates(subset=['cv_val1', 'cv_val2'])
                        elif cv_algo == 'ocsvm':
                            mq_sort_oc_df = mq_params_df.sort_values(['cv_val1', 'cv_val2'],
                                                                     ascending=[False, True]
                                                                     )
                            mq_params_df = mq_sort_oc_df.drop_duplicates(subset=['cv_val1', 'cv_val2'])
                    nc_cv_param1 = nc_params_df['cv_param1'].values[0]
                    nc_cv_param2 = nc_params_df['cv_param2'].values[0]
                    nc_cv_val1 = nc_params_df['cv_val1'].values[0]
                    nc_cv_val2 = nc_params_df['cv_val2'].values[0]
                    nc_best_df = type_df.query("cv_param1 == @nc_cv_param1 &"
                                               "cv_param2 == @nc_cv_param2 &"
                                               "cv_val1 == @nc_cv_val1 &"
                                               "cv_val2 == @nc_cv_val2"
                                               )
                    mq_cv_param1 = mq_params_df['cv_param1'].values[0]
                    mq_cv_param2 = mq_params_df['cv_param2'].values[0]
                    mq_cv_val1 = mq_params_df['cv_val1'].values[0]
                    mq_cv_val2 = mq_params_df['cv_val2'].values[0]
                    mq_best_df = type_df.query("cv_param1 == @mq_cv_param1 &"
                                               "cv_param2 == @mq_cv_param2 &"
                                               "cv_val1 == @mq_cv_val1 &"
                                               "cv_val2 == @mq_cv_val2"
                                               )
                    # SUM These df then add expected calc for grouping level
                    nc_group_df = nc_best_df.groupby(['cv_algo', 'algo', 'level', 'cv_param1',
                                                      'cv_param2', 'cv_val1', 'cv_val2']
                                                     )['nc_cnt'].sum().reset_index()
                    mq_group_df = mq_best_df.groupby(['cv_algo', 'algo', 'level', 'cv_param1',
                                                      'cv_param2', 'cv_val1', 'cv_val2']
                                                     )['mq_cnt'].sum().reset_index()
                    nc_group_df['group_val'] = val
                    mq_group_df['group_val'] = val
                    nc_group_df['grouping'] = group
                    mq_group_df['grouping'] = group
                    nc_group_df['nc_expected'] = nc_expected
                    mq_group_df['mq_expected'] = mq_expected
                    nc_group_df['nc_exp_r'] = nc_group_df['nc_cnt'] / nc_group_df['nc_expected']
                    mq_group_df['mq_exp_r'] = mq_group_df['mq_cnt'] / mq_group_df['mq_expected']
                    nc_agg_list.append(nc_group_df)
                    mq_agg_list.append(mq_group_df)
    nc_agg_df = pd.concat(nc_agg_list)
    mq_agg_df = pd.concat(mq_agg_list)
    nc_agg_df.to_csv(os.path.join(working_dir, "NC_agg_params.tsv"), sep='\t', index=False)
    mq_agg_df.to_csv(os.path.join(working_dir, "MQ_agg_params.tsv"), sep='\t', index=False)
    nc_group_df = nc_agg_df.groupby(['grouping', 'cv_algo', 'algo', 'level']
                                    )[['nc_cnt', 'nc_expected']].sum().reset_index()
    mq_group_df = mq_agg_df.groupby(['grouping', 'cv_algo', 'algo', 'level']
                                    )[['mq_cnt', 'mq_expected']].sum().reset_index()
    nc_group_df['nc_r'] = nc_group_df['nc_cnt'] / nc_group_df['nc_expected']
    mq_group_df['mq_r'] = mq_group_df['mq_cnt'] / mq_group_df['mq_expected']
    nc_group_df.to_csv(os.path.join(working_dir, "NC_agg_stats.tsv"), sep='\t', index=False)
    mq_group_df.to_csv(os.path.join(working_dir, "MQ_agg_stats.tsv"), sep='\t', index=False)

    return nc_agg_df, mq_agg_df


def best_match_params(real_dir, clust_all_df):
    real_df = pd.read_csv(os.path.join(real_dir, 'cluster_clean.tsv'), sep='\t', header=0)
    real_df = real_df[['sample_type', 'sample_id', 'best_match', 'best_cluster']]
    # best_match
    nc_max_df = clust_all_df.groupby(['sample_id', 'cv_algo', 'algo', 'level',
                                      'cv_param1', 'cv_param2', 'cv_val1', 'cv_val2']
                                     )[['nc_cnt']].max().reset_index()
    mq_max_df = clust_all_df.groupby(['sample_id', 'cv_algo', 'algo', 'level',
                                      'cv_param1', 'cv_param2', 'cv_val1', 'cv_val2']
                                     )[['mq_cnt']].max().reset_index()
    nc_hdb_df = nc_max_df.query("cv_algo == 'hdbscan' & (algo == 'hdbscan' | algo == 'denovo')")
    nc_ocs_df = nc_max_df.query("cv_algo == 'ocsvm' & algo == 'ocsvm'")
    mq_hdb_df = mq_max_df.query("cv_algo == 'hdbscan' & (algo == 'hdbscan' | algo == 'denovo')")
    mq_ocs_df = mq_max_df.query("cv_algo == 'ocsvm' & algo == 'ocsvm'")
    nc_hdb_sort_df = nc_hdb_df.sort_values(['nc_cnt', 'cv_val1', 'cv_val2'],
                                           ascending=[False, True, True]
                                           )
    nc_ocs_sort_df = nc_ocs_df.sort_values(['nc_cnt', 'cv_val1', 'cv_val2'],
                                           ascending=[False, False, True]
                                           )
    mq_hdb_sort_df = mq_hdb_df.sort_values(['mq_cnt', 'cv_val1', 'cv_val2'],
                                           ascending=[False, True, True]
                                           )
    mq_ocs_sort_df = mq_ocs_df.sort_values(['mq_cnt', 'cv_val1', 'cv_val2'],
                                           ascending=[False, False, True]
                                           )
    nc_hdb_dup_df = nc_hdb_sort_df.drop_duplicates(subset=['sample_id', 'cv_algo',
                                                           'algo', 'level']
                                                   )
    nc_ocs_dup_df = nc_ocs_sort_df.drop_duplicates(subset=['sample_id', 'cv_algo',
                                                           'algo', 'level']
                                                   )
    mq_hdb_dup_df = mq_hdb_sort_df.drop_duplicates(subset=['sample_id', 'cv_algo',
                                                           'algo', 'level']
                                                   )
    mq_ocs_dup_df = mq_ocs_sort_df.drop_duplicates(subset=['sample_id', 'cv_algo',
                                                           'algo', 'level']
                                                   )
    best_nc_hdb_df = real_df.merge(nc_hdb_dup_df, left_on='best_match',
                                   right_on='sample_id', how='left'
                                   )
    best_nc_ocs_df = real_df.merge(nc_ocs_dup_df, left_on='best_match',
                                   right_on='sample_id', how='left'
                                   )
    best_mq_hdb_df = real_df.merge(mq_hdb_dup_df, left_on='best_match',
                                   right_on='sample_id', how='left'
                                   )
    best_mq_ocs_df = real_df.merge(mq_ocs_dup_df, left_on='best_match',
                                   right_on='sample_id', how='left'
                                   )
    best_nc_hdb_df['mq_nc'] = 'nc'
    best_nc_ocs_df['mq_nc'] = 'nc'
    best_mq_hdb_df['mq_nc'] = 'mq'
    best_mq_ocs_df['mq_nc'] = 'mq'
    keep_cols = ['sample_type', 'sample_id_x', 'best_match', 'cv_algo', 'algo', 'level',
                 'cv_param1', 'cv_param2', 'cv_val1', 'cv_val2', 'mq_nc'
                 ]
    best_cat_df = pd.concat([best_nc_hdb_df[keep_cols], best_nc_ocs_df[keep_cols],
                             best_mq_hdb_df[keep_cols], best_mq_ocs_df[keep_cols]
                             ])
    best_dup_df = best_cat_df.drop_duplicates()
    best_dup_df.to_csv(os.path.join(real_dir, 'best_match_params.tsv'), sep='\t', index=False)

    return best_dup_df, real_df


def best_cluster_params(real_dir, nc_agg_df, mq_agg_df, real_df):
    nc_clust_df = nc_agg_df.query("grouping == 'best_cluster'")
    mq_clust_df = mq_agg_df.query("grouping == 'best_cluster'")
    nc_max_df = nc_clust_df.groupby(['group_val', 'cv_algo', 'algo', 'level',
                                     'cv_param1', 'cv_param2', 'cv_val1', 'cv_val2']
                                    )['nc_cnt'].max().reset_index()
    mq_max_df = mq_clust_df.groupby(['group_val', 'cv_algo', 'algo', 'level',
                                     'cv_param1', 'cv_param2', 'cv_val1', 'cv_val2']
                                    )['mq_cnt'].max().reset_index()
    nc_hdb_df = nc_max_df.query("cv_algo == 'hdbscan' & (algo == 'hdbscan' | algo == 'denovo')")
    nc_ocs_df = nc_max_df.query("cv_algo == 'ocsvm' & algo == 'ocsvm'")
    mq_hdb_df = mq_max_df.query("cv_algo == 'hdbscan' & (algo == 'hdbscan' | algo == 'denovo')")
    mq_ocs_df = mq_max_df.query("cv_algo == 'ocsvm' & algo == 'ocsvm'")
    nc_hdb_sort_df = nc_hdb_df.sort_values(['nc_cnt', 'cv_val1', 'cv_val2'],
                                           ascending=[False, True, True]
                                           )
    nc_ocs_sort_df = nc_ocs_df.sort_values(['nc_cnt', 'cv_val1', 'cv_val2'],
                                           ascending=[False, False, True]
                                           )
    mq_hdb_sort_df = mq_hdb_df.sort_values(['mq_cnt', 'cv_val1', 'cv_val2'],
                                           ascending=[False, True, True]
                                           )
    mq_ocs_sort_df = mq_ocs_df.sort_values(['mq_cnt', 'cv_val1', 'cv_val2'],
                                           ascending=[False, False, True]
                                           )
    nc_hdb_dup_df = nc_hdb_sort_df.drop_duplicates(subset=['group_val', 'cv_algo',
                                                           'algo', 'level']
                                                   )
    nc_ocs_dup_df = nc_ocs_sort_df.drop_duplicates(subset=['group_val', 'cv_algo',
                                                           'algo', 'level']
                                                   )
    mq_hdb_dup_df = mq_hdb_sort_df.drop_duplicates(subset=['group_val', 'cv_algo',
                                                           'algo', 'level']
                                                   )
    mq_ocs_dup_df = mq_ocs_sort_df.drop_duplicates(subset=['group_val', 'cv_algo',
                                                           'algo', 'level']
                                                   )
    best_nc_hdb_df = real_df.merge(nc_hdb_dup_df, left_on='best_cluster',
                                   right_on='group_val', how='left'
                                   )
    best_nc_ocs_df = real_df.merge(nc_ocs_dup_df, left_on='best_cluster',
                                   right_on='group_val', how='left'
                                   )
    best_mq_hdb_df = real_df.merge(mq_hdb_dup_df, left_on='best_cluster',
                                   right_on='group_val', how='left'
                                   )
    best_mq_ocs_df = real_df.merge(mq_ocs_dup_df, left_on='best_cluster',
                                   right_on='group_val', how='left'
                                   )
    best_nc_hdb_df['mq_nc'] = 'nc'
    best_nc_ocs_df['mq_nc'] = 'nc'
    best_mq_hdb_df['mq_nc'] = 'mq'
    best_mq_ocs_df['mq_nc'] = 'mq'
    keep_cols = ['sample_type', 'sample_id', 'best_cluster', 'cv_algo', 'algo', 'level',
                 'cv_param1', 'cv_param2', 'cv_val1', 'cv_val2', 'mq_nc'
                 ]
    best_cat_df = pd.concat([best_nc_hdb_df[keep_cols], best_nc_ocs_df[keep_cols],
                             best_mq_hdb_df[keep_cols], best_mq_ocs_df[keep_cols]
                             ])
    best_dup_df = best_cat_df.drop_duplicates()
    best_dup_df.to_csv(os.path.join(real_dir, 'cluster_params.tsv'), sep='\t', index=False)

    return best_dup_df


def majority_rule_params(real_dir, nc_agg_df, mq_agg_df, real_df):
    nc_major_df = nc_agg_df.query("grouping == 'majority_rule'")
    mq_major_df = mq_agg_df.query("grouping == 'majority_rule'")
    nc_max_df = nc_major_df.groupby(['group_val', 'cv_algo', 'algo', 'level',
                                     'cv_param1', 'cv_param2', 'cv_val1', 'cv_val2']
                                    )['nc_cnt'].max().reset_index()
    mq_max_df = mq_major_df.groupby(['group_val', 'cv_algo', 'algo', 'level',
                                     'cv_param1', 'cv_param2', 'cv_val1', 'cv_val2']
                                    )['mq_cnt'].max().reset_index()
    nc_hdb_df = nc_max_df.query("cv_algo == 'hdbscan' & (algo == 'hdbscan' | algo == 'denovo')")
    nc_ocs_df = nc_max_df.query("cv_algo == 'ocsvm' & algo == 'ocsvm'")
    mq_hdb_df = mq_max_df.query("cv_algo == 'hdbscan' & (algo == 'hdbscan' | algo == 'denovo')")
    mq_ocs_df = mq_max_df.query("cv_algo == 'ocsvm' & algo == 'ocsvm'")
    nc_hdb_sort_df = nc_hdb_df.sort_values(['nc_cnt', 'cv_val1', 'cv_val2'],
                                           ascending=[False, True, True]
                                           )
    nc_ocs_sort_df = nc_ocs_df.sort_values(['nc_cnt', 'cv_val1', 'cv_val2'],
                                           ascending=[False, False, True]
                                           )
    mq_hdb_sort_df = mq_hdb_df.sort_values(['mq_cnt', 'cv_val1', 'cv_val2'],
                                           ascending=[False, True, True]
                                           )
    mq_ocs_sort_df = mq_ocs_df.sort_values(['mq_cnt', 'cv_val1', 'cv_val2'],
                                           ascending=[False, False, True]
                                           )
    nc_hdb_dup_df = nc_hdb_sort_df.drop_duplicates(subset=['group_val', 'cv_algo',
                                                           'algo', 'level']
                                                   )
    nc_ocs_dup_df = nc_ocs_sort_df.drop_duplicates(subset=['group_val', 'cv_algo',
                                                           'algo', 'level']
                                                   )
    mq_hdb_dup_df = mq_hdb_sort_df.drop_duplicates(subset=['group_val', 'cv_algo',
                                                           'algo', 'level']
                                                   )
    mq_ocs_dup_df = mq_ocs_sort_df.drop_duplicates(subset=['group_val', 'cv_algo',
                                                           'algo', 'level']
                                                   )
    nc_hdb_dup_df['tmp_key'] = 0
    nc_ocs_dup_df['tmp_key'] = 0
    mq_hdb_dup_df['tmp_key'] = 0
    mq_ocs_dup_df['tmp_key'] = 0
    real_df['tmp_key'] = 0
    best_nc_hdb_df = real_df.merge(nc_hdb_dup_df, on='tmp_key', how='left')
    best_nc_ocs_df = real_df.merge(nc_ocs_dup_df, on='tmp_key', how='left')
    best_mq_hdb_df = real_df.merge(mq_hdb_dup_df, on='tmp_key', how='left')
    best_mq_ocs_df = real_df.merge(mq_ocs_dup_df, on='tmp_key', how='left')
    best_nc_hdb_df['mq_nc'] = 'nc'
    best_nc_ocs_df['mq_nc'] = 'nc'
    best_mq_hdb_df['mq_nc'] = 'mq'
    best_mq_ocs_df['mq_nc'] = 'mq'
    keep_cols = ['sample_type', 'sample_id', 'cv_algo', 'algo', 'level',
                 'cv_param1', 'cv_param2', 'cv_val1', 'cv_val2', 'mq_nc'
                 ]
    best_cat_df = pd.concat([best_nc_hdb_df[keep_cols], best_nc_ocs_df[keep_cols],
                             best_mq_hdb_df[keep_cols], best_mq_ocs_df[keep_cols]
                             ])
    best_dup_df = best_cat_df.drop_duplicates()
    best_dup_df.to_csv(os.path.join(real_dir, 'majority_rule_params.tsv'), sep='\t', index=False)

    return best_dup_df


def run_param_match(working_dir, real_dir):
    sample_list = glob.glob(os.path.join(working_dir, "*/CV_*.tsv"))
    config_list = [['hdbscan', 'denovo'], ['hdbscan', 'hdbscan'], ['ocsvm', 'ocsvm']]
    cv_concat_df, cv_all_df = build_cv_tables(working_dir, sample_list, config_list)
    clust_cv_df, clust_all_df = merge_cluster_table(working_dir, cv_concat_df, cv_all_df)
    nc_agg_df, mq_agg_df = aggregate_best_params(working_dir, config_list, clust_all_df)
    # Need to get best params for best_match, best_cluster, and majority_rule
    # Assign the best match params to the SI data
    best_match_df, real_df = best_match_params(real_dir, clust_all_df)
    # best_cluster
    best_cluster_df = best_cluster_params(real_dir, nc_agg_df, mq_agg_df, real_df)
    # majority_rule
    majority_rule_df = majority_rule_params(real_dir, nc_agg_df, mq_agg_df, real_df)

    return best_match_df, best_cluster_df, majority_rule_df


############################################################################################
ref_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                       'renyi_entropy/SCARAB_CV/'
                       )
real_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                        'renyi_entropy/SI/'
                        )
best_match_df, best_cluster_df, majority_rule_df = run_param_match(ref_dir, real_dir)

'''
cv_file = sys.argv[1]
cv_df = pd.read_csv(cv_file, sep='\t', header=0)
sv_dir = os.path.split(cv_file)[0] + '/'
for cv_algo in cv_df['cv_algo'].unique():
    cv_algo_df = cv_df.query('cv_algo == @cv_algo')
    for algo in cv_algo_df['algo'].unique():
        algo_df = cv_algo_df.query('algo == @algo')
        if ((cv_algo == 'hdbscan' and (algo == 'hdbscan' or algo == 'denovo')) or
                (cv_algo == 'ocsvm' and algo == 'ocsvm')
        ):
            for level in algo_df['level'].unique():
                level_df = algo_df.query('level == @level')
                fig = go.Figure(data=go.Heatmap(
                    z=level_df['mq_cnt'],
                    x=[str(x) for x in level_df['cv_val1']],
                    y=[str(x) for x in level_df['cv_val2']],
                    hoverongaps=False, colorscale='Viridis'))
                sv_file = sv_dir + '_'.join([cv_algo, algo, level]) + '.MQ_BINS.png'
                fig.write_image(sv_file)
                fig2 = go.Figure(data=go.Heatmap(
                    z=level_df['nc_cnt'],
                    x=[str(x) for x in level_df['cv_val1']],
                    y=[str(x) for x in level_df['cv_val2']],
                    hoverongaps=False, colorscale='Viridis'))
                sv2_file = sv_dir + '_'.join([cv_algo, algo, level]) + '.NC_BINS.png'
                fig2.write_image(sv2_file)
'''
