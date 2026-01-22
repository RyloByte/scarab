import glob
import os

import hdbscan
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import umap
from scarab.entropy import renyi_entropy
from scipy import sparse
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.preprocessing import StandardScaler


def calc_ref_entropy(sample_list, ref_dir, rerun_ref):
    entropy_file = os.path.join(ref_dir, 'entropy_table.tsv')
    if rerun_ref:
        # Calculate entropy for all references
        entropy_list = []
        for samp_file in sample_list:
            samp_id = samp_file.split('/')[-1].rsplit('.', 1)[0]
            print(samp_id)
            if (('entropy' not in samp_id) and ('table' not in samp_id)
                    and ('cluster' not in samp_id)
            ):
                if samp_id.rsplit('_', 1)[1][0].isdigit():
                    samp_label = samp_id.rsplit('_', 1)[0]
                    samp_rep = samp_id.rsplit('_', 1)[1]
                else:
                    samp_label = samp_id.rsplit('_', 1)[0]
                    samp_rep = samp_id.rsplit('_', 1)[1][-1]
                cov_df = pd.read_csv(samp_file, sep='\t', header=0)
                depth_sum = cov_df['totalAvgDepth'].sum()
                relative_depth = [x / depth_sum for x in cov_df['totalAvgDepth']]
                probs = relative_depth
                q_list = [0, 1, 2, 4, 8, 16, 32, np.inf]
                for q in q_list:
                    r_ent = renyi_entropy(probs, q)
                    print(samp_id, q, r_ent)
                    entropy_list.append([samp_id, samp_label, samp_rep, q, r_ent])
                '''
                c = 0
                for col in cov_df.columns:
                    if col.split('.')[-1] == 'bam':
                        samp_id_c = samp_id # + '_' + str(c)
                        depth_sum = cov_df[col].sum()
                        relative_depth = [x/depth_sum for x in cov_df[col]]
                        probs = relative_depth
                        q_list = [0, 1, 2, 4, 8, 16, 32, np.inf]
                        for q in q_list:
                            r_ent = renyi_entropy(probs, q)
                            print(samp_id_c, q, r_ent)
                            entropy_list.append([samp_id_c, samp_label, samp_rep, q, r_ent])
                        c += 1
                '''
        ent_df = pd.DataFrame(entropy_list, columns=['sample_id', 'sample_type',
                                                     'sample_rep', 'alpha',
                                                     'Renyi_Entropy'
                                                     ])
        # Have to replace the np.inf with a real value for plotting
        ent_df['alpha_int'] = ent_df['alpha'].replace({4: 3, 8: 4, 16: 5, 32: 6, np.inf: 7}).astype(int)
        x_labels = {0: 'Richness (a=0)', 1: 'Shannon (a=1)',
                    2: 'Simpson (a=2)', 3: '4', 4: '8', 5: '16',
                    6: '32', 7: 'Berger–Parker (a=inf)'
                    }
        ent_df['x_labels'] = ent_df['alpha_int'].map(x_labels)
        ent_df.to_csv(entropy_file, sep='\t', index=False)
    else:
        ent_df = pd.read_csv(entropy_file, sep='\t', header=0)

    return ent_df


def entropy_cluster(working_dir, ent_df):
    cluster_table = os.path.join(working_dir, 'cluster_table.tsv')
    samp2type = {x: y for x, y in zip(ent_df['sample_id'], ent_df['sample_type'])}
    piv_df = ent_df.pivot(index='sample_id', columns='alpha', values='Renyi_Entropy')
    scaler = StandardScaler()
    scale_fit = scaler.fit(piv_df)
    scale_emb = scale_fit.transform(piv_df)
    scale_df = pd.DataFrame(scale_emb, index=piv_df.index.values,
                            columns=piv_df.columns
                            )
    umap_fit = umap.UMAP(n_neighbors=2, min_dist=0.0, n_components=2,
                         random_state=42
                         ).fit(scale_df)
    umap_emb = umap_fit.transform(scale_df)
    umap_df = pd.DataFrame(umap_emb, index=scale_df.index.values, columns=['u0', 'u1'])
    clusterer = hdbscan.HDBSCAN(min_cluster_size=2, allow_single_cluster=False,
                                prediction_data=True
                                ).fit(umap_df)

    cluster_labels = clusterer.labels_
    cluster_probs = clusterer.probabilities_
    cluster_outlier = clusterer.outlier_scores_

    umap_df['sample_type'] = [samp2type[x] for x in umap_df.index.values]
    umap_df['sample_id'] = umap_df.index.values
    umap_df['cluster'] = cluster_labels
    umap_df['probabilities'] = cluster_probs
    umap_df['outlier_scores'] = cluster_outlier

    ent_umap_df = ent_df.merge(umap_df, on=['sample_id', 'sample_type'], how='left')
    ent_best_df = find_best_match(piv_df, ent_umap_df)
    ent_best_df.to_csv(cluster_table, sep='\t', index=False)

    return umap_df, ent_best_df, piv_df, umap_fit, clusterer, scale_fit


def find_best_match(piv_df, ent_umap_df):
    # Closest ref sample methods
    outlier_list = list(ent_umap_df.query("cluster == -1")['sample_id'])
    cmpr_list = []
    for r1, row1 in piv_df.iterrows():
        keep_diff = [r1, '', np.inf]
        for r2, row2 in piv_df.iterrows():
            euc_d = np.linalg.norm(row1 - row2)
            if euc_d < keep_diff[2] and r1 != r2 and r2 not in outlier_list:
                keep_diff = [r1, r2, euc_d]
        cmpr_list.append(keep_diff)
    cmpr_df = pd.DataFrame(cmpr_list, columns=['sample_id', 'best_match', 'euc_d'])
    ent_best_df = ent_umap_df.merge(cmpr_df, on='sample_id', how='left')

    return ent_best_df


def real_best_match(piv_df, real_piv_df, real_umap_df, working_dir):
    # Closest ref sample methods
    r_cmpr_list = []
    for r1, row1 in real_piv_df.iterrows():
        keep_diff = [r1, '', np.inf]
        for r2, row2 in piv_df.iterrows():
            euc_d = np.linalg.norm(row1 - row2)
            if euc_d < keep_diff[2] and r1 != r2:
                keep_diff = [r1, r2, euc_d]
        r_cmpr_list.append(keep_diff)
    r_cmpr_df = pd.DataFrame(r_cmpr_list, columns=['sample_id', 'best_match', 'euc_d'])
    best_df = real_umap_df.merge(r_cmpr_df, on='sample_id', how='left')
    best_df.to_csv(os.path.join(working_dir, 'cluster_table.tsv'), sep='\t', index=False)
    return best_df


def real_cluster(clusterer, real_df, umap_fit, scale_fit):
    # Assign real data to clusters
    real_piv_df = real_df.pivot(index='sample_id', columns='alpha', values='Renyi_Entropy')
    scale_emb = scale_fit.transform(real_piv_df)
    scale_df = pd.DataFrame(scale_emb, index=real_piv_df.index.values, columns=real_piv_df.columns)
    umap_emb = umap_fit.transform(scale_df)
    umap_df = pd.DataFrame(umap_emb, index=scale_df.index.values, columns=['u0', 'u1'])
    test_labels, strengths = hdbscan.approximate_predict(clusterer, umap_df)
    cluster_labels = test_labels
    cluster_probs = strengths
    samp2type = {x: y for x, y in zip(real_df['sample_id'], real_df['sample_type'])}
    umap_df['sample_type'] = [samp2type[x] for x in umap_df.index.values]
    umap_df['sample_id'] = umap_df.index.values
    umap_df['cluster'] = cluster_labels
    umap_df['probabilities'] = cluster_probs
    return real_piv_df, umap_df


def calc_real_entrophy(mba_cov_list, working_dir):
    entropy_list = []
    for samp_file in mba_cov_list:
        samp_id = samp_file.split('/')[-1].rsplit('.', 1)[0]
        print(samp_id)
        if samp_id.rsplit('_', 1)[1][0].isdigit():
            samp_label = samp_id.rsplit('_', 1)[0]
            samp_rep = samp_id.rsplit('_', 1)[1]
        else:
            samp_label = samp_id
            samp_rep = 0
        cov_df = pd.read_csv(samp_file, sep='\t', header=0)
        depth_sum = cov_df['totalAvgDepth'].sum()
        cov_df['relative_depth'] = [x / depth_sum for x in cov_df['totalAvgDepth']]
        probs = cov_df['relative_depth'].tolist()
        q_list = [0, 1, 2, 4, 8, 16, 32, np.inf]
        for q in q_list:
            r_ent = renyi_entropy(probs, q)
            print(q, r_ent)
            entropy_list.append([samp_id, samp_label, samp_rep, q, r_ent])
    real_df = pd.DataFrame(entropy_list, columns=['sample_id', 'sample_type',
                                                  'sample_rep', 'alpha',
                                                  'Renyi_Entropy'
                                                  ])
    # Have to replace the np.inf with a real value for plotting
    real_df['alpha_int'] = real_df['alpha'].replace({4: 3, 8: 4, 16: 5, 32: 6, np.inf: 7}).astype(int)
    x_labels = {0: 'Richness (a=0)', 1: 'Shannon (a=1)',
                2: 'Simpson (a=2)', 3: '4', 4: '8', 5: '16',
                6: '32', 7: 'Berger–Parker (a=inf)'
                }
    real_df['x_labels'] = real_df['alpha_int'].map(x_labels)
    real_df.to_csv(os.path.join(working_dir, 'entropy_table.tsv'), sep='\t', index=False)

    return real_df


def remove_outliers(ent_best_df, real_merge_df, umap_fit, scale_fit):
    # If labeled as an outlier, take the closest match
    keep_cols = ['sample_id', 'sample_type', 'alpha', 'Renyi_Entropy', 'alpha_int',
                 'x_labels', 'u0', 'u1', 'cluster', 'probabilities', 'best_match', 'euc_d'
    best_merge_df = pd.concat([ent_best_df[keep_cols],
                               real_merge_df[keep_cols]]
                              )
    # Calculate centroids for cluster
    real_id_list = real_merge_df['sample_id']
    cent_merge_df = calc_centroid(best_merge_df,
                                  best_merge_df[['sample_id', 'alpha',
                                                 'Renyi_Entropy',
                                                 'cluster']],
                                  real_id_list, umap_fit, scale_fit
                                  )

    return cent_merge_df


def calc_centroid(best_df, ren_df, real_ids, umap_fit, scale_fit):
    ren_piv_df = ren_df.pivot(index=['sample_id', 'cluster'], columns='alpha', values='Renyi_Entropy').reset_index()
    sample_ren_df = ren_piv_df.query('cluster == -1').drop(['cluster'], axis=1)
    sample_ren_df.set_index('sample_id', inplace=True)
    samp_scale_emb = scale_fit.transform(sample_ren_df)
    samp_scale_df = pd.DataFrame(samp_scale_emb, index=sample_ren_df.index.values, columns=sample_ren_df.columns)
    sample_ren_emb = umap_fit.transform(samp_scale_df)
    sample_umap_df = pd.DataFrame(sample_ren_emb, index=samp_scale_df.index.values, columns=['u0', 'u1'])
    cluster_ren_df = ren_piv_df.query('cluster != -1 and sample_id not in @real_ids'
                                      ).drop(['sample_id'], axis=1)
    cluster_ren_df.set_index('cluster', inplace=True)
    clust_scale_emb = scale_fit.transform(cluster_ren_df)
    clust_scale_df = pd.DataFrame(clust_scale_emb, index=cluster_ren_df.index.values, columns=cluster_ren_df.columns)
    cluster_ren_emb = umap_fit.transform(clust_scale_df)
    cluster_umap_df = pd.DataFrame(cluster_ren_emb, index=clust_scale_df.index.values, columns=['u0', 'u1'])
    median_df = cluster_umap_df.groupby(by=cluster_umap_df.index).median()
    out_dict = {}
    for r1, row1 in sample_umap_df.iterrows():
        keep_diff = [r1, '', np.inf]
        for r2, row2 in median_df.iterrows():
            euc_d = np.linalg.norm(row1 - row2)
            if euc_d < keep_diff[2] and r1 != r2:
                keep_diff = [r1, r2, euc_d]
        out_dict[r1] = keep_diff
    best_df['best_cluster'] = [out_dict[x][1] if x in
                                                 out_dict.keys() else y for x, y
                               in zip(best_df['sample_id'],
                                      best_df['cluster']
                                      )]

    return best_df


def plot_renyi_entropy(working_dir, ent_df, real_type):
    # Build all the reference plots and run clustering
    type_list = ent_df['sample_type'].unique()
    cpal = {x: y for x, y in zip(type_list, sns.color_palette(n_colors=len(type_list)))}
    cpal[real_type] = 'black'
    sns.set(rc={'figure.figsize': (12, 8)})
    sns.set_style("white")
    p = sns.catplot(x="x_labels", y="Renyi_Entropy", hue="sample_type", col="sample_type",
                    data=ent_df, palette=cpal, col_wrap=3, legend_out=True
                    )
    p.map_dataframe(sns.boxplot, x="x_labels", y="Renyi_Entropy",
                    data=ent_df, boxprops={'facecolor': 'None'},
                    whiskerprops={'linewidth': 1},
                    showfliers=False
                    )
    p.set_xticklabels(rotation=45)
    lgd_dat = p._legend_data
    p.add_legend(legend_data={x: lgd_dat[x] for x in lgd_dat.keys() if x in cpal.keys()})
    p.savefig(os.path.join(working_dir, 'entropy_plot.png'), bbox_inches='tight')
    plt.clf()
    plt.close()

    return cpal


def plot_euc_dist(working_dir, ent_df, cpal):
    # Calculate euclidean distances and plot
    samp2type = {x: y for x, y in zip(ent_df['sample_id'], ent_df['sample_type'])}
    piv_df = ent_df.pivot(index='sample_id', columns='alpha', values='Renyi_Entropy')
    A_sparse = sparse.csr_matrix(piv_df)
    similarities = euclidean_distances(A_sparse)
    euc_df = pd.DataFrame(similarities, index=piv_df.index.values, columns=piv_df.index.values)
    euc_unstack_df = euc_df.unstack().reset_index()
    euc_unstack_df.columns = ['sample1', 'sample2', 'euclidean_distances']
    euc_unstack_df['samp1_label'] = [samp2type[x] for x in euc_unstack_df['sample1']]
    euc_unstack_df['samp2_label'] = [samp2type[x] for x in euc_unstack_df['sample2']]
    sns.set(rc={'figure.figsize': (12, 8)})
    sns.set_style("white")
    b = sns.catplot(x="samp1_label", y="euclidean_distances", hue="samp2_label", col="samp2_label",
                    data=euc_unstack_df, palette=cpal, col_wrap=3, legend_out=True
                    )
    b.map_dataframe(sns.boxplot, x="samp1_label", y="euclidean_distances",
                    data=euc_unstack_df, boxprops={'facecolor': 'None'},
                    whiskerprops={'linewidth': 1},
                    showfliers=False
                    )
    b.set_xticklabels(rotation=45)
    lgd_dat = b._legend_data
    b.add_legend(legend_data={x: lgd_dat[x] for x in lgd_dat.keys() if x in cpal.keys()})
    b.savefig(os.path.join(working_dir, 'similarity_plot.png'), bbox_inches='tight')
    plt.clf()
    plt.close()


def plot_ent_clust(working_dir, ent_umap_df, cpal):
    sns.set(rc={'figure.figsize': (12, 8)})
    sns.set_style("white")
    mark_list = ['.', 'v', '^', '<', '>', '8', 's', 'p', 'P', '*', 'H', 'X', 'D', 'd', 'o']
    ent_umap_df['0_jit'] = jitter(ent_umap_df['u0'], 0)
    ent_umap_df['1_jit'] = jitter(ent_umap_df['u1'], 0)
    ent_dup_df = ent_umap_df.drop_duplicates(subset=['sample_id'])
    b = sns.scatterplot(x='0_jit', y='1_jit', hue="sample_type", style='best_cluster',
                        data=ent_dup_df, palette=cpal,
                        markers=mark_list[:len(ent_dup_df['best_cluster'].unique())], s=200
                        )
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    fig = b.get_figure()
    fig.savefig(os.path.join(working_dir, 'entropy_clusters.png'), bbox_inches='tight')
    plt.clf()
    plt.close()

    sns.set(rc={'figure.figsize': (12, 8)})
    sns.set_style("white")
    p = sns.catplot(x="x_labels", y="Renyi_Entropy", hue="sample_type",
                    col="best_cluster", data=ent_umap_df, col_wrap=4, palette=cpal,
                    legend=False
                    )
    p.map_dataframe(sns.boxplot, x="x_labels", y="Renyi_Entropy",
                    data=ent_umap_df, boxprops={'facecolor': 'None'},
                    whiskerprops={'linewidth': 1},
                    showfliers=False
                    )
    p.set_xticklabels(rotation=45)
    lgd_dat = p._legend_data
    p.add_legend(legend_data={x: lgd_dat[x] for x in lgd_dat.keys() if x in cpal.keys()})
    p.savefig(os.path.join(working_dir, 'clustent_plot.png'), bbox_inches='tight')
    plt.clf()
    plt.close()


def jitter(values, j):
    return values + np.random.normal(j, 0.5, values.shape)


def calc_entropy(working_dir, mba_cov_list):
    rerun_ref = True  # to re-calc reference profiles, set to True
    make_plots = True  # if you want plots, set to True

    # Calculate entropy for all references
    ref_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           'renyi_entropy/references/'
                           )
    sample_list = glob.glob(os.path.join(ref_dir, "*.tsv"))
    ent_df = calc_ref_entropy(sample_list, ref_dir, rerun_ref)
    ent_results = entropy_cluster(ref_dir, ent_df)
    umap_df = ent_results[0]
    ent_best_df = ent_results[1]
    piv_df = ent_results[2]
    umap_fit = ent_results[3]
    clusterer = ent_results[4]
    scale_fit = ent_results[5]

    # Cluster real samples
    real_df = calc_real_entrophy(mba_cov_list, working_dir)
    real_piv_df, real_umap_df = real_cluster(clusterer, real_df, umap_fit, scale_fit)
    real_best_df = real_best_match(piv_df, real_piv_df, real_umap_df, working_dir)
    real_merge_df = real_df.merge(real_best_df, on=['sample_id', 'sample_type'], how='left')

    # Replace outliers with best match
    best_merge_df = remove_outliers(ent_best_df, real_merge_df, umap_fit, scale_fit)

    # Build plots
    real_type = real_df['sample_type'].values[0]
    if make_plots == True:
        cpal = plot_renyi_entropy(ref_dir, best_merge_df, real_type)
        plot_euc_dist(ref_dir, best_merge_df, cpal)
        plot_ent_clust(ref_dir, best_merge_df, cpal)

    # Save final files
    ref_only_df = best_merge_df.query('sample_type != @real_type')
    ref_clean = os.path.join(ref_dir, 'cluster_clean.tsv')
    ref_only_df.to_csv(ref_clean, sep='\t', index=False)
    real_only_df = best_merge_df.query('sample_type == @real_type')
    real_clean = os.path.join(working_dir, 'cluster_clean.tsv')
    real_only_df.to_csv(real_clean, sep='\t', index=False)

    return real_only_df


###############################################################################################
###############################################################################################
working_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                           'renyi_entropy/SI/'
                           )
mba_cov_list = glob.glob(os.path.join(working_dir, "SI*.tsv"))
best_params_df = calc_entropy(working_dir, mba_cov_list)
