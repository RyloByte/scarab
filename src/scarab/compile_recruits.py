import logging
from os.path import join as o_join
import subprocess
import pandas as pd
import scarab.utilities as s_utils

logger = logging.getLogger(__name__)


def run_combine_recruits(save_dirs_dict, mg_file, clusters,
                         trusted_list, mode, threads
                         ):
    denovo_clust_df = clusters[0]
    denovo_sv_path = save_dirs_dict['denovo']
    mode_path = save_dirs_dict[mode]

    logger.info('Combining recruits.')
    mg_contigs_dict = s_utils.get_seqs(mg_file)
    mg_contigs = tuple([(r[0], r[1]) for r in mg_contigs_dict])
    mg_contigs_df = pd.DataFrame(mg_contigs, columns=['contig_id', 'seq'])
    trust_dict = {t[0]: t[1] for t in trusted_list}
    
    # De Novo Bins
    denovo_set = list(set(denovo_clust_df['best_label']))
    for best_label in denovo_set:
        sub_merge_df = denovo_clust_df[['best_label', 'contig_id']
        ].query('best_label == @best_label').drop_duplicates()
        logger.info(
            'Recruited %s contigs from de novo analysis for %s.',
            sub_merge_df.shape[0],
            best_label,
        )
        final_rec = o_join(denovo_sv_path, str(best_label) + '.denovo.fasta')
        with open(final_rec, 'w') as final_out:
            contig_list = list(set(sub_merge_df['contig_id']))
            mg_sub_filter_df = mg_contigs_df.query('contig_id in @contig_list')
            final_mgsubs_list = ['\n'.join(['>' + x[0], x[1]]) for x in
                                 zip(mg_sub_filter_df['contig_id'],
                                     mg_sub_filter_df['seq']
                                     )
                                 ]
            final_out.write('\n'.join(final_mgsubs_list))

    ######################################
    ########## ANCHORED BINNING ##########
    ######################################

    trusted_clust_df = clusters[1]
    ocsvm_clust_df = clusters[2]
    inter_clust_df = clusters[3]

    hdbscan_sv_path = save_dirs_dict['hdbscan']
    ocsvm_sv_path = save_dirs_dict['ocsvm']
    inter_sv_path = save_dirs_dict['intersect']
    xpg_sv_path = save_dirs_dict['xpgs']

    if isinstance(inter_clust_df, pd.DataFrame):
        # HDBSCAN Bins
        hdbscan_set = list(set(trusted_clust_df['best_label']))
        logger.info('Found %s HDBSCAN bins.', len(hdbscan_set))
        for best_label in hdbscan_set:
            sub_merge_df = trusted_clust_df[['best_label', 'contig_id']
            ].query('best_label == @best_label').drop_duplicates()
            logger.info(
                'Recruited %s contigs from HDBSCAN anchored analysis for %s.',
                sub_merge_df.shape[0],
                best_label,
            )
            final_rec = o_join(hdbscan_sv_path, str(best_label) + '.hdbscan.fasta')
            with open(final_rec, 'w') as final_out:
                contig_list = list(set(sub_merge_df['contig_id']))
                mg_sub_filter_df = mg_contigs_df.query('contig_id in @contig_list')
                final_mgsubs_list = ['\n'.join(['>' + x[0], x[1]]) for x in
                                     zip(mg_sub_filter_df['contig_id'],
                                         mg_sub_filter_df['seq']
                                         )
                                     ]
                final_out.write('\n'.join(final_mgsubs_list))
        # Combine final recruits and reference trusted contigs
        logger.info('Running BBtools dedup on HDBSCAN bins.')
        for t_id in s_utils.progress(hdbscan_set):
            t_file = trust_dict[t_id]
            concat_file = o_join(xpg_sv_path, t_id + '.hdbscan.concat.fasta')
            with open(concat_file, 'w') as cat_out:
                data = []
                with open(t_file, 'r') as t_in:
                    data.extend(t_in.readlines())
                hdbscan_bin = o_join(hdbscan_sv_path, t_id + '.hdbscan.fasta')
                with open(hdbscan_bin, 'r') as r_file:
                    data.extend(r_file.readlines())
                join_data = '\n'.join(data).replace('\n\n', '\n')
                cat_out.write(join_data)

            # Use BBTools dedupe.sh to deduplicate the extend SAG file
            dedupe_fa = o_join(xpg_sv_path, t_id + '.hdbscan.xPG.fasta')
            dedupe_cmd = ['dedupe.sh', 'in=' + concat_file, 'out=' + dedupe_fa,
                          'threads=' + str(threads), 'minidentity=97', 'overwrite=true']
            subprocess.run(dedupe_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        # OC-SVM Bins
        ocsvm_set = list(set(ocsvm_clust_df['best_label']))
        for best_label in ocsvm_set:
            sub_merge_df = ocsvm_clust_df[['best_label', 'contig_id']
            ].query('best_label == @best_label').drop_duplicates()
            logger.info(
                'Recruited %s contigs from OC-SVM anchored analysis for %s.',
                sub_merge_df.shape[0],
                best_label,
            )
            final_rec = o_join(ocsvm_sv_path, str(best_label) + '.ocsvm.fasta')
            with open(final_rec, 'w') as final_out:
                contig_list = list(set(sub_merge_df['contig_id']))
                mg_sub_filter_df = mg_contigs_df.query('contig_id in @contig_list')
                final_mgsubs_list = ['\n'.join(['>' + x[0], x[1]]) for x in
                                     zip(mg_sub_filter_df['contig_id'],
                                         mg_sub_filter_df['seq']
                                         )
                                     ]
                final_out.write('\n'.join(final_mgsubs_list))
        # Combine final recruits and reference trusted contigs
        logger.info('Running BBtools dedup on OC-SVM bins.')
        for t_id in s_utils.progress(ocsvm_set):
            t_file = trust_dict[t_id]
            concat_file = o_join(xpg_sv_path, t_id + '.ocsvm.concat.fasta')
            with open(concat_file, 'w') as cat_out:
                data = []
                with open(t_file, 'r') as t_in:
                    data.extend(t_in.readlines())
                ocsvm_bin = o_join(ocsvm_sv_path, t_id + '.ocsvm.fasta')
                with open(ocsvm_bin, 'r') as r_file:
                    data.extend(r_file.readlines())
                join_data = '\n'.join(data).replace('\n\n', '\n')
                cat_out.write(join_data)

            # Use BBTools dedupe.sh to deduplicate the extend SAG file
            dedupe_fa = o_join(xpg_sv_path, t_id + '.ocsvm.xPG.fasta')
            dedupe_cmd = ['dedupe.sh', 'in=' + concat_file, 'out=' + dedupe_fa,
                          'threads=' + str(threads), 'minidentity=97', 'overwrite=true']
            subprocess.run(dedupe_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)    
        # Combined Bins
        inter_set = list(set(inter_clust_df['best_label']))
        if isinstance(inter_clust_df, pd.DataFrame):
            for best_label in inter_set:
                sub_merge_df = inter_clust_df.query('best_label == @best_label')
                logger.info(
                    'Recruited %s contigs from intersection analysis for %s.',
                    sub_merge_df.shape[0],
                    best_label,
                )
                final_rec = o_join(inter_sv_path, str(best_label) + '.intersect.fasta')
                with open(final_rec, 'w') as final_out:
                    contig_list = list(set(sub_merge_df['contig_id']))
                    mg_sub_filter_df = mg_contigs_df.query('contig_id in @contig_list')
                    final_mgsubs_list = ['\n'.join(['>' + x[0], x[1]]) for x in
                                         zip(mg_sub_filter_df['contig_id'],
                                             mg_sub_filter_df['seq']
                                             )
                                         ]
                    final_out.write('\n'.join(final_mgsubs_list))
            # Combine final recruits and reference trusted contigs
            logger.info('Running BBtools dedup on intersection bins.')
            for t_id in s_utils.progress(inter_set):
                t_file = trust_dict[t_id]
                concat_file = o_join(xpg_sv_path, t_id + '.intersect.concat.fasta')
                with open(concat_file, 'w') as cat_out:
                    data = []
                    with open(t_file, 'r') as t_in:
                        data.extend(t_in.readlines())
                    recruit_bin = o_join(inter_sv_path, t_id + '.intersect.fasta')
                    with open(recruit_bin, 'r') as r_file:
                        data.extend(r_file.readlines())
                    join_data = '\n'.join(data).replace('\n\n', '\n')
                    cat_out.write(join_data)

                # Use BBTools dedupe.sh to deduplicate the extend SAG file
                dedupe_fa = o_join(xpg_sv_path, t_id + '.intersect.xPG.fasta')
                dedupe_cmd = ['dedupe.sh', 'in=' + concat_file, 'out=' + dedupe_fa,
                              'threads=' + str(threads), 'minidentity=97', 'overwrite=true']
                subprocess.run(dedupe_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # Clean up the directory
    logger.info('Cleaning up intermediate files.')
    #s_utils.runCleaner(mode_path, "denovo")
    #s_utils.runCleaner(mode_path, "hdbscan")
    #s_utils.runCleaner(mode_path, "ocsvm")
    #s_utils.runCleaner(mode_path, "intersect")
    s_utils.runCleaner(xpg_sv_path, "*.concat.fasta")
