"""Command handlers for the SCARAB CLI."""

__author__ = 'Ryan J McLaughlin'

import logging
import os

logger = logging.getLogger(__name__)


def _python_dependency_versions():
    """Return installed Python package versions as a name -> version mapping."""
    from pip._internal.operations import freeze

    return {x.split('==')[0]: x.split('==')[1] for x in freeze.freeze()}


def _parse_recruit_args(sys_args):
    """Parse CLI arguments for the recruit subcommand."""
    import scarab.s_args as s_args

    parser = s_args.ScarabArgumentParser(description='Recruit environmental reads to reference contigs.')
    parser.add_recruit_args()
    return parser.parse_args(sys_args)


def _configure_recruit_base(args):
    """Create and populate a SCARAB base object with recruit settings."""
    import scarab.classy as s_class

    recruit_s = s_class.ScarabBase('recruit')
    recruit_s.trust_path = args.trust_path
    recruit_s.mg_file = args.mg_file
    recruit_s.mg_raw_file_list = args.mg_raw_file_list
    recruit_s.pacbio = args.pacbio
    recruit_s.save_path = args.save_path
    recruit_s.max_contig_len = int(args.max_contig_len)
    recruit_s.overlap_len = int(args.overlap_len)
    recruit_s.min_len = int(args.min_len)
    recruit_s.kmer_size = int(args.kmer_size)
    recruit_s.jaccard = float(args.jaccard)
    recruit_s.nthreads = int(args.nthreads)
    recruit_s.force = args.force
    recruit_s.denovo_min_clust, recruit_s.denovo_min_samp = args.denovo_min_clust, args.denovo_min_samp
    recruit_s.anchor_min_clust, recruit_s.anchor_min_samp = args.anchor_min_clust, args.anchor_min_samp
    recruit_s.nu, recruit_s.gamma = args.nu, args.gamma
    recruit_s.vr, recruit_s.r = args.vr_params, args.r_params
    recruit_s.s, recruit_s.vs = args.s_params, args.vs_params
    recruit_s.a = args.auto_params
    return recruit_s


def _select_mode(recruit_s):
    """Select a clustering mode from available flags."""
    mode_list = [recruit_s.vr, recruit_s.r, recruit_s.s, recruit_s.vs]
    mode = 'algo_defaults'
    for mode_candidate in mode_list:
        if mode_candidate:
            mode = mode_candidate
    return mode


def _build_metagenome_subcontigs(recruit_s):
    """Build subcontigs for the metagenome input."""
    import scarab.utilities as s_utils

    mg_id = os.path.splitext(os.path.basename(recruit_s.mg_file))[0]
    mg_file = (mg_id, recruit_s.mg_file)
    mg_sub_file = s_utils.build_subcontigs(
        'Metagenomes',
        [recruit_s.mg_file],
        recruit_s.save_path,
        recruit_s.max_contig_len,
        recruit_s.overlap_len,
        recruit_s.min_len,
    )[0]
    return mg_file, mg_sub_file


def _build_trusted_assets(recruit_s, mg_file):
    """Build trusted contig assets and minhash signatures when configured."""
    if not recruit_s.trust_path:
        return tuple(), False
    import scarab.minhash_recruiter as mhr
    import scarab.utilities as s_utils

    tc_list = s_utils.get_SAGs(recruit_s.trust_path)
    trust_files = tuple((os.path.splitext(os.path.basename(x))[0], x) for x in tc_list)
    trust_subs = s_utils.build_subcontigs(
        'SAGs',
        tc_list,
        recruit_s.save_path,
        recruit_s.max_contig_len,
        recruit_s.overlap_len,
        recruit_s.min_len,
    )
    minhash_df_dict = mhr.run_minhash_recruiter(
        recruit_s.save_path,
        recruit_s.save_path,
        trust_subs,
        mg_file,
        recruit_s.nthreads,
        recruit_s.min_len,
        recruit_s.kmer_size,
    )
    return trust_files, minhash_df_dict


def _build_abundance_tables(recruit_s, mg_sub_file):
    """Generate abundance table outputs for clustering."""
    import scarab.abundance_recruiter as abr

    return abr.runAbundRecruiter(
        recruit_s.save_path,
        recruit_s.save_path,
        mg_sub_file,
        recruit_s.mg_raw_file_list,
        recruit_s.pacbio,
        recruit_s.nthreads,
    )


def _build_tetranuc_table(recruit_s, mg_sub_file):
    """Generate the tetranucleotide frequency table."""
    import scarab.tetranuc_recruiter as tra

    return tra.run_tetra_recruiter(recruit_s.save_path, mg_sub_file)


def _set_clustering_params(recruit_s, abund_raw_file):
    """Resolve clustering parameters for the recruit run."""
    import scarab.utilities as s_utils

    return s_utils.set_clust_params(
        recruit_s.denovo_min_clust,
        recruit_s.denovo_min_samp,
        recruit_s.anchor_min_clust,
        recruit_s.anchor_min_samp,
        recruit_s.nu,
        recruit_s.gamma,
        recruit_s.vr,
        recruit_s.r,
        recruit_s.s,
        recruit_s.vs,
        recruit_s.a,
        abund_raw_file,
        recruit_s.save_path,
    )


def _run_clustering(recruit_s, mg_sub_file, save_dirs_dict, abund_scale_file, tetra_file, minhash_df_dict):
    """Run clustering for the recruit pipeline."""
    import scarab.clusterer as clst

    mg_id = mg_sub_file[0]
    return clst.runClusterer(
        mg_id,
        save_dirs_dict[recruit_s.set],
        save_dirs_dict[recruit_s.set],
        abund_scale_file,
        tetra_file,
        minhash_df_dict,
        recruit_s.params_dict['d_min_clust'],
        recruit_s.params_dict['d_min_samp'],
        recruit_s.params_dict['a_min_clust'],
        recruit_s.params_dict['a_min_samp'],
        recruit_s.params_dict['nu'],
        recruit_s.params_dict['gamma'],
        recruit_s.jaccard,
        recruit_s.nthreads,
    )


def info(sys_args):
    """Write version and dependency information for SCARAB."""
    import scarab.s_args as s_args
    import scarab.logger as s_log
    from scarab.__init__ import version

    parser = s_args.ScarabArgumentParser(description='Return package and executable information.')
    args = parser.parse_args(sys_args)
    s_log.prep_logging()

    logger.info('SCARAB version %s.', version)

    py_deps = _python_dependency_versions()
    logger.info(
        'Python package dependency versions:\n\t'
        + '\n\t'.join([k + ': ' + v for k, v in py_deps.items()])
    )

    if args.verbose:
        pass

    return


def recruit(sys_args):
    """Recruit environmental reads to reference contigs."""
    import scarab.logger as s_log
    import scarab.compile_recruits as com
    import scarab.utilities as s_utils

    args = _parse_recruit_args(sys_args)
    recruit_s = _configure_recruit_base(args)
    recruit_s.mode = _select_mode(recruit_s)
    s_log.prep_logging(os.path.join(recruit_s.save_path, 'SCARAB_log.txt'), verbosity=args.verbose)

    mg_file, mg_sub_file = _build_metagenome_subcontigs(recruit_s)
    trust_files, minhash_df_dict = _build_trusted_assets(recruit_s, mg_file)
    abund_scale_file, abund_raw_file = _build_abundance_tables(recruit_s, mg_sub_file)
    tetra_file = _build_tetranuc_table(recruit_s, mg_sub_file)

    recruit_s.mode, recruit_s.set, recruit_s.params_dict = _set_clustering_params(recruit_s, abund_raw_file)
    save_dirs_dict = s_utils.check_out_dirs(recruit_s.save_path, recruit_s.mode, recruit_s.set)
    clusters = _run_clustering(
        recruit_s,
        mg_sub_file,
        save_dirs_dict,
        abund_scale_file,
        tetra_file,
        minhash_df_dict,
    )
    com.run_combine_recruits(
        save_dirs_dict,
        recruit_s.mg_file,
        clusters,
        trust_files,
        recruit_s.set,
        recruit_s.nthreads,
    )

    return
