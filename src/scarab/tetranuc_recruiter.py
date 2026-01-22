import logging
import warnings
from os.path import isfile
from os.path import join as o_join

import scarab.utilities as s_utils

warnings.simplefilter(action='ignore', category=FutureWarning)

logger = logging.getLogger(__name__)


def run_tetra_recruiter(tra_path, mg_sub_file):
    logger.info('Starting tetranucleotide data transformation.')
    mg_id = mg_sub_file[0]
    if isfile(o_join(tra_path, mg_id + '.tetras.tsv')):
        logger.info('Loading tetramer Hz matrix for %s.', mg_id)
        mg_tetra_file = o_join(tra_path, mg_id + '.tetras.tsv')
    else:
        logger.info('Calculating tetramer Hz matrix for %s.', mg_id)
        mg_subcontigs = s_utils.get_seqs(mg_sub_file[1])
        mg_tetra_df = s_utils.tetra_cnt(mg_subcontigs)
        mg_tetra_df.to_csv(o_join(tra_path, mg_id + '.tetras.tsv'),
                           sep='\t'
                           )
        mg_tetra_file = o_join(tra_path, mg_id + '.tetras.tsv')

    return mg_tetra_file
