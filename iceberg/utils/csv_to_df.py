"""
This Python module is used to read icebergs sites csv file and convert it to pandas dataframe properly.
"""

import pandas as pd
import logging


def get_single_experiment_df(icebergs_df_path):
    """
    Reads the experiment (treatment or control) icebergs (reads clustered by their chromosome index)
    from csv file and convert them to dataframe properly.

    :param icebergs_df_path: The absolute path to the icebergs csv file.

    :return: An experiment (control or treatment) icebergs dataframe.
    """
    try:
        converters = {'cut-position-candidates': eval,
                      'orientations': eval}

        dtype = {'chromosome': 'str',
                 'start index': 'int',
                 'end index': 'int',
                 'mapq': 'int',
                 'read count': 'int',
                 'iceberg depth': 'int',
                 'iceberg length': 'int'}

        icebergs_df = pd.read_csv(icebergs_df_path, converters=converters, dtype=dtype)
        return icebergs_df
    except Exception as e:
        logging.error(f'Failed to read csv file properly :{icebergs_df_path}', exc_info=False)
        raise e


def get_merged_experiments_df(icebergs_df_path,
                              csv_contains_alignments_fields=False):
    """
    Reads the experiment (treatment and control) icebergs (reads clustered by their chromosome index)
    from csv file and convert them to dataframe properly.

    :param icebergs_df_path: The absolute path to the icebergs csv file.
    :param csv_contains_alignments_fields: True if csv file contains GuideRNA alignments fields, else, False.
    :return: An experiment (control or treatment) icebergs dataframe.
    """
    try:
        converters = {'cut-position-candidates m': eval,
                      'cut-position-candidates t': eval,
                      'orientations m': eval,
                      'orientations t': eval}

        dtype = {'chromosome': 'str',
                 'start index': 'int',
                 'end index': 'int',
                 'iceberg length': 'int',
                 'chromosome t': 'str',
                 'chromosome m': 'str',
                 'start index t': 'int',
                 'start index m': 'int',
                 'end index t': 'int',
                 'end index m': 'int',
                 'mapq t': 'int',
                 'mapq m': 'int',
                 'read count t': 'int',
                 'read count m': 'int',
                 'iceberg depth t': 'int',
                 'iceberg depth m': 'int',
                 'iceberg length t': 'int',
                 'iceberg length m': 'int'}
        if csv_contains_alignments_fields:
            dtype.update({'cut position': 'int',
                          'strand': 'str',
                          'peak t': 'int',
                          'gRNA alignment': 'str',
                          'short site alignment': 'str',
                          'Hamming distance': 'int',
                          'Levenshtein distances': 'int'})

        icebergs_df = pd.read_csv(icebergs_df_path, converters=converters, dtype=dtype)
        # icebergs_df = pd.read_csv(icebergs_df_path,
        #                           converters={'cut-position-candidates m': eval,
        #                                       'cut-position-candidates t': eval,
        #                                       'orientations m': eval,
        #                                       'orientations t': eval},
        #                           dtype={'chromosome': 'str',
        #                                  'start index': 'int',
        #                                  'end index': 'int',
        #                                  'iceberg length': 'int',
        #                                  'chromosome t': 'str',
        #                                  'chromosome m': 'str',
        #                                  'start index t': 'int',
        #                                  'start index m': 'int',
        #                                  'end index t': 'int',
        #                                  'end index m': 'int',
        #                                  'mapq t': 'int',
        #                                  'mapq m': 'int',
        #                                  'read count t': 'int',
        #                                  'read count m': 'int',
        #                                  'iceberg depth t': 'int',
        #                                  'iceberg depth m': 'int',
        #                                  'iceberg length t': 'int',
        #                                  'iceberg length m': 'int'})
        #
        # if csv_contains_alignments_fields:
        #     icebergs_df = pd.read_csv(icebergs_df_path,
        #                               converters={'cut-position-candidates m': eval,
        #                                           'cut-position-candidates t': eval,
        #                                           'orientations m': eval,
        #                                           'orientations t': eval},
        #                               dtype={'chromosome': 'str',
        #                                      'start index': 'int',
        #                                      'end index': 'int',
        #                                      'iceberg length': 'int',
        #                                      'chromosome t': 'str',
        #                                      'chromosome m': 'str',
        #                                      'start index t': 'int',
        #                                      'start index m': 'int',
        #                                      'end index t': 'int',
        #                                      'end index m': 'int',
        #                                      'mapq t': 'int',
        #                                      'mapq m': 'int',
        #                                      'read count t': 'int',
        #                                      'read count m': 'int',
        #                                      'iceberg depth t': 'int',
        #                                      'iceberg depth m': 'int',
        #                                      'iceberg length t': 'int',
        #                                      'iceberg length m': 'int',
        #                                      'cut position': 'int',
        #                                      'strand': 'str',
        #                                      'peak t': 'int',
        #                                      'gRNA alignment': 'str',
        #                                      'short site alignment': 'str',
        #                                      'Hamming distance': 'int',
        #                                      'Levenshtein distances': 'int'})

        return icebergs_df
    except Exception as e:
        logging.error(f'Failed to read csv file properly:{icebergs_df_path}', exc_info=False)
        raise e