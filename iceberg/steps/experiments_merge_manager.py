"""
This Python module is used to run the MERGE_TREATMENT_AND_CONTROL_ICEBERGS_BY_LOCUS step on treatment and control
experiments.

Merge the csv files of treatment and control to one csv files where iceberg sites from different experiments
with the same chromosome and overlapping indexes are in the same row in the csv, experiment suffix is added to
all iceberg field, t for treatment and m for control/mock.

If in specific location there is iceberg only in one experiment, then empty iceberg is added to the row of
the present iceberg with the needed suffix.
"""

import settings
from utils import generators, csv_to_df, logs, validation

import numpy as np
import pandas as pd
import logging
from pathlib import Path

STEP_NAME = 'MERGE_TREATMENT_AND_CONTROL_ICEBERGS_BY_LOCUS'
STEP_OUT_FOLDER = 'MERGE_TREATMENT_AND_CONTROL_ICEBERGS_BY_LOCUS_RESULTS'
STEP_LOG_FILE = 'MERGE_TREATMENT_AND_CONTROL_ICEBERGS_BY_LOCUS.log'


class ExperimentsMergeManager(object):

    def __init__(self,
                 # distance: int,
                 in_folder: Path,
                 out_folder: Path
                 ) -> None:
        """
        Initiate the ExperimentsMergeManager object with the distance thresholds for the
        treatment and control merging process.

        # :param distance: The distance (in bps) for unite close icebergs.
        :param in_folder: The folder where the input csv files are stored.
        :param out_folder: The folder where the output csv file will be saved.
        """
        # self.distance = distance
        self.in_folder = in_folder
        self.out_folder = out_folder

    def merge_sites(self,
                    t_iceberg: pd.Series,
                    m_iceberg: pd.Series) -> pd.Series:
        """
        Merge two icebergs from different experiments (treatment and control).
        If one of the icebergs is None, an empty iceberg will be merged instead.

        :param t_iceberg: Iceberg from treatment (None if there is no treatment iceberg in the site)
        :param m_iceberg: Iceberg from control (None if there is no control iceberg in the site)

        :return: icebergs site.
        """
        merged_sites = pd.Series()
        if m_iceberg is None:
            merged_sites['chromosome'] = t_iceberg['chromosome']
            merged_sites['start index'] = t_iceberg['start index']
            merged_sites['end index'] = t_iceberg['end index']
            merged_sites['iceberg length'] = int(merged_sites['end index'] - merged_sites['start index'])

            merged_sites = merged_sites.append(pd.Series(data=np.array(t_iceberg),
                                                         index=[i + ' t' for i in t_iceberg.index]))

            merged_sites = merged_sites.append(pd.Series(data=[0
                                                               if k in ['read count', 'iceberg length', 'iceberg depth',
                                                                        'start index', 'end index']
                                                               else -1
                                                               if k in ['mapq']
                                                               else {}
                                                               if k in ['cut-position-candidates', 'orientations']
                                                               else '*'
                                                               for k in t_iceberg.index],
                                                         index=[i + ' m' for i in t_iceberg.index]))
        elif t_iceberg is None:
            merged_sites['chromosome'] = m_iceberg['chromosome']
            merged_sites['start index'] = m_iceberg['start index']
            merged_sites['end index'] = m_iceberg['end index']
            merged_sites['iceberg length'] = int(merged_sites['end index'] - merged_sites['start index'])

            merged_sites = merged_sites.append(pd.Series(data=[0
                                                               if k in ['read count', 'iceberg length', 'iceberg depth',
                                                                        'start index', 'end index']
                                                               else -1
                                                               if k in ['mapq']
                                                               else {}
                                                               if k in ['cut-position-candidates', 'orientations']
                                                               else '*'
                                                               for k in m_iceberg.index],
                                                         index=[i + ' t' for i in m_iceberg.index]))

            merged_sites = merged_sites.append(pd.Series(data=np.array(m_iceberg),
                                                         index=[i + ' m' for i in m_iceberg.index]))
        else:
            merged_sites['chromosome'] = t_iceberg['chromosome']
            merged_sites['start index'] = int(min(t_iceberg['start index'], m_iceberg['start index']))
            merged_sites['end index'] = int(max(t_iceberg['end index'], m_iceberg['end index']))
            merged_sites['iceberg length'] = int(merged_sites['end index'] - merged_sites['start index'])

            merged_sites = merged_sites.append(pd.Series(data=np.array(t_iceberg),
                                                         index=[i + ' t' for i in t_iceberg.index]))

            merged_sites = merged_sites.append(pd.Series(data=np.array(m_iceberg),
                                                         index=[i + ' m' for i in m_iceberg.index]))

        return merged_sites

    def merge_treatment_and_control_chromosome_icebergs_by_locus(self,
                                                                 t_icebergs_df: pd.DataFrame = None,
                                                                 m_icebergs_df: pd.DataFrame = None) -> pd.DataFrame:
        """
        Merge close treatment and control clusters/icebergs sites based on location within a chromosome.

        :param t_icebergs_df: A dataframe of sorted icebergs from specific chromosome of treatment experiment.
        :param m_icebergs_df: A dataframe of sorted icebergs from specific chromosome of control experiment.

        :return: A dataframe which contains both experiments icebergs joint by chromosome index and scored for the
                 classification process.
        """
        i = 0
        j = 0
        merged_df = []
        if t_icebergs_df is not None and m_icebergs_df is not None:
            while i < t_icebergs_df.shape[0] and j < m_icebergs_df.shape[0]:
                t_iceberg = t_icebergs_df.iloc[i, :].copy()
                m_iceberg = m_icebergs_df.iloc[j, :].copy()
                if m_iceberg['start index'] > t_iceberg['end index']:  # + self.distance:
                    merged_df.append(self.merge_sites(t_iceberg, None))
                    i += 1
                elif t_iceberg['start index'] > m_iceberg['end index']:  # + self.distance:
                    merged_df.append(self.merge_sites(None, m_iceberg))
                    j += 1
                else:
                    merged_df.append(self.merge_sites(t_iceberg.copy(), m_iceberg.copy()))
                    i += 1
                    j += 1

        if t_icebergs_df is not None:
            while i < t_icebergs_df.shape[0]:
                t_iceberg = t_icebergs_df.iloc[i, :]
                merged_df.append(self.merge_sites(t_iceberg, None))
                i += 1
        if m_icebergs_df is not None:
            while j < m_icebergs_df.shape[0]:
                m_iceberg = m_icebergs_df.iloc[j, :]
                merged_df.append(self.merge_sites(None, m_iceberg))
                j += 1
        merged_df = pd.DataFrame(merged_df)
        return merged_df

    def merge_treatment_and_control_icebergs_by_locus(self,
                                                      csv1_file_name: str,
                                                      csv2_file_name: str) -> (Path, Path, Path):
        """
        Merge close icebergs from treatment and control experiments based on locus (chromosome and location).

        :param csv1_file_name: The control icebergs csv file name (stored in in_dir, the output of read_uniter.py).
        :param csv2_file_name: The treatment icebergs csv file name (stored in in_dir, the output of read_uniter.py).

        :return: The three csv files path (CRISPR activity csv, spontaneous breaks csv and noise csv).
        """
        merged_file_name = csv1_file_name.replace('.icebergs.csv', '') + '-' + \
            csv2_file_name.replace('.icebergs.csv', '') + '.icebergs.csv'

        t_icebergs_df = csv_to_df.get_single_experiment_df(self.in_folder / csv1_file_name)
        m_icebergs_df = csv_to_df.get_single_experiment_df(self.in_folder / csv2_file_name)

        pbar = logs.get_progress_bar(step_name=STEP_NAME,
                                     total=t_icebergs_df.loc[:, 'read count'].map(lambda x: float(x)).sum() +
                                     m_icebergs_df.loc[:, 'read count'].map(lambda x: float(x)).sum(),
                                     out_of="read")

        merged_df = pd.DataFrame(columns=[col + ' t' for col in t_icebergs_df.columns] +
                                         [col + ' m' for col in m_icebergs_df.columns])

        chromosomes = generators.get_icebergs_by_chromosomes(t_icebergs_df,
                                                             m_icebergs_df)
        # iterate over all chromosomes names.
        merge = self.merge_treatment_and_control_chromosome_icebergs_by_locus
        for chromosome in chromosomes:

            if 'm df' not in chromosome:
                merged_chromosome_df = merge(t_icebergs_df=chromosome['t df'])
                pbar.update(chromosome['t df']['read count'].map(lambda x: float(x)).sum())
            elif 't df' not in chromosome:
                merged_chromosome_df = merge(m_icebergs_df=chromosome['m df'])
                pbar.update(chromosome['m df']['read count'].map(lambda x: float(x)).sum())

            else:
                merged_chromosome_df = merge(t_icebergs_df=chromosome['t df'],
                                             m_icebergs_df=chromosome['m df'])

                pbar.update(chromosome['t df']['read count'].map(lambda x: float(x)).sum() +
                            chromosome['m df']['read count'].map(lambda x: float(x)).sum())

            merged_df = pd.concat([merged_df, merged_chromosome_df])

        merged_df['start index'] = merged_df['start index'].astype(int)
        merged_df['end index'] = merged_df['end index'].astype(int)
        merged_df['iceberg length'] = merged_df['iceberg length'].astype(int)
        merged_df['mapq m'] = merged_df['mapq m'].astype(int)
        merged_df['mapq t'] = merged_df['mapq t'].astype(int)

        merged_df.to_csv(path_or_buf=self.out_folder / merged_file_name, index=False)

        return merged_file_name


def validate_arguments(args: dict) -> None:
    """
    Validates the step arguments from the given arguments dict.

    :param args: The arguments dict (see main function description).
    """
    try:
        validation.validate_path_existence(Path(args['curr_in_dir'], args['csv1']))
        validation.validate_path_existence(Path(args['curr_in_dir'], args['csv2']))
    except ValueError as e:
        logging.error(f'not all necessary arguments for {STEP_NAME} step where giving.')
        raise e


def prepare_step(args: dict) -> None:
    """
    Prepares the step by - setting the logging configurations, setting the output directory and
    validate the step arguments.

    :param args: A dict of arguments (as keys) and their values (see main function description)
    """
    logs.set_logging_config(STEP_NAME, args['OUTPUT_FOLDER_PATH'] / settings.LOGS_FOLDER_RELATIVE_PATH / STEP_LOG_FILE)
    args['curr_out_dir'] = validation.set_curr_out_dir(args['OUTPUT_FOLDER_PATH'], STEP_OUT_FOLDER)
    validate_arguments(args)


def main(args: dict) -> None:
    """
    Run the MERGE_TREATMENT_AND_CONTROL_ICEBERGS_BY_LOCUS step on treatment and control.

    :param args: The pipeline arguments dictionary - dictionary of arguments (as keys) and their values.
    """
    try:
        emm = ExperimentsMergeManager(  # args['HYPERPARAMATERS']['MAX_ICEBERG_DISTANCE'],
                                      in_folder=args['curr_in_dir'],
                                      out_folder=args['curr_out_dir'])

        args['tx_and_control_merged_icebergs'] = emm.merge_treatment_and_control_icebergs_by_locus(args['csv1'],
                                                                                                   args['csv2'])

        args['curr_in_dir'] = validation.set_curr_in_dir(args['OUTPUT_FOLDER_PATH'], STEP_OUT_FOLDER)
        logging.info(f'{STEP_NAME} step done.')
    except Exception as e:
        logging.error(f'Failed to run {STEP_NAME} step, please read the following for more information.',
                      exc_info=False)
        # sys.exit()
        raise e
