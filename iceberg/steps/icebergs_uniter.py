"""
This Python module is used to run the UNITE_CLOSE_ICEBERGS_SITES step on treatment and control experiments.

The step unite close icebergs sites on the genome to one large icebergs sites
"""

import settings
from utils import generators, csv_to_df, logs, validation

from collections import Counter
import pandas as pd
import logging
from pathlib import Path


STEP_NAME = 'UNITE_CLOSE_ICEBERGS_SITES'
STEP_OUT_FOLDER = 'UNITE_CLOSE_ICEBERGS_SITES_RESULTS'
STEP_LOG_FILE = 'UNITE_CLOSE_ICEBERGS_SITES.log'


class IcebergsUniter(object):

    def __init__(self,
                 distance: int,
                 in_folder: Path,
                 out_folder: Path) -> None:
        """
        Initiate the IcebergsUniter object with the thresholds for the unification process.

        :param distance: The distance (in bps) for unite nearby icebergs sites.
        :param in_folder: The folder where the input csv file is stored.
        :param out_folder: The folder where the output csv file will be saved.
        """
        self.distance = distance
        self.in_folder = in_folder
        self.out_folder = out_folder

    def unite_icebergs_sites(self,
                             iceberg: pd.Series,
                             iceberg_to_add: pd.Series) -> pd.Series:
        """
        Unites close iceberg by modified iceberg fields.

        :param iceberg: The first iceberg to unite.
        :param iceberg_to_add: The second iceberg to unite.
        """
        mock_columns_end_symbol = 'm'
        treatment_columns_end_symbol = 't'
        united_sites = self.unite_experiment_icebergs_sites(iceberg, iceberg_to_add,
                                                            mock_columns_end_symbol)

        united_sites = self.unite_experiment_icebergs_sites(united_sites, iceberg_to_add,
                                                            treatment_columns_end_symbol)
        # logging.info(f"{iceberg}")
        # logging.info(f"{iceberg_to_add}")

        # logging.info(f"start t: {united_sites[f'start index {treatment_columns_end_symbol}']}")
        # logging.info(f"start control: {united_sites[f'start index {mock_columns_end_symbol}']}")
        #
        # logging.info(f"end t: {united_sites[f'end index {treatment_columns_end_symbol}']}")
        # logging.info(f"end control: {united_sites[f'end index {mock_columns_end_symbol}']}")
        if united_sites['chromosome t'] == '*':
            united_sites[f'start index'] = united_sites[f'start index {mock_columns_end_symbol}']
            united_sites[f'end index'] = united_sites[f'end index {mock_columns_end_symbol}']

        elif united_sites['chromosome m'] == '*':
            united_sites[f'start index'] = united_sites[f'start index {treatment_columns_end_symbol}']
            united_sites[f'end index'] = united_sites[f'end index {treatment_columns_end_symbol}']

        else:
            united_sites[f'start index'] = min(united_sites[f'start index {treatment_columns_end_symbol}'],
                                               united_sites[f'start index {mock_columns_end_symbol}'])
            united_sites[f'end index'] = max(united_sites[f'end index {treatment_columns_end_symbol}'],
                                             united_sites[f'end index {mock_columns_end_symbol}'])
        return united_sites

    def unite_experiment_icebergs_sites(self,
                                        iceberg: pd.Series,
                                        iceberg_to_add: pd.Series,
                                        exp_suffix: str) -> pd.Series:
        """
        Unites close iceberg by modified iceberg fields.

        :param iceberg: The first iceberg to unite.
        :param iceberg_to_add: The second iceberg to unite.
        :param exp_suffix: 't' or 'm'.
        """
        if iceberg[f'chromosome {exp_suffix}'] == '*' or iceberg_to_add[f'chromosome {exp_suffix}'] == '*':
            if iceberg[f'chromosome {exp_suffix}'] == '*':
                experiment_indexes = iceberg.index.str.endswith(exp_suffix)
                iceberg[experiment_indexes] = iceberg_to_add[experiment_indexes]
        else:
            gap_size = iceberg_to_add[f'start index {exp_suffix}'] - iceberg[f'end index {exp_suffix}'] - 1
            iceberg[f'end index {exp_suffix}'] = iceberg_to_add[f'end index {exp_suffix}']
            iceberg[f'iceberg length {exp_suffix}'] += gap_size + iceberg_to_add[f'iceberg length {exp_suffix}']

            iceberg[f'mapq {exp_suffix}'] = (iceberg[f'mapq {exp_suffix}'] * iceberg[f'read count {exp_suffix}'] +
                                             iceberg_to_add[f'mapq {exp_suffix}'] *
                                             iceberg_to_add[f'read count {exp_suffix}'])

            iceberg[f'mapq {exp_suffix}'] /= (iceberg[f'read count {exp_suffix}'] +
                                              iceberg_to_add[f'read count {exp_suffix}'])

            iceberg[f'mapq {exp_suffix}'] = round(iceberg[f'mapq {exp_suffix}'])

            iceberg[f'read count {exp_suffix}'] += iceberg_to_add[f'read count {exp_suffix}']
            iceberg[f'cut-position-candidates {exp_suffix}'].update(iceberg_to_add[f'cut-position-candidates {exp_suffix}'])
            iceberg[f'orientations {exp_suffix}'] = dict(Counter(iceberg[f'orientations {exp_suffix}']) +
                                                         Counter(iceberg_to_add[f'orientations {exp_suffix}']))

            iceberg[f'iceberg depth {exp_suffix}'] = max(iceberg[f'iceberg depth {exp_suffix}'],
                                                         iceberg_to_add[f'iceberg depth {exp_suffix}'])

        return iceberg

    def unite_chromosome_close_icebergs_sites(self,
                                              icebergs_df: pd.DataFrame) -> pd.DataFrame:
        """
        Unite close icebergs sites of a chromosome based on sites locus (after merging treatment and control based
        on icebergs locus).

        The function handles the cases where there are icebergs sites within the chromosome only in one of the
        experiments (treatment or control).

        :param icebergs_df: A dataframe of sorted icebergs from specific chromosome of treatment experiment.

        :return: A dataframe which contains both experiments icebergs joint by chromosome index and scored for the
                 classification process.
        """
        i = 0
        united_df = []
        while i < icebergs_df.shape[0]:
            united_sites = icebergs_df.loc[i, :].copy()
            while i + 1 < icebergs_df.shape[0]:
                if icebergs_df.loc[i + 1, 'start index'] <= united_sites['end index'] + self.distance:
                    united_sites = self.unite_icebergs_sites(united_sites.copy(), icebergs_df.loc[i + 1, :])
                    # united_sites = self.unite_icebergs_sites(united_sites.copy(), icebergs_df.loc[i + 1, :])
                    i += 1
                else:
                    break
            united_df.append(united_sites)
            i += 1

        united_df = pd.DataFrame(united_df)
        return united_df

    def unite_close_icebergs_sites(self,
                                   merged_csv_file_name: str) -> (Path, Path, Path):
        """
        Unite close icebergs sites based on sites locus (after merging treatment and control based
        on icebergs locus).

        :param merged_csv_file_name: The merged icebergs csv file name

        :return: The three csv files path (CRISPR activity csv, spontaneous breaks csv and noise csv).
        """
        res_file_name = merged_csv_file_name.replace('.icebergs.csv', '.united.icebergs.csv')
        icebergs_df = csv_to_df.get_merged_experiments_df(self.in_folder / merged_csv_file_name)
        united_df = pd.DataFrame(columns=[col for col in icebergs_df.columns])
        pbar = logs.get_progress_bar(step_name=STEP_NAME,
                                     total=icebergs_df.shape[0],
                                     out_of='icebergs-site')

        chromosomes = generators.get_icebergs_by_chromosomes_merged(icebergs_df)
        # iterate over all chromosomes names.
        for chromosome in chromosomes:
            united_df = pd.concat([united_df, self.unite_chromosome_close_icebergs_sites(icebergs_df=chromosome['df'])])
            pbar.update(chromosome['df'].shape[0])
        united_df.to_csv(path_or_buf=self.out_folder / res_file_name, index=False)
        return res_file_name


def validate_arguments(args: dict) -> None:
    """
    Validates the step arguments from the given arguments dict.

    :param args: The arguments dict (see main function description).
    """
    try:
        validation.validate_path_existence(args['curr_in_dir'])
        if args['tx_and_control_merged_icebergs'] != '':
            validation.validate_path_existence(Path(args['curr_in_dir'], args['tx_and_control_merged_icebergs']))
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
    Run the UNITE_CLOSE_ICEBERGS_SITES step on treatment and control.

    :param args: The pipeline arguments dictionary - dictionary of arguments (as keys) and their values.
    """
    try:
        iu = IcebergsUniter(args['HYPERPARAMATERS']['MAX_ICEBERG_DISTANCE'],
                            args['curr_in_dir'],
                            args['curr_out_dir'])

        args['tx_and_control_merged_icebergs'] = iu.unite_close_icebergs_sites(args['tx_and_control_merged_icebergs'])

        args['curr_in_dir'] = validation.set_curr_in_dir(args['OUTPUT_FOLDER_PATH'], STEP_OUT_FOLDER)
        logging.info(f'{STEP_NAME} step done.')
    except Exception as e:
        logging.error(f'Failed to run {STEP_NAME} step, please read the following for more information.',
                      exc_info=False)
        # sys.exit()
        raise e
