"""
This Python module is used to run the CALCULATE_ICEBERGS_SITES_PROFILE step on treatment and control experiments.

Calculate the iceberg sites profile from the cut position candidates dictionary for each site in the
CRISPR activity sites.

For each site two html files are generated - one for the html table, which contains the amount of each paired end reads
from each experiment library (forward and reverse) for each cut position candidate, and one for the html
histogram which visualize this table.

the html file names are chromosome-start_index-end_index-table.html and chromosome-start_index-end_index-histogram.html

"""

import settings
from utils import generators, csv_to_df, logs, validation

import pandas as pd
import logging
from pathlib import Path
import plotly.express as px

STEP_NAME = 'CALCULATE_ICEBERGS_SITES_PROFILE'
STEP_OUT_FOLDER = 'CALCULATE_ICEBERGS_SITES_PROFILE_RESULTS'
STEP_LOG_FILE = 'CALCULATE_ICEBERGS_SITES_PROFILE.log'


class SitesProfileCalculator(object):

    def __init__(self,
                 in_folder: Path,
                 out_folder: Path) -> None:
        """
        Initiate the sitesProfileCalculator object with the thresholds for the unification process.

        :param in_folder: The folder where the input csv file is stored.
        :param out_folder: The folder where the output files directory will be saved.
        """
        self.in_folder = in_folder
        self.out_folder = out_folder

    def generate_site_profile_histogram_df(self,
                                           iceberg_cut_position_candidates: dict):
        """
        Generate the site profile dataframe.

        :param iceberg_cut_position_candidates: The iceberg cut position candidates data.

        :return: The site profile data in dataframe.
        """
        if iceberg_cut_position_candidates == {}:
            return None

        df = []
        for cut_position_candidate in range(min(iceberg_cut_position_candidates.keys()),
                                            max(iceberg_cut_position_candidates.keys())+1):
            if cut_position_candidate in iceberg_cut_position_candidates:
                if 'forward-primer' in iceberg_cut_position_candidates[cut_position_candidate]:
                    hist_data = {'cut candidate position': cut_position_candidate,
                                 'library': 'forward-primer',
                                 'count': iceberg_cut_position_candidates[cut_position_candidate]['forward-primer']}

                    df.append(hist_data)
                if 'reverse-primer' in iceberg_cut_position_candidates[cut_position_candidate]:
                    hist_data = {'cut candidate position': cut_position_candidate,
                                 'library': 'reverse-primer',
                                 'count': iceberg_cut_position_candidates[cut_position_candidate]['reverse-primer']}

                    df.append(hist_data)
            else:
                hist_data = {'cut candidate position': cut_position_candidate,
                             'library': 'forward-primer',
                             'count': 0}

                df.append(hist_data)

                hist_data = {'cut candidate position': cut_position_candidate,
                             'library': 'reverse-primer',
                             'count': 0}

                df.append(hist_data)

        if len(df) > 0:
            df = pd.DataFrame(df)
            df['cut candidate position'] = df['cut candidate position'].astype(int)
            df = df.sort_values(by=['cut candidate position'])
        else:
            df = None

        return df

    def generate_site_profile_histogram_html_file(self,
                                                  iceberg_cut_position_candidates: dict,
                                                  out_folder: Path,
                                                  out_file_name: str):
        """
        Generate the site profile histogram as html file.

        :param iceberg_cut_position_candidates: The iceberg cut position candidates data.
        :param out_folder: The folder where the output files will be saved.
        :param out_file_name: The output file name.
        """
        histogram_df = self.generate_site_profile_histogram_df(iceberg_cut_position_candidates)
        if histogram_df is not None:
            fig = px.histogram(histogram_df,
                               x=[str(x) for x in histogram_df["cut candidate position"]],
                               y="count",
                               color='library',
                               barmode='group',
                               nbins=len(histogram_df),
                               # marginal="box",  # or violin, rug
                               height=400)
            fig.write_html(str(out_folder / out_file_name))

    def generate_site_profile_table_html_file(self,
                                              iceberg_cut_position_candidates: dict,
                                              out_folder: Path,
                                              out_file_name: str):
        """
        Generate the site profile table as html file.

        :param iceberg_cut_position_candidates: The iceberg cut position candidates data.
        :param out_folder: The folder where the output files directory will be saved.
        :param out_file_name: The output file name.



        """
        if iceberg_cut_position_candidates == {}:
            return

        highest_candidate = max(iceberg_cut_position_candidates,
                                key=lambda candidate: sum(iceberg_cut_position_candidates[candidate][primer]
                                                          for primer in iceberg_cut_position_candidates[candidate]))

        cut_position_candidates = dict(sorted(iceberg_cut_position_candidates.items()))

        df = pd.DataFrame(cut_position_candidates).T

        df["distance from highest cut position candidate"] = [cut_position_candidate - highest_candidate for
                                                              cut_position_candidate in cut_position_candidates]

        df = df.fillna(0)

        # write html to file
        text_file = open(out_folder / out_file_name, "w")
        text_file.write(df.to_html(classes='table table-stripped'))
        text_file.close()

    def generate_chromosome_sites_profile_html_files(self,
                                                     icebergs_df: pd.DataFrame,
                                                     out_folder: Path):
        """
        Generate chromosome sites profile html files (histogram and table)

        :param icebergs_df: A dataframe of sorted icebergs from specific chromosome of treatment experiment.
        :param out_folder: The folder where the output files will be saved.

        :return: A dataframe which contains both experiments icebergs joint by chromosome index and scored for the
                 classification process.
        """
        for i in icebergs_df.index:
            cut_position_candidates = dict([(int(k), v)
                                            for k, v in icebergs_df.loc[i, 'cut-position-candidates t'].items()])
            chromosome = icebergs_df.loc[i, "chromosome"]
            start_index = int(icebergs_df.loc[i, "start index"])
            end_index = int(icebergs_df.loc[i, "end index"])
            out_file_name = f'{chromosome}-{start_index}-{end_index}'
            self.generate_site_profile_histogram_html_file(cut_position_candidates,
                                                           out_folder,
                                                           f'{out_file_name}-histogram.html')

            self.generate_site_profile_table_html_file(cut_position_candidates,
                                                       out_folder,
                                                       f'{out_file_name}-table.html')

    def calculate_sites_profile(self,
                                merged_csv_file_name: str) -> (Path, Path, Path):
        """
        Calculate the iceberg sites profile from the cut position candidates dictionary.
        :param merged_csv_file_name: The treatment icebergs csv file name.

        :return: The three csv files path (CRISPR activity csv, spontaneous breaks csv and noise csv).
        """
        icebergs_df = csv_to_df.get_merged_experiments_df(self.in_folder / merged_csv_file_name)

        pbar = logs.get_progress_bar(step_name=STEP_NAME,
                                total=icebergs_df.shape[0],
                                out_of="icebergs-site")

        chromosomes = generators.get_icebergs_by_chromosomes_merged(icebergs_df)
        # iterate over all chromosomes names.
        _ = validation.set_curr_out_dir(self.out_folder, 'sites-profiles')
        for chromosome in chromosomes:
            self.generate_chromosome_sites_profile_html_files(icebergs_df=chromosome['df'],
                                                              out_folder=self.out_folder / 'sites-profiles')
            pbar.update(chromosome['df'].shape[0])

        return self.out_folder / 'sites-profiles'


def validate_arguments(args: dict) -> None:
    """
    Validates the step arguments from the given arguments dict.

    :param args: The arguments dict (see main function description).
    """
    try:
        validation.validate_path_existence(args['curr_in_dir'])
        if args['tx_and_control_merged_icebergs'] != '':
            validation.validate_path_existence(
                Path(args['curr_in_dir'], args['tx_and_control_merged_icebergs']))
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
    Run the CALCULATE_ICEBERGS_SITES_PROFILE step on treatment and control.

    :param args: The pipeline arguments dictionary - dictionary of arguments (as keys) and their values.
    """
    try:
        spc = SitesProfileCalculator(Path(args['OUTPUT_FOLDER_PATH']),
                                     args['curr_out_dir'])

        args['sites_profiles_dir'] = spc.calculate_sites_profile(args['csv_ca'])

        args['curr_in_dir'] = validation.set_curr_in_dir(args['OUTPUT_FOLDER_PATH'],
                                                               STEP_OUT_FOLDER)
        logging.info(f'{STEP_NAME} step done.')
    except Exception as e:
        logging.error(f'Failed to run {STEP_NAME} step, please read the following for more information.',
                      exc_info=False)
        # sys.exit()
        raise e
