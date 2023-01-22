"""
This Python module is used to run the BREAKS_CLASSIFY step on treatment and control experiments.

First, filter icebergs in treatment with reads only from one side of the injected tag or if the
signal of one side of the tag is low.

Secondly, filter low signal sites by control icebergs.

The sites which filtered classified to the noise csv.

Finally, scores and classifies each of the remaining iceberg in treatment to CRISPR activity csv or
spontaneous breaks by iterate the chromosomes and considering the signals of icebergs in the same location but
from different experiment, that is, The signal (read count) in treatment clusters/icebergs vs the signals
(read count) in control clusters/icebergs with normalization of control to treatment.
Icebergs sites can be classified into one of the following three csv files -
CRISPR activity csv, spontaneous breaks csv and noise csv.
"""

import settings
from utils import generators, csv_to_df, logs, validation

import pandas as pd
import logging
from pathlib import Path

STEP_NAME = 'BREAKS_CLASSIFY'
STEP_OUT_FOLDER = 'BREAKS_CLASSIFY_RESULTS'
STEP_LOG_FILE = 'BREAKS_CLASSIFY.log'


class BreakSitesClassifier(object):

    def __init__(self,
                 crispr_activity_threshold: float,
                 noise_bins_and_control_mapq_percentile: list,
                 treatment_read_count: int,
                 control_read_count: int,
                 in_folder: Path,
                 out_folder: Path) -> None:
        """
        Initiate the BreakSitesClassifier object for the classification process.

        :param crispr_activity_threshold: The threshold for separate CRISPR activities sites from spontaneous breaks
         sites (0-1).
        :param noise_bins_and_control_mapq_percentile: lists of bins, where we defined bin as
         [[bin-min-icebergs-mapq, bin-max-icebergs-mapq], bin-icebergs-control-percentile, bin-icebergs-min-read-count].
        :param treatment_read_count:
        :param control_read_count:
        :param in_folder: The folder where the input csv file is stored.
        :param out_folder: The folder where the output csv files will be saved.
        """
        self.crispr_activity_threshold = crispr_activity_threshold
        self.noise_bins_and_control_mapq_percentile = noise_bins_and_control_mapq_percentile
        self.treatment_total_read_count = float(treatment_read_count)
        self.control_total_read_count = float(control_read_count)
        self.in_folder = in_folder
        self.out_folder = out_folder

    def get_chromosome_statistics(self,
                                  chromosome_df: pd.DataFrame,
                                  experiment_name: str) -> str:
        """
        Returns chromosome statistic information.
        The report is based on the chromosome dataframe (of treatment or control experiment).

        :param chromosome_df: The chromosome dataframe.
        :param experiment_name: The experiment name.

        :return: string with the statistic information about the chromosome in the experiment.
        """
        experiment_suffix = 't' if experiment_name == 'treatment' else 'm'
        longest_index = chromosome_df[f'iceberg length {experiment_suffix}'].idxmax()
        largest_index = chromosome_df[f'read count {experiment_suffix}'].idxmax()
        return settings.EXPERIMENT_CHROMOSOME_STAT_TEMPLATE.format(
           experiment_name=experiment_name,
           icebergs_count=chromosome_df.shape[0],
           reads_count=chromosome_df[f'read count {experiment_suffix}'].sum(),
           avg_length=chromosome_df[f'iceberg length {experiment_suffix}'].mean(),
           avg_reads_count=chromosome_df[f'read count {experiment_suffix}'].mean(),
           longest_length=chromosome_df.loc[longest_index, f'iceberg length {experiment_suffix}'],
           longest_reads_count=chromosome_df.loc[longest_index, f'read count {experiment_suffix}'],
           longest_start=chromosome_df.loc[longest_index, f'start index {experiment_suffix}'],
           longest_end=chromosome_df.loc[longest_index, f'end index {experiment_suffix}'],
           largest_length=chromosome_df.loc[largest_index, f'iceberg length {experiment_suffix}'],
           largest_reads_count=chromosome_df.loc[largest_index, f'read count {experiment_suffix}'],
           largest_start=chromosome_df.loc[largest_index, f'start index {experiment_suffix}'],
           largest_end=chromosome_df.loc[largest_index, f'end index {experiment_suffix}']
        )

    def filter_spontaneous_breaks_sites(self,
                                        scored_df: pd.DataFrame) -> (pd.DataFrame, pd.DataFrame):
        """
        Filters only the sites which are classified as spontaneous breaks sites.
        Those sites are the sites with "crispr activity score" < crispr_activity_threshold.

        :param scored_df: The merged icebergs dataframe (of treatment and control) after scoring all sites.

        :return: A merged icebergs dataframe which contain only the sites that classified as spontaneous break sites.
        """
        df = scored_df.copy()
        ca_df = df[df['crispr activity score'] >= self.crispr_activity_threshold]
        sb_df = df[df['crispr activity score'] < self.crispr_activity_threshold]

        ca_df = ca_df.sort_values(by=['read count t'], ascending=False)
        sb_df = sb_df.sort_values(by=['read count t'], ascending=False)

        logging.info(f'amount of CRISPR activities (on/off target) sites detected: {ca_df.shape[0]}')
        logging.info(f'amount of spontaneous breaks detected: {sb_df.shape[0]}')

        return sb_df, ca_df

    def filter_noise_sites(self,
                           scored_df: pd.DataFrame) -> (pd.DataFrame, pd.DataFrame):
        """
        Filters only the sites which are flagged as noise sites.
        This is done by looking at the "iceberg filtered t" column (1 if the treatment iceberg was filtered, else 0).

        :param scored_df: The joint icebergs dataframe (from treatment and control) after scoring
         each site in all chromosomes.

        :return: The joint icebergs dataframe which contain only the sites which are classified as noise sites.
        """
        noise_df = scored_df[scored_df['iceberg filtered'] == 1]
        rest_df = scored_df[scored_df['iceberg filtered'] == 0]

        logging.info(f'amount of noise icebergs detected: {noise_df.shape[0]}')

        noise_df = noise_df.sort_values(by=['read count t', 'read count m'],
                                        ascending=False)
        return noise_df, rest_df

    # def confidence_interval(self, data, percentile, confidence):
    #     # confidence interval for percentile
    #
    #     res = st.mstats.mquantiles_cimj(data=data, prob=percentile, alpha=1-confidence)
    #     return res[0][0], res[1][0]

    def get_quantile(self, icebergs_depth, quantile):
        # confidence interval for percentile
        return icebergs_depth.quantile(quantile) * (self.treatment_total_read_count / self.control_total_read_count)

    def flag_treatment_sites_with_pair_end_reads_from_one_tag_side_only(self,
                                                                        icebergs_df: pd.DataFrame) -> pd.DataFrame:
        """
        Flag treatment sites with less then two pair end reads in both sides of the tag.

        :param icebergs_df: The merged iceberg dataframe of treatment and control experiments.

        :return: The flagged merged iceberg dataframe.
        """
        logging.info("Flag treatment sites with less then two reads libraries (forward and reverse libraries)")
        flagged_icebergs_df = icebergs_df.copy()

        empty_iceberg_sites_treatment_indexes = flagged_icebergs_df[flagged_icebergs_df['chromosome t'] == '*'].index
        flagged_icebergs_df.loc[empty_iceberg_sites_treatment_indexes, 'iceberg filtered'] = 1

        for i in flagged_icebergs_df.index:
            orientations = flagged_icebergs_df.loc[i, "orientations t"]
            # if not all(orientation in orientations for orientation in ['tag_left_side', 'tag_right_side']):
            pair_end_reads_from_both_sides_of_the_tag = all(orientation in orientations
                                                            for orientation in ['tag_left_side', 'tag_right_side'])

            if pair_end_reads_from_both_sides_of_the_tag:
                low_signal_from_tag_sides = orientations['tag_left_side'] < 2 or orientations['tag_right_side'] < 2
                if low_signal_from_tag_sides:
                    flagged_icebergs_df.loc[i, 'iceberg filtered'] = 1
            else:
                flagged_icebergs_df.loc[i, 'iceberg filtered'] = 1

        return flagged_icebergs_df

    # def flag_treatment_sites_with_less_then_two_reads_libraries(self,
    #                                                             icebergs_df: pd.DataFrame) -> pd.DataFrame:
    #     """
    #
    #     """
    #     logging.info("Flag treatment sites with less then two reads libraries (forward and reverse libraries)")
    #     flagged_icebergs_df = icebergs_df.copy()
    #     flagged_icebergs_df['iceberg filtered'] = 0
    #     for i in flagged_icebergs_df.index:
    #         cut_position_candidates = flagged_icebergs_df.loc[i, "cut-position-candidates t"]
    #         site_libraries = set(primer for candidate in cut_position_candidates.keys() for primer in cut_position_candidates[candidate].keys())
    #         if not all(primer in site_libraries for primer in ['forward-primer', 'reverse-primer']):
    #             flagged_icebergs_df.loc[i, 'iceberg filtered'] = 1
    #
    #     return flagged_icebergs_df

    # def flag_treatment_sites_with_less_then_two_reads_libraries(self,
    #                                                             icebergs_df: pd.DataFrame) -> pd.DataFrame:
    #     """
    #
    #     """
    #     logging.info("Flag treatment sites with less then two reads libraries (forward and reverse libraries)")
    #     flagged_icebergs_df = icebergs_df.copy()
    #     flagged_icebergs_df['iceberg filtered'] = 0
    #     for i in flagged_icebergs_df.index:
    #         cut_position_candidates = flagged_icebergs_df.loc[i, "cut-position-candidates t"]
    #         if cut_position_candidates == {}:
    #             flagged_icebergs_df.loc[i, 'iceberg filtered'] = 1
    #             continue
    #
    #         highest_candidate = max(cut_position_candidates,
    #                                 key=lambda candidate: sum(cut_position_candidates[candidate][primer]
    #                                                           for primer in cut_position_candidates[candidate]))
    #         highest_candidate_data = cut_position_candidates[highest_candidate]
    #         if sum(highest_candidate_data[primer] for primer in highest_candidate_data.keys()) < 5:
    #             flagged_icebergs_df.loc[i, 'iceberg filtered'] = 1
    #
    #     return flagged_icebergs_df

    # def flag_treatment_noise_sites_by_control(self,
    #                                           icebergs_df: pd.DataFrame) -> pd.DataFrame:
    #     """
    #     flag_treatment_noise_sites_by_control:
    #     flags clusters/icebergs sites in treatment to noise by control clusters/icebergs sites.
    #     This process done by looking at the icebergs read count after dividing the icebergs sites
    #     to three bins by their mapq (rounded median mapq of their reads) - [0, 1),[1, 50),[50, 61),
    #     in each bin, treatment icebergs with less then noise_threshold percentile of the control icebergs read count
    #     (in the corresponding bin) will consider as noise site.
    #
    #     :param icebergs_df: The merged iceberg dataframe of treatment and control experiments.
    #
    #     :return: The treatment icebergs data frame with new column - iceberg filtered (1 if the iceberg was filtered, else 0)
    #     """
    #     logging.info("Flag treatment low signal sites by control (by bins)")
    #
    #     flagged_icebergs_df = icebergs_df.copy()
    #     # flagged_icebergs_df['iceberg filtered'] = 0
    #
    #     empty_t_sites = flagged_icebergs_df[flagged_icebergs_df['chromosome t'] == '*'].index
    #     flagged_icebergs_df.loc[empty_t_sites, 'iceberg filtered'] = 1
    #
    #     noise_sites_count = 0
    #     for bin_data in self.noise_bins_and_control_mapq_percentile:
    #         bin_mapq_limits, percentile, min_read_count = bin_data
    #         bin_icebergs_t = icebergs_df[icebergs_df['mapq t'] >= bin_mapq_limits[0]]
    #         bin_icebergs_t = bin_icebergs_t[bin_icebergs_t['mapq t'] < bin_mapq_limits[1]]
    #
    #         bin_icebergs_m = icebergs_df[icebergs_df['mapq m'] >= bin_mapq_limits[0]]
    #         bin_icebergs_m = bin_icebergs_m[bin_icebergs_m['mapq m'] < bin_mapq_limits[1]]
    #
    #         if bin_icebergs_m.shape[0] > 0:
    #             fallback = False
    #             # icebergs_read_counts_precentile_ci = self.confidence_interval(data=bin_icebergs_m['read count m'],
    #             #                                                               percentile=percentile,
    #             #                                                               confidence=0.95)
    #             icebergs_read_counts_precentile_ci = self.get_quantile(data=bin_icebergs_m['read count m'],
    #                                                                    quantile=percentile)
    #             icebergs_read_counts_precentile_ci = [2, min(icebergs_read_counts_precentile_ci, min_read_count)]
    #         else:
    #             fallback = True
    #             icebergs_read_counts_precentile_ci = [2,  min_read_count]
    #
    #         fall_back_str = f'     !!!!!!! default threshold: 4 - no icebergs sites in control !!!!!!!' + '\n' if fallback else ''
    #
    #         noise_icebergs_index = bin_icebergs_t[bin_icebergs_t['read count t'] <= icebergs_read_counts_precentile_ci[1]].index
    #         flagged_icebergs_df.loc[noise_icebergs_index, 'iceberg filtered'] = 1
    #         noise_sites_count += len(noise_icebergs_index)
    #
    #         logging.info(f'\n\n' +
    #                      f'bin {bin_mapq_limits[0]}-{bin_mapq_limits[1]} filtering low signal sites info:' + f'\n{fall_back_str}'
    #                      f'     bin treatment total icebergs: {bin_icebergs_t.shape[0]}.' + '\n'
    #                      f'     bin control total icebergs: {bin_icebergs_m.shape[0]}.' + '\n'+ '\n'
    #                      f'     bin control icebergs read-count percentile: '
    #                      f'     {percentile}% confidence=0.95% (percentile conf-interval: [{icebergs_read_counts_precentile_ci[0]}, {icebergs_read_counts_precentile_ci[1]}])' + '\n'
    #                      f'     bin treatment icebergs read-count min threshold: {icebergs_read_counts_precentile_ci[1]}'+ '\n'
    #                      f'     bin treatment filtered icebergs: {len(noise_icebergs_index)}.\n\n')
    #
    #     logging.info(f'treatment total icebergs: {sum(1 for chrom in icebergs_df["chromosome t"] if chrom != "*")}.')
    #     logging.info(f'treatment noise icebergs: {noise_sites_count}.')
    #     return flagged_icebergs_df

    def flag_treatment_noise_sites_by_control(self,
                                              icebergs_df: pd.DataFrame) -> pd.DataFrame:
        """
        Flags clusters/icebergs sites in treatment to noise by control clusters/icebergs sites.
        This process done by looking at the icebergs read count after dividing the icebergs sites
        to three bins by their mapq (rounded median mapq of their reads) - [0, 1),[1, 50),[50, 61),
        in each bin, treatment icebergs with less then noise_threshold percentile of the control icebergs read count
        (in the corresponding bin) will consider as noise site.

        :param icebergs_df: The merged iceberg dataframe of treatment and control experiments.

        :return: The treatment icebergs data frame with new column - iceberg filtered
         (1 if the iceberg was filtered, else 0)
        """
        logging.info("Flag treatment low signal sites by control (by bins)")

        flagged_icebergs_df = icebergs_df.copy()
        # flagged_icebergs_df['iceberg filtered'] = 0

        noise_sites_count = 0
        for bin_data in self.noise_bins_and_control_mapq_percentile:
            bin_mapq_limits, quantile, min_read_count = bin_data
            bin_icebergs_t = icebergs_df[icebergs_df['mapq t'] >= bin_mapq_limits[0]]
            bin_icebergs_t = bin_icebergs_t[bin_icebergs_t['mapq t'] < bin_mapq_limits[1]]

            bin_icebergs_m = icebergs_df[icebergs_df['mapq m'] >= bin_mapq_limits[0]]
            bin_icebergs_m = bin_icebergs_m[bin_icebergs_m['mapq m'] < bin_mapq_limits[1]]

            if bin_icebergs_m.shape[0] > 0:
                # icebergs_read_counts_percentile_ci = self.confidence_interval(data=bin_icebergs_m['read count m'],
                #                                                               percentile=percentile,
                #                                                               confidence=0.95)
                icebergs_depth_quantile_value = self.get_quantile(icebergs_depth=bin_icebergs_m['iceberg depth m'],
                                                                  quantile=quantile)
                icebergs_depth_quantile_value = max(icebergs_depth_quantile_value, min_read_count)
                fall_back_str = ''
            else:
                icebergs_depth_quantile_value = min_read_count
                fall_back_str = f'      default threshold: {min_read_count} - no icebergs sites in control.'

            noise_sites_index = bin_icebergs_t[bin_icebergs_t['iceberg depth t'] <= icebergs_depth_quantile_value].index
            flagged_icebergs_df.loc[noise_sites_index, 'iceberg filtered'] = 1
            noise_sites_count += len(noise_sites_index)

            logging.info(f'\n\n' +
                         f'bin {bin_mapq_limits[0]}-{bin_mapq_limits[1]} filtering low signal sites info:' + '\n' 
                         f'     {fall_back_str}'
                         f'     bin treatment total icebergs: {bin_icebergs_t.shape[0]}.' + '\n'
                         f'     bin control total icebergs: {bin_icebergs_m.shape[0]}.' + '\n' + '\n'
                         f'     bin control icebergs read-count quantile: '
                         f'     {quantile}% (quantile value: {icebergs_depth_quantile_value})' + '\n'
                         f'     bin treatment filtered icebergs: {len(noise_sites_index)}.\n\n')

        logging.info(f'treatment total icebergs: {sum(1 for chrom in icebergs_df["chromosome t"] if chrom != "*")}.')
        logging.info(f'treatment noise icebergs: {noise_sites_count}.')
        return flagged_icebergs_df

    def classify(self, merged_csv_file_name: str) -> (Path, Path, Path):
        """
        First, filter icebergs in treatment with reads only from one side of the injected tag or if the
        signal of one side of the tag is low.
        Secondly, filter low signal sites by control icebergs - the sites which filtered classified to the noise csv.
        Finally, scores and classifies each of the remaining iceberg in treatment to CRISPR activity csv or
        spontaneous breaks by iterate the chromosomes and considering the signals of icebergs in the same location but
        from different experiment, that is, The signal (read count) in treatment clusters/icebergs vs the signals
        (read count) in control clusters/icebergs with normalization of control to treatment
        Icebergs sites can be classified into one of the following three csv files -
        CRISPR activity csv, spontaneous breaks csv and noise csv.

        :param merged_csv_file_name: The control icebergs csv file name.

        :return: The three csv files path (CRISPR activity csv, spontaneous breaks csv and noise csv).
        """
        noise_df_name = merged_csv_file_name.replace('.icebergs.csv', '-noise.icebergs.csv')
        crispr_activity_df_name = merged_csv_file_name.replace('.icebergs.csv', '-crispr-activity.icebergs.csv')
        spontaneous_breaks_df_name = merged_csv_file_name.replace('.icebergs.csv', '-spontaneous-breaks.icebergs.csv')
        icebergs_df = csv_to_df.get_merged_experiments_df(self.in_folder / merged_csv_file_name)
        icebergs_df['crispr activity score'] = None
        icebergs_df['iceberg filtered'] = 0

        pbar = logs.get_progress_bar(step_name=STEP_NAME,
                                     total=icebergs_df.shape[0],
                                     out_of="icebergs-site")

        icebergs_df = self.flag_treatment_sites_with_pair_end_reads_from_one_tag_side_only(icebergs_df)
        noise_df, rest_df = self.filter_noise_sites(icebergs_df)
        pbar.update(noise_df.shape[0])

        rest_df = self.flag_treatment_noise_sites_by_control(rest_df)
        noise_df_by_control, rest_df = self.filter_noise_sites(rest_df)
        pbar.update(noise_df_by_control.shape[0])

        noise_df = noise_df.append(noise_df_by_control, ignore_index=True)

        scored_df = pd.DataFrame()
        chromosomes = generators.get_icebergs_by_chromosomes_merged(rest_df)
        # iterate over all chromosomes names.
        for chromosome in chromosomes:
            stat_string = settings.UPPER_FRAME
            stat_string += settings.CHROMOSOME_TITLE_TEMPLATE.format(chromosome_name=chromosome['chromosome_name'])

            if all(chromosome['df']['chromosome t'] == '*'):
                stat_string += self.get_chromosome_statistics(chromosome['df'], 'control')
                chromosome['df']['crispr activity score'] = 0

            elif all(chromosome['df']['chromosome m'] == '*'):
                stat_string += self.get_chromosome_statistics(chromosome['df'], 'treatment')
                chromosome['df']['crispr activity score'] = 1
            else:
                stat_string += self.get_chromosome_statistics(chromosome['df'], 'treatment')
                stat_string += self.get_chromosome_statistics(chromosome['df'], 'control')

                for i in chromosome['df'].index:
                    if chromosome['df'].loc[i, 'chromosome t'] == '*':
                        chromosome['df'].loc[i, 'crispr activity score'] = 0
                    elif chromosome['df'].loc[i, 'chromosome m'] == '*':
                        chromosome['df'].loc[i, 'crispr activity score'] = 1
                    else:
                        depth_t = chromosome['df'].loc[i, 'iceberg depth t']
                        depth_m = chromosome['df'].loc[i, 'iceberg depth m']
                        if depth_t == depth_m == 0:
                            chromosome['df'].loc[i, 'crispr activity score'] = 0
                        else:
                            depth_normalized_m = depth_m
                            depth_normalized_m *= (self.treatment_total_read_count / self.control_total_read_count)
                            numerator = depth_t - depth_normalized_m
                            denominator = ((depth_t ** 2 + depth_normalized_m ** 2) ** 0.5)
                            chromosome['df'].loc[i, 'crispr activity score'] = numerator / denominator

            stat_string += settings.LOWER_FRAME
            logging.info('\n' + stat_string)

            scored_df = scored_df.append(chromosome['df'])

            pbar.update(chromosome['df'].shape[0])

        spontaneous_breaks_df, crispr_activity_df = self.filter_spontaneous_breaks_sites(scored_df)

        scored_df.to_csv(path_or_buf=self.out_folder / merged_csv_file_name)
        noise_df.to_csv(path_or_buf=self.out_folder / noise_df_name)
        crispr_activity_df.to_csv(path_or_buf=self.out_folder / crispr_activity_df_name)
        spontaneous_breaks_df.to_csv(path_or_buf=self.out_folder / spontaneous_breaks_df_name)

        return Path(STEP_OUT_FOLDER, crispr_activity_df_name),\
            Path(STEP_OUT_FOLDER, spontaneous_breaks_df_name),\
            Path(STEP_OUT_FOLDER, noise_df_name)


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
    logs.set_logging_config(STEP_NAME, Path(args['OUTPUT_FOLDER_PATH'], 'DEBUG', 'LOGS', STEP_LOG_FILE))
    args['curr_out_dir'] = validation.set_curr_out_dir(args['OUTPUT_FOLDER_PATH'], STEP_OUT_FOLDER)
    validate_arguments(args)


def main(args: dict) -> None:
    """
    Run the BREAKS_CLASSIFY step on treatment and control.

    :param args: The pipeline arguments dictionary - dictionary of arguments (as keys) and their values.
    """
    try:
        previous_step_name = next(reversed(args['remaining_reads']))
        bc = BreakSitesClassifier(args['HYPERPARAMATERS']['CRISPR_ACTIVITY_THRESHOLD'],
                                  args['HYPERPARAMATERS']['NOISE_BINS_AND_CONTROL_MAPQ_PERCENTILE'],
                                  args['remaining_reads'][previous_step_name][0][1],  # treatment total read count.
                                  args['remaining_reads'][previous_step_name][1][1],  # control total read count.
                                  args['curr_in_dir'],
                                  args['curr_out_dir'])

        args['csv_ca'], args['csv_sb'], args['csv_noise'] = bc.classify(args['tx_and_control_merged_icebergs'])

        args['curr_in_dir'] = validation.set_curr_in_dir(args['OUTPUT_FOLDER_PATH'],
                                                         STEP_OUT_FOLDER)
        logging.info(f'{STEP_NAME} step done.')
    except Exception as e:
        logging.error(f'Failed to run {STEP_NAME} step, please read the following for more information.',
                      exc_info=False)
        # sys.exit()
        raise e
