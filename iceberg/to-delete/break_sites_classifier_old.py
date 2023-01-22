import scipy.stats as st
import numpy as np
import pandas as pd
import logging
from pathlib import Path
import settings
import utils.validation
from utils import directories_operations
from utils.logs import set_logging_config
classifier_dir = 'CLASSIFY-RESULTS'


class BreakSitesClassifier:

    def __init__(self, crispr_activity_threshold: float, noise_threshold:float, distance: int) -> None:
        """
        __init__:
        Initiate the breaks_classifier object with the thresholds for the classification process.

        :param crispr_activity_threshold: The threshold for separate CRISPR activities sites from spontaneous breaks sites (0-1).
        :param noise_threshold: The threshold for separate noise sites from the rest (0-1)
        :param distance: The distance (in bps) for unite nearby icebergs.
        """
        self.crispr_activity_threshold = crispr_activity_threshold
        self.noise_threshold = noise_threshold
        self.distance = distance

    # def calc_scores(self, t_iceberg: pd.Series, m_iceberg: pd.Series, t_sum_reads: int) -> (float, float, int):
    #     """
    #     calc_scores:
    #     Calculates the CRISPR activity score, noise score and max read count (old - noise score and max read count are not in used)
    #     of a site using the read count of each iceberg (the icebergs index should intersect, if there is only
    #     iceberg site in one of the experiment then the second experiment iceberg should be None).
    #
    #     The function handles the cases where there is iceberg from one experiment (treatment or mock) only in the site.
    #
    #     :param t_iceberg: The iceberg from treatment (None if there only mock iceberg in the site)
    #     :param m_iceberg: The iceberg from mock (None if there only treatment iceberg in the site)
    #     :param t_sum_reads: The read amount in treatment experiment.
    #
    #     :return: CRISPR activity score, noise score and max read count of the site
    #     (either read count of treatment or read count of mock).
    #     """
    #     t_read_count = 0 if t_iceberg is None or t_iceberg['iceberg filtered'] == 1 else t_iceberg["read count"]
    #     m_read_count = 0 if m_iceberg is None else m_iceberg["read count"]
    #     t_read_count_d = 1 if t_read_count == m_read_count == 0 else t_read_count
    #     crispr_activity_score = (t_read_count - m_read_count) / ((t_read_count_d ** 2 + m_read_count ** 2) ** 0.5)
    #     max_read_counts = max(t_read_count, m_read_count)
    #     noise_score = max_read_counts / t_sum_reads
    #     return crispr_activity_score, noise_score, max_read_counts

    # new calc_scores

    def calc_crispr_activity_score(self, t_iceberg: pd.Series, m_iceberg: pd.Series) -> (float, float, int):
        """
        calc_scores:
        Calculates the CRISPR activity score and max read count (old - noise score and max read count are not in used)
        of a site using the read count of each iceberg (the icebergs index should intersect, if there is only
        iceberg site in one of the experiment then the second experiment iceberg should be None).

        The function handles the cases where there is iceberg from one experiment (treatment or control) only in the site.

        :param t_iceberg: The iceberg from treatment (None if there only control iceberg in the site)
        :param m_iceberg: The iceberg from control (None if there only treatment iceberg in the site)

        :return: CRISPR activity score, noise score and max read count of the site
        (either read count of treatment or read count of control).
        """
        t_read_count = 0 if t_iceberg is None or t_iceberg['iceberg filtered'] == 1 else t_iceberg["read count"]
        m_read_count = 0 if m_iceberg is None else m_iceberg["read count"]
        t_read_count_d = 1 if t_read_count == m_read_count == 0 else t_read_count
        crispr_activity_score = (t_read_count - m_read_count) / ((t_read_count_d ** 2 + m_read_count ** 2) ** 0.5)
        return crispr_activity_score

    def score_site(self, t_iceberg: pd.Series, m_iceberg: pd.Series) -> pd.Series:
        """
        score_site:
        Calculates the crispr activity score of a site using the read count
        of each iceberg and prepare the new row for the scored dataframe (the icebergs index should intersect,
        if there is only one iceberg in the site then the second iceberg should be None).

        :param t_iceberg: Iceberg from treatment (None if there only control iceberg in the site)
        :param m_iceberg: Iceberg from control (None if there only treatment iceberg in the site)

        :return: Crispr activity score, noise score and max read count of the site.
        """
        crispr_activity_score = self.calc_crispr_activity_score(t_iceberg, m_iceberg)
        score_info = pd.Series([crispr_activity_score], index=['crispr activity score'])
        if m_iceberg is None:
            score_info = score_info.append(pd.Series(data=np.array(t_iceberg),
                                                     index=[k + ' t' for k in t_iceberg.index]))

            score_info = score_info.append(pd.Series(data=['*' if k not in ['read count', 'mapq','peak'] else 0 for k in t_iceberg.index],
                                                     index=[k + ' m' for k in t_iceberg.index]))
        elif t_iceberg is None:
            score_info = score_info.append(pd.Series(['*' if k not in ['read count', 'mapq', 'peak'] else 0 for k in m_iceberg.index],
                                                 index=[k + ' t' for k in m_iceberg.index]))
            score_info['iceberg filtered t'] = 1
            score_info = score_info.append(pd.Series(np.array(m_iceberg), index=[k + ' m' for k in m_iceberg.index]))
        else:
            score_info = score_info.append(pd.Series(np.array(t_iceberg), index=[i + ' t' for i in t_iceberg.index]))
            score_info = score_info.append(pd.Series(np.array(m_iceberg), index=[i + ' m' for i in m_iceberg.index]))
        return score_info

    def unite_icebergs(self, iceberg: pd.Series, iceberg_to_add: pd.Series) -> pd.Series:
        """
        unite_icebergs:
        Unites the given pair of iceberg into one larger iceberg

        :param iceberg: The first iceberg
        :param iceberg_to_add: The next iceberg (from the same experiment) that need to be unite with the iceberg (previous param).

        :return: The united iceberg.
        """
        gap_size = iceberg_to_add['start index'] - iceberg['end index'] - 1
        iceberg['end index'] = iceberg_to_add['end index']
        iceberg['iceberg length'] += gap_size + iceberg_to_add['iceberg length']
        iceberg['mapq'] = (iceberg['mapq'] * iceberg['read count'] + iceberg_to_add['mapq'] * iceberg_to_add['read count'])
        iceberg['mapq'] /= (iceberg['read count'] + iceberg_to_add['read count'])
        iceberg['mapq'] = round(iceberg['mapq'])
        iceberg['read count'] += iceberg_to_add['read count']
        if 'iceberg filtered' in iceberg.index:
            iceberg['iceberg filtered'] = iceberg['iceberg filtered'] and iceberg_to_add['iceberg filtered']
        return iceberg

    def unite_icebergs_if_needed(self, t_icebergs_df: pd.DataFrame, m_icebergs_df: pd.DataFrame, t_iceberg: pd.Series, m_iceberg: pd.Series, i: int, j: int) -> (int, int):
        """
        unite_icebergs_if_needed:
        given each of the experiment icebergs dataframe of specific chromosome (treatment and control), and the indices of
        overlapping icebergs (one from treatment and one from control), the function unite the iceberg from each experiment
        with the following iceberg in the same experiment only if the distance of following iceberg from the
        current iceberg is smaller then icebergs_distance.
        The process is done until the condition is no longer true.
        The function return the new indices.

        :param t_icebergs_df: A dataframe of sorted icebergs from specific chromosome of treatment experiment.
        :param m_icebergs_df: A dataframe of sorted icebergs from specific chromosome of control experiment.
        :param t_iceberg: The current iceberg from treatment.
        :param m_iceberg: The current iceberg from control.
        :param i: The t_iceberg index in t_icebergs_df.
        :param j: The m_iceberg index in m_icebergs_df.

        :return: The modified i and j if unification was done, else same i,j.
        """
        unification_check_needed = True
        while unification_check_needed and (i < t_icebergs_df.shape[0] - 1 or j < m_icebergs_df.shape[0] - 1):
            unification_check_needed_t, unification_check_needed_m = False, False
            while i < t_icebergs_df.shape[0] - 1:
                t_next_iceberg = t_icebergs_df.iloc[i + 1, :]
                if t_next_iceberg['start index'] <= max(m_iceberg['end index'], t_iceberg['end index']) + self.distance:
                    unification_check_needed_m = True
                    t_iceberg = self.unite_icebergs(t_iceberg, t_next_iceberg)
                    i += 1
                else:
                    unification_check_needed_t = False
                    break

            while j < m_icebergs_df.shape[0] - 1:
                m_next_iceberg = m_icebergs_df.iloc[j + 1, :]
                if m_next_iceberg['start index'] <= max(t_iceberg['end index'], m_iceberg['end index']) + self.distance:
                    unification_check_needed_t = True
                    m_iceberg = self.unite_icebergs(m_iceberg, m_next_iceberg)
                    j += 1
                else:
                    unification_check_needed_m = False
                    break

            unification_check_needed = unification_check_needed_t or unification_check_needed_m
        return i, j

    def score_chromosome_sites(self, t_icebergs_df: pd.DataFrame = None, m_icebergs_df: pd.DataFrame = None) -> pd.DataFrame:
        """
        score_chromosome_sites:

        Scores each of the remaining cluster/iceberg sites in treatment for the later classification to CRISPR activity csv
        or spontaneous breaks csv, this is done by - first we check in each experiment if there are some others nearby
        icebergs with bps-distance < icebergs_distance from the overlapping icebergs, if so, we unite them to one iceberg.
        Then, we consider the signals of clusters/icebergs in the same location but from different experiment, that is,
        given two overlap icebergs - one from treatment and one from control (after we unite each as explained above), we
        consider The signal (read count) in treatment clusters/icebergs vs the the signals (read count) in control clusters/icebergs
        and normalize those signals to get a score value.

        The icebergs that classified as noise in filter_treatment_by_control will treat as they have read-count = 0 and mapq = 0,
        unless they been unite with at least one iceberg that is not classified as noise.

        The function handles the cases where there are icebergs sites within the chromosome only in one of the experiments
        (treatment or control).

        :param t_sum_reads: The read amount in treatment experiment.
        :param t_icebergs_df: A dataframe of sorted icebergs from specific chromosome of treatment experiment.
        :param m_icebergs_df: A dataframe of sorted icebergs from specific chromosome of control experiment.

        :return: A dataframe which contains both experiments icebergs joint by chromosome index and scored for the
                 classification process.
        """
        i = 0
        j = 0
        scored_df = []
        if t_icebergs_df is not None and m_icebergs_df is not None:
            while i < t_icebergs_df.shape[0] and j < m_icebergs_df.shape[0]:
                t_iceberg = t_icebergs_df.iloc[i, :].copy()
                m_iceberg = m_icebergs_df.iloc[j, :].copy()
                if m_iceberg['start index'] > t_iceberg['end index'] + self.distance:
                    scored_df.append(self.score_site(t_iceberg, None))
                    i += 1
                elif t_iceberg['start index'] > m_iceberg['end index'] + self.distance:
                    scored_df.append(self.score_site(None, m_iceberg))
                    j += 1
                else:
                    i, j = self.unite_icebergs_if_needed(t_icebergs_df, m_icebergs_df, t_iceberg, m_iceberg, i, j)
                    scored_df.append(self.score_site(t_iceberg.copy(), m_iceberg.copy()))
                    i += 1
                    j += 1

        if t_icebergs_df is not None:
            while i < t_icebergs_df.shape[0]:
                t_iceberg = t_icebergs_df.iloc[i, :]
                scored_df.append(self.score_site(t_iceberg, None))
                i += 1
        if m_icebergs_df is not None:
            while j < m_icebergs_df.shape[0]:
                m_iceberg = m_icebergs_df.iloc[j, :]
                scored_df.append(self.score_site(None, m_iceberg))
                j += 1
        scored_df = pd.DataFrame(scored_df)
        return scored_df

    def set_df_types(self, icebergs_df_path):
        """
        set_df_types:
        Reads the experiment (treatment or control) icebergs (reads clustered by their chromosome index)
        from csv file and convert them to dataframe properly.

        :param icebergs_df_path: The absolute path to the icebergs csv file.

        :return: An experiment (control or treatment) icebergs dataframe.
        """
        try:
            icebergs_df = pd.read_csv(icebergs_df_path)
            icebergs_df['chromosome'] = icebergs_df['chromosome'].apply(lambda x: str(x))
            icebergs_df['mapq'] = icebergs_df['mapq'].astype(float)
            icebergs_df['start index'] = icebergs_df['start index'].astype(float)
            icebergs_df['end index'] = icebergs_df['end index'].astype(float)
            icebergs_df['read count'] = icebergs_df['read count'].astype(float)
            return icebergs_df
        except Exception as e:
            logging.error(f'cant read csv file properly :{icebergs_df_path}', exc_info=False)
            raise e

    def get_chromosome_statistics(self, chromosome_df: pd.DataFrame, experiment_name: str) -> str:
        """
        add_chromosome_statistics:
        returns chromosome statistic information.
        The report is based on the chromosome dataframe (of treatment or control experiment).

        :param chromosome_df: The chromosome dataframe.
        :param experiment_name: The experiment name.

        :return: string with the statistic information about the chromosome in the experiment.
        """
        longest_index = chromosome_df['iceberg length'].idxmax()
        largest_index = chromosome_df['read count'].idxmax()
        return settings.EXPERIMENT_CHROMOSOME_STAT_TEMPLATE.format(experiment_name=experiment_name,
                                                                   icebergs_count=chromosome_df.shape[0],
                                                                   reads_count=chromosome_df['read count'].sum(),
                                                                   avg_length=chromosome_df['iceberg length'].mean(),
                                                                   avg_reads_count=chromosome_df['read count'].mean(),
                                                                   longest_length=chromosome_df.loc[longest_index, 'iceberg length'],
                                                                   longest_reads_count=chromosome_df.loc[longest_index, 'read count'],
                                                                   longest_start=chromosome_df.loc[longest_index, 'start index'],
                                                                   longest_end=chromosome_df.loc[longest_index, 'end index'],
                                                                   largest_length=chromosome_df.loc[largest_index, 'iceberg length'],
                                                                   largest_reads_count=chromosome_df.loc[largest_index, 'read count'],
                                                                   largest_start=chromosome_df.loc[largest_index, 'start index'],
                                                                   largest_end=chromosome_df.loc[largest_index, 'end index'])

    def get_crispr_activity_icebergs(self, scored_df: pd.DataFrame) -> pd.DataFrame:
        """
        get_crispr_activity_icebergs:
        Filters only the sites which are classified as crispr activity sites.
        Those sites are the sites with "crispr activity score" >= crispr_activity_threshold.

        :param scored_df: The joint icebergs dataframe (from treatment and control) after scoring each site in all chromosomes.

        :return: The joint icebergs dataframe which contain only the sites which are classified as crispr activity sites.
        """
        df = scored_df.copy()
        df = df[df['crispr activity score'] >= self.crispr_activity_threshold]
        logging.info(f'amount of CRISPR activities (on/off target) sites detected: {df.shape[0]}')
        df = df.sort_values(by=['read count t'], ascending=False)
        return df

    def get_spontaneous_breaks_icebergs(self, scored_df: pd.DataFrame) -> pd.DataFrame:
        """
        get_spontaneous_breaks_icebergs:
        Filters only the sites which are classified as spontaneous breaks sites.
        Those sites are the sites with "crispr activity score" < crispr_activity_threshold.

        :param scored_df: The joint icebergs dataframe (from treatment and control) after scoring each site in all chromosomes.

        :return: The joint icebergs dataframe which contain only the sites which are classified as spontaneous breaks sites.
        """
        df = scored_df.copy()
        df = df[df['crispr activity score'] < self.crispr_activity_threshold]
        logging.info(f'amount of spontaneous breaks detected: {df.shape[0]}')
        # df['peak t'] = [peaks[-1] for peaks in list(df['peak t'])]
        # df['peak m'] = [peaks[-1] for peaks in list(df['peak m'])]

        df = df.sort_values(by=['read count t'], ascending=False)
        return df

    def get_noise_icebergs(self, scored_df: pd.DataFrame) -> (pd.DataFrame, pd.DataFrame):
        """
        get_noise_icebergs:
        Filters only the sites which are classified as noise sites.
        This is done by looking at the "iceberg filtered t" column (1 if the treatment iceberg was filtered, else 0).

        :param scored_df: The joint icebergs dataframe (from treatment and control) after scoring each site in all chromosomes.

        :return: The joint icebergs dataframe which contain only the sites which are classified as noise sites.
        """
        df = scored_df.copy()
        noise_df = scored_df[scored_df['iceberg filtered t'].apply(lambda x: float(x) if x != '*' else 1)==1]
        rest_df = scored_df[scored_df['iceberg filtered t'].apply(lambda x:  float(x) if x != '*' else 1)==0]
        logging.info(f'amount of noise icebergs detected: {noise_df.shape[0]}')
        # noise_df['peak t'] = [peaks if peaks == 0 else peaks[-1] for peaks in list(noise_df['peak t'])]
        # noise_df['peak m'] = [peaks if peaks == 0 else peaks[-1] for peaks in list(noise_df['peak m'])]

        noise_df = noise_df.sort_values(by=['read count t', 'read count m'],
                                        ascending=False)
        return noise_df, rest_df

    def split_icebergs_by_chromosomes(self, t_icebergs_df: pd.DataFrame, m_icebergs_df: pd.DataFrame) -> dict:
        """
        split_icebergs_by_chromosomes:
        Creates dictionary where the key is the chromosome name and the value
        is also a dictionary which contains at most two keys,'t df' and 'm df', where the corresponding values are
        the icebergs from treatment and control respectively at the specific chromosome.

        :param t_icebergs_df: The icebergs dataframe of the treatment experiment.
        :param m_icebergs_df: The icebergs dataframe of the control experiment.

        :return: The created dictionary.

        """
        chromosomes = {}
        names = t_icebergs_df['chromosome'].unique().tolist()
        for chromosome_name in names:
            chromosomes[chromosome_name] = {'t df': t_icebergs_df.loc[t_icebergs_df.chromosome == chromosome_name]}

        names = m_icebergs_df['chromosome'].unique().tolist()
        for chromosome_name in names:
            if chromosome_name not in chromosomes:
                chromosomes[chromosome_name] = {}
            chromosomes[chromosome_name]['m df'] = m_icebergs_df.loc[m_icebergs_df.chromosome == chromosome_name]
        return chromosomes

    def classify_sites_by_scores(self, df: pd.DataFrame, out_dir: Path, res_file_name: str) -> (str, str, str):
        """
        classify_sites_by_scores:
        After the outer-join of treatment and control dataframe and after the sites scoring, split the clusters/icebergs
        sites in to three csv files - crispr activity csv, spontaneous breaks csv and noise csv.

        :param df: The outer-joint scored icebergs dataframe .
        :param out_dir: The directory where the three classified csv files will be saved.
        :param res_file_name: The wanted base name for the scored icebergs dataframe csv file.

        :return: Three csv files name.
        """
        crispr_activity_df_name = res_file_name.replace('.icebergs.csv',
                                                        '-crispr-activity.icebergs.csv')

        spontaneous_breaks_df_name = res_file_name.replace('.icebergs.csv',
                                                           '-spontaneous-breaks.icebergs.csv')

        noise_df_name = res_file_name.replace('.icebergs.csv',
                                              '-noise.icebergs.csv')

        df = df[[x for x in df.columns if x not in ['CIGAR m', 'CIGAR t']]]
        df.to_csv(path_or_buf=out_dir / res_file_name)

        noise_df, rest_df = self.get_noise_icebergs(df)
        noise_df.to_csv(path_or_buf=out_dir / noise_df_name)

        crispr_activity_df = self.get_crispr_activity_icebergs(rest_df)
        crispr_activity_df.to_csv(path_or_buf=out_dir / crispr_activity_df_name)

        spontaneous_breaks_df = self.get_spontaneous_breaks_icebergs(rest_df)
        spontaneous_breaks_df.to_csv(path_or_buf=out_dir / spontaneous_breaks_df_name)

        return crispr_activity_df_name, spontaneous_breaks_df_name, noise_df_name

    def confidence_interval(self, data, percentile, confidence):
        # confidence interval for percentile
        res = st.mstats.mquantiles_cimj(data=data, prob=percentile, alpha=1-confidence)
        return res[0][0], res[1][0]

    def filter_treatment_by_control(self, t_icebergs_df: pd.DataFrame, m_icebergs_df: pd.DataFrame) -> pd.DataFrame:
        """
        filter_treatment_by_control:
        Classifies clusters/icebergs sites in treatment to noise by control clusters/icebergs sites.
        This process done by looking at the icebergs read count after dividing the icebergs sites
        to three bins by their mapq (rounded median mapq of their reads) - [0, 1),[1, 50),[50, 61),
        in each bin, treatment icebergs with less then noise_threshold percentile of the control icebergs read count
        (in the corresponding bin) will consider as noise site.

        :param t_icebergs_df: The icebergs dataframe of the treatment experiment.
        :param m_icebergs_df: The icebergs dataframe of the control experiment.

        :return: The treatment icebergs data frame with new column - iceberg filtered (1 if the iceberg was filtered, else 0)
        """
        t_filtered_icebergs_df = t_icebergs_df.copy()
        t_filtered_icebergs_df['iceberg filtered'] = 0
        count = 0
        for bin in [(0, 1), (1, 50), (50, 61)]: #[(0, 10), (10, 20), (20, 30), (30, 40), (40, 50), (50, 61)]:#[(0, 50), (50, 61)]:
            bin_icebergs_t = t_icebergs_df[t_icebergs_df['mapq'] >= bin[0]]
            bin_icebergs_t = bin_icebergs_t[bin_icebergs_t['mapq'] < bin[1]]

            bin_icebergs_m = m_icebergs_df[m_icebergs_df['mapq'] >= bin[0]]
            bin_icebergs_m = bin_icebergs_m[bin_icebergs_m['mapq'] < bin[1]]

            # bin_icebergs_mean_read_counts_m = bin_icebergs_m['read count'].median()
            # bin_icebergs_std_read_counts_m = bin_icebergs_m['read count'].std()
            # icebergs_read_counts_precentile = bin_icebergs_mean_read_counts_m + 60 * bin_icebergs_std_read_counts_m

            #icebergs_read_counts_precentile = bin_icebergs_m['read count'].quantile(q=self.noise_threshold)
            # if bin == (50, 61):
            icebergs_read_counts_precentile_ci = self.confidence_interval(data=bin_icebergs_m['read count'],
                                                                          percentile=0.99,
                                                                          confidence=0.95)
                # icebergs_read_counts_precentile = bin_icebergs_m['read count'].quantile(q=0.90)
            # else:
            #     icebergs_read_counts_precentile_ci = self.confidence_interval(data=bin_icebergs_m['read count'],
            #                                                                percentile=0.95,
            #                                                                confidence=0.95)
                # icebergs_read_counts_precentile = bin_icebergs_m['read count'].quantile(q=self.noise_threshold)
            filtered_icebergs_index = bin_icebergs_t[bin_icebergs_t['read count'] <= icebergs_read_counts_precentile_ci[1]].index
            t_filtered_icebergs_df.loc[filtered_icebergs_index, 'iceberg filtered'] = 1
            count += len(filtered_icebergs_index)

            logging.info(f'bin {bin[0]}-{bin[1]} info:' + '\n'
                         f'treatment total icebergs: {bin_icebergs_t.shape[0]}.' + '\n'
                         f'treatment filtered icebergs: {len(filtered_icebergs_index)}.' + '\n'
                         f'treatment-read-count-threshold: control-{self.noise_threshold}-read-count-percentile-ci = ' + '\n'
                         f'{icebergs_read_counts_precentile_ci[0],icebergs_read_counts_precentile_ci[1]}.')

        logging.info(f'treatment total icebergs: {t_icebergs_df.shape[0]}.')
        logging.info(f'treatment filtered icebergs: {count}.')
        return t_filtered_icebergs_df

    def classify(self, in_dir: Path, out_dir: Path, df1_file_name: str, df2_file_name: str) -> (Path, Path, Path):
        """
        classify:
        First classifies clusters/icebergs in treatment to noise by control clusters/icebergs
        (see filter_treatment_by_control function description), then scores and classifies each of the remaining
        cluster/iceberg in treatment to CRISPR activity csv or spontaneous breaks by iterate the chromosomes
        and considering the signals of clusters/icebergs in the same location but from different experiment, that is,
        The signal (read count) in treatment clusters/icebergs vs the the signals (read count) in control clusters/icebergs.
        Clusters/Icebergs sites can be classified into one of the following three csv files -
        CRISPR activity csv, spontaneous breaks csv and noise csv.

        :param in_dir: The directory where the icebergs csv files from treatment and control are stored.
        :param out_dir: The directory where the three classified csv files will be saved.
        :param df1_file_name: The control icebergs csv file name (stored within in_dir, the output of read_uniter.py).
        :param df2_file_name: The treatment icebergs csv file name (stored within in_dir, the output of read_uniter.py).

        :return: The three csv files path (CRISPR activity csv, spontaneous breaks csv and noise csv).
        """
        m_icebergs_df = self.set_df_types(in_dir / df2_file_name)
        t_icebergs_df = self.set_df_types(in_dir / df1_file_name)

        t_icebergs_df = self.filter_treatment_by_control(t_icebergs_df,
                                                         m_icebergs_df)

        scored_df = pd.DataFrame(columns=['crispr activity score'] +
                                         [col + ' t' for col in t_icebergs_df.columns] +
                                         [col + ' m' for col in m_icebergs_df.columns])

        t_sum_read_count = t_icebergs_df['read count'].sum()
        chromosomes = self.split_icebergs_by_chromosomes(t_icebergs_df,
                                                         m_icebergs_df)
        # iterate over all chromosomes names.
        for chrom in chromosomes:
            stat_string = settings.UPPER_FRAME
            stat_string += settings.CHROMOSOME_TITLE_TEMPLATE.format(chromosome_name=chrom)
            if 'm df' not in chromosomes[chrom]:
                scored_df = pd.concat([scored_df, self.score_chromosome_sites(t_icebergs_df=chromosomes[chrom]['t df'])])
                stat_string += self.get_chromosome_statistics(chromosomes[chrom]['t df'], 'treatment')
            elif 't df' not in chromosomes[chrom]:
                scored_df = pd.concat([scored_df, self.score_chromosome_sites(m_icebergs_df=chromosomes[chrom]['m df'])])
                stat_string += self.get_chromosome_statistics(chromosomes[chrom]['m df'], 'control')
            else:
                stat_string += self.get_chromosome_statistics(chromosomes[chrom]['t df'], 'treatment')
                stat_string += self.get_chromosome_statistics(chromosomes[chrom]['m df'], 'control')
                print('bbb')
                scored_df = pd.concat([scored_df, self.score_chromosome_sites(t_icebergs_df=chromosomes[chrom]['t df'],
                                                                              m_icebergs_df=chromosomes[chrom]['m df'])])
            stat_string += settings.LOWER_FRAME
            logging.info('\n' + stat_string)
        res_file_name = df1_file_name.replace('.icebergs.csv', '') + '-' + \
                        df2_file_name.replace('.icebergs.csv', '') + '.icebergs.csv'

        crispr_activity_df_name, spontaneous_breaks_df_name, noise_df_name = self.classify_sites_by_scores(scored_df,
                                                                                                           out_dir,
                                                                                                           res_file_name)
        return Path(classifier_dir, crispr_activity_df_name),\
               Path(classifier_dir, spontaneous_breaks_df_name),\
               Path(classifier_dir, noise_df_name)


def validate_arguments(args: dict) -> None:
    """
    validate_arguments:
    Validates the step arguments from the given arguments dict.

    :param args: The arguments dict (see main function description).
    """
    try:
        utils.validation.validate_path_existence(args['curr_in_dir'])
        if args['csv1'] != '':
            utils.validation.validate_path_existence(Path(args['curr_in_dir'], args['csv1']))
        if args['csv2'] != '':
                utils.validation.validate_path_existence(Path(args['curr_in_dir'], args['csv2']))
        if args['csv1'] == '' or args['csv1'] == '' and args['csv2'] == '':
            raise ValueError()
    except ValueError as e:
        logging.error('not all necessary arguments for classify step where giving.')
        raise e


def prepare_step(args: dict) -> None:
    """
    prepare_step:
    Prepares the step by - setting the logging configurations, setting the output directory and
    validate the step arguments.

    :param args: A dict of arguments (as keys) and their values (see main function description)
    """
    set_logging_config('classify', Path(args['OUTPUT_FOLDER_PATH'], 'DEBUG', 'LOGS', 'classify.log'))
    args['curr_out_dir'] = utils.validation.set_curr_out_dir(args['OUTPUT_FOLDER_PATH'], classifier_dir)
    validate_arguments(args)


def main(args: dict) -> None:
    """
    main:
    Classifies treatment clusters/icebergs sites based on the control clusters/icebergs sites.

    :param args: A dict of arguments (as keys) and their values.
     Can be one of the following:
     1.The pipeline arguments dict - if you run analyzer_old.py.
     2.The step arguments dict - if you run this class directly.
    """
    try:
        bc = BreakSitesClassifier(args['HYPERPARAMATERS']['CRISPR_ACTIVITY_THRESHOLD'],
                                  0,
                                  args['HYPERPARAMATERS']['MAX_ICEBERG_DISTANCE'])

        csv_ca, csv_sb, csv_noise = bc.classify(args['curr_in_dir'],
                                                args['curr_out_dir'],
                                                args['csv1'],
                                                args['csv2'])

        args['csv_ca'] = csv_ca
        args['csv_sb'] = csv_sb
        args['csv_noise'] = csv_noise
        args['curr_in_dir'] = utils.validation.set_curr_in_dir(args['OUTPUT_FOLDER_PATH'],
                                                               classifier_dir)
        logging.info('classify step done.')
    except Exception as e:
        logging.error('error occurred while running classify step, please read the following for more information.',
                      exc_info=False)
        # sys.exit()
        raise e
