"""
This Python module is used to run the REPORTS step on treatment and control experiments.

The step generates the iceberg main report and report for each sites class/type, in addition
it gather the sites-guideRNA-alignments and the sites-profile html files in the REPORTS folder
for the reports links.

Each report contains up to 1000 sites s.t the sites with the max read count are chosen.
"""

import settings
from utils import csv_to_df, logs, validation

import os
import logging
import pandas as pd
import subprocess
from pathlib import Path
import sys
import plotly.express as px
import shutil
import re
import codecs

STEP_NAME = 'REPORTS'
STEP_OUT_FOLDER = 'REPORTS_RESULTS'
STEP_LOG_FILE = 'REPORTS.log'


class ReportsManager(object):

    def __init__(self,):
        """
        Initiate the ReportsManager object.
        """

    def get_bed_line(self,
                     icebergs_site: pd.Series,
                     html_alignments_dir_name: str = None) -> str:
        """
        Giving a site which contains icebergs (one iceberg from both treatment and control or one only from one of them)
        generates a line for the bed file for the create_report command. Also, adds the
        guide alignment fields to the line if html_alignments_dir_name is not None.

        :param icebergs_site: The site sequence.
        :param html_alignments_dir_name: The name of the directory that contains the alignments html files,
         the directory should be located in the GUIDERNA_ALIGNMENT_RESULTS directory (default is None).        """
        locus = "{0}\t{1}\t{2}\tviewport={0}:{1}-{2};" \
                "Reads_Count_Treatment={3};" \
                "Reads_Count_Control={4};" \
                "Depth_Treatment={5};" \
                "Depth_Control={6};" \
                "Iceberg_Length_Treatment={7};" \
                "Iceberg_Length_Control={8};" \
                "Mapq_Treatment={9};" \
                "Mapq_Control={10}"
        if html_alignments_dir_name is not None:
            locus += ";strand={11};" \
                     "cut_position={12};" \
                     "guideRNA_alignment={13}<br>{14};" \
                     "Hamming_distance={15};" \
                     "Levenshtein_distances={16}"
            locus = locus.format(icebergs_site.loc['chromosome'],
                                 icebergs_site.loc['start index'],
                                 icebergs_site.loc['end index'],
                                 icebergs_site.loc['read count t'],
                                 icebergs_site.loc['read count m'],
                                 icebergs_site.loc['iceberg depth t'],
                                 icebergs_site.loc['iceberg depth m'],
                                 icebergs_site.loc['iceberg length t'],
                                 icebergs_site.loc['iceberg length m'],
                                 icebergs_site.loc['mapq t'],
                                 icebergs_site.loc['mapq m'],
                                 icebergs_site.loc['strand'],
                                 icebergs_site.loc['cut position'],
                                 icebergs_site.loc['gRNA alignment'],
                                 icebergs_site.loc['short site alignment'],
                                 icebergs_site.loc['Hamming distance'],
                                 icebergs_site.loc['Levenshtein distances'])
        else:
            locus = locus.format(icebergs_site.loc['chromosome'],
                                 icebergs_site.loc['start index'],
                                 icebergs_site.loc['end index'],
                                 icebergs_site.loc['read count t'],
                                 icebergs_site.loc['read count m'],
                                 icebergs_site.loc['iceberg depth t'],
                                 icebergs_site.loc['iceberg depth m'],
                                 icebergs_site.loc['iceberg length t'],
                                 icebergs_site.loc['iceberg length m'],
                                 icebergs_site.loc['mapq t'],
                                 icebergs_site.loc['mapq m'])
        return locus

    def generate_bed_file(self,
                          icebergs_df: pd.DataFrame,
                          out_dir: Path,
                          bed_out_file_name: str,
                          html_alignments_dir_name: str = None) -> Path:
        """
        Generates the bed file for the html report from the given icebergs DataFrame
        (crispr activity, spontaneous break or noise).

        :param icebergs_df: The given icebergs DataFrame (crispr activity, spontaneous break or noise).
        :param out_dir: The directory where the bed file will be saved.
        :param bed_out_file_name: The wanted name for the html report file.
        :param html_alignments_dir_name: The name of the directory that contains the alignments html files,
         the directory should be located in the GUIDERNA_ALIGNMENT_RESULTS directory (default is None).
        """
        bed_file = open(out_dir / bed_out_file_name, 'w')
        with bed_file as f:
            for i in range(icebergs_df.shape[0]-1):
                iceberg = icebergs_df.iloc[i, :]
                f.write(self.get_bed_line(iceberg, html_alignments_dir_name) + '\n')
            iceberg = icebergs_df.iloc[-1, :]
            f.write(self.get_bed_line(iceberg, html_alignments_dir_name))
        bed_file.close()
        return out_dir / bed_out_file_name

    def add_features_to_report(self,
                               html_report_path: Path,
                               html_alignments_dir_name: Path = None) -> None:
        """
        Adds some features (titles and components) to the html report created by create_report (command by igv).
        if images_directory is not None embeds the alignment html to the report.

        :param html_report_path: The path to the html report created by create_report.
        :param html_alignments_dir_name: The name of the directory that contains the alignments html files,
         the directory should be located in the GUIDERNA_ALIGNMENT_RESULTS directory (default is None).
        """
        try:
            pattern_fix_header = re.compile(r"<h1>.*</h1><!--title-->")
            pattern_add_image = re.compile('<div id="igvContainer">')
            pattern_change_image_on_click = re.compile('const t2 = t.1..split."-".;')

            filename = html_report_path
            new_filename = str(html_report_path).replace(".html", "-new.html")

            with codecs.open(str(filename), 'r') as readfile, codecs.open(str(new_filename), 'w') as writefile:
                lines = readfile.readlines()
                for line in lines:
                    if pattern_fix_header.search(line):
                        writefile.write(settings.HTML_LOGO +
                                        settings.HTML_CENTER_HEADER_TEMPLATE.format(line=line))

                    elif html_alignments_dir_name is not None and pattern_add_image.search(line):
                        writefile.write(line)
                        writefile.write(settings.HTML_ALIGNMENTS_IFRAME)
                    elif html_alignments_dir_name is not None and pattern_change_image_on_click.search(line):
                        writefile.write(settings.HTML_CHANGE_ALIGNMENTS_IFRAME_ON_CLICK_TEMPLATE.
                                        format(os_path_sep=os.sep, html_alignments_dir_name=html_alignments_dir_name))
                        writefile.write(line)
                    else:
                        writefile.write(line)
        except Exception:
            logging.warning('Failed to add features to report.', exc_info=False)
            sys.exit()

    def generate_html_report(self,
                             bam1: Path,
                             icebergs_df: pd.DataFrame,
                             out_dir: Path,
                             html_out_name: str,
                             genome: Path,
                             bam2: Path = None,
                             html_alignments_dir_name: str = None) -> None:
        """
        Generates an html report from the given icebergs data frame.
        If html_alignments_dir_name is given adds features as cut position, GuideRNA alignment visualization,
        sites profiles tables and histograms to the report.

        The function should be call for each site type separately (CRISPR activity, Spontaneous Break and Noise sites.)

        :param bam1: The path to the sorted and indexed bam file of the treatment experiment
        :param icebergs_df: The given icebergs dataframe (crispr activity, spontaneous break or noise).
        :param out_dir: Directory where the html report will be saved.
        :param html_out_name: The wanted name for the html report file.
        :param genome: Path to the genome fastq file.
        :param bam2: The path to the sorted and indexed bam file of the control experiment
        :param html_alignments_dir_name: The name of the directory that contains the alignments html files,
         the directory should be located in the GUIDERNA_ALIGNMENT_RESULTS directory (default is None).
        """
        try:
            if len(icebergs_df) > 0:
                if html_alignments_dir_name is None:
                    guiderna_alignment = ''
                else:
                    guiderna_alignment = 'cut_position strand guideRNA_alignment Hamming_distance Levenshtein_distances'
                if icebergs_df.shape[0] > 1000:
                    icebergs_df = icebergs_df.iloc[:1000, :]
                # try:
                # if bam2 is None:
                #     bedfile = self.generate_bed_file(icebergs_df,
                #                                      out_dir,
                #                                      html_out_name.replace('-report.html', '-locations.bed'),
                #                                      html_alignments_dir_name)
                #
                #     cmd = settings.CREATE_REPORT_CMD_TEMPLATE.format(bedfile=bedfile,
                #                                                      genome=genome,
                #                                                      tracks=f'"{bam1}"',
                #                                                      out_file=os.path.join(out_dir, html_out_name),
                #                                                      guideRNA_alignment=guiderna_alignment,
                #                                                      html_title=html_out_name)
                #
                #     subprocess.check_call(cmd, shell=True, env=os.environ.copy())
                # else:
                bedfile = self.generate_bed_file(icebergs_df,
                                                 out_dir,
                                                 html_out_name.replace('-report.html', '-locations.bed'),
                                                 html_alignments_dir_name)

                cmd = settings.CREATE_REPORT_CMD_TEMPLATE.format(bedfile=bedfile,
                                                                 genome=genome,
                                                                 tracks=f'"{bam1}" "{bam2}"',
                                                                 out_file=os.path.join(out_dir, html_out_name),
                                                                 guideRNA_alignment=guiderna_alignment,
                                                                 html_title=html_out_name)

                subprocess.check_call(cmd, shell=True, env=os.environ.copy())
                self.add_features_to_report(out_dir / html_out_name, html_alignments_dir_name)
        except Exception:
            logging.warning(f'Failed to create HTML report - {html_out_name}, please check files paths and make sure '
                            'igv_reports is install correctly (by typing create_report in console).',
                            exc_info=False)
            sys.exit()

    def generate_reads_count_pie_chart(self,
                                       count_ca: int,
                                       count_sb: int,
                                       count_n: int,
                                       out_file_path: Path,
                                       kind: str) -> None:
        """
        Generates the read-count vs mapq scatter plot as html file (in treatment or control, depends on the inputs)

        :param count_ca: The CRISPR activities count (reads count or sites count).
        :param count_sb: The spontaneous breaks count (reads count or sites count).
        :param count_n: The noise count (reads count or sites count).
        :param out_file_path: A path to the generated html file (including the html file name as base name).
        :param kind: 'reads' or 'sites' (depends on count_ca, count_sb and count_n).
        """
        # Pie chart, where the slices will be ordered and plotted counter-clockwise:
        labels = [f'CRISPR Activities {kind} amount', f'Spontaneous Breaks {kind} amount', f'Noise {kind} count']
        values = [count_ca, count_sb, count_n]
        # explode = (0.1, 0, 0)  # only "explode" the 2nd slice (i.e. 'Hogs')
        fig = px.pie(values=values, names=labels, color=labels,
                     color_discrete_map={f'CRISPR Activities {kind} amount': 'green',
                                         f'Spontaneous Breaks {kind} amount': 'red',
                                         f'Noise {kind} count': 'blue'})
        fig.write_html(str(out_file_path))

    def generate_reads_count_mapq_scatter_plot(self,
                                               icebergs_ca_df: pd.DataFrame,
                                               icebergs_sb_df: pd.DataFrame,
                                               icebergs_noise_df: pd.DataFrame,
                                               out_file_path: Path,
                                               experiment_name: str) -> None:
        """
        Generates the read-count vs mapq scatter plot as html file.

        :param icebergs_ca_df: The CRISPR activities dataframe.
        :param icebergs_sb_df: The spontaneous breaks dataframe.
        :param icebergs_noise_df: The noise dataframe.
        :param out_file_path: A path to the generated html file (including the html file name as base name).
        :param experiment_name: The experiment name (treatment or control)
        """
        icebergs_ca_df['class'] = 'CRISPR activity'
        icebergs_sb_df['class'] = 'Spontaneous Breaks'
        icebergs_noise_df['class'] = 'Noise'
        df_all = pd.concat((icebergs_ca_df, icebergs_sb_df, icebergs_noise_df))
        fig = px.scatter(df_all, x=f'mapq {experiment_name}', y=f'read count {experiment_name}', color='class',
                         color_discrete_sequence=["green", "red", "blue"])
        fig.write_html(str(out_file_path))

    def generate_sites_amount_in_each_class(self,
                                            icebergs_ca_df: pd.DataFrame,
                                            icebergs_sb_df: pd.DataFrame,
                                            icebergs_noise_df: pd.DataFrame,
                                            out_file_path: Path) -> None:

        text = settings.HTML_SITES_AMOUNTS_TEMPLATE.format(ca_sites_count=len(icebergs_ca_df),
                                                           sb_sites_count=len(icebergs_sb_df),
                                                           n_sites_count=len(icebergs_noise_df))
        with open(out_file_path, "w") as file:
            file.write(text)

    def generate_html_main_report(self,
                                  csv_ca: Path,
                                  csv_sb: Path,
                                  csv_noise: Path,
                                  icebergs_ca_df: pd.DataFrame,
                                  icebergs_sb_df: pd.DataFrame,
                                  icebergs_noise_df: pd.DataFrame,
                                  main_report_path: Path,
                                  experiments_overview_output_folder: Path) -> None:
        """
        Generates the main report and additional html files which embedded in to the main report.

        :param csv_ca: The path to the CRISPR activities csv file.
        :param csv_sb: The path to spontaneous csv file.
        :param csv_noise: The path to noise csv file.
        :param icebergs_ca_df: The CRISPR activities dataframe.
        :param icebergs_sb_df: The spontaneous breaks dataframe.
        :param icebergs_noise_df: The noise dataframe.
        :param main_report_path: A path to the generated main report (including the report file name).
        :param experiments_overview_output_folder: A path to a directory where the
         additional embedded html files will be saved.
        """
        try:

            html_main_report = settings.HTML_MAIN_REPORT_TEMPLATE.\
                                format(csv_ca=os.path.basename(csv_ca).replace(".csv", "-report-new.html"),
                                       csv_sb=os.path.basename(csv_sb).replace(".csv", "-report-new.html"),
                                       csv_noise=os.path.basename(csv_noise).replace(".csv", "-report-new.html"))

            with codecs.open(str(main_report_path), 'w') as writefile:
                writefile.write(html_main_report)

            self.generate_reads_count_pie_chart(icebergs_ca_df.loc[:, 'read count t'].sum(),
                                                icebergs_sb_df.loc[:, 'read count t'].sum(),
                                                icebergs_noise_df.loc[:, 'read count t'].sum(),
                                                experiments_overview_output_folder / 'treatment-reads-pie-chart.html',
                                                'reads')

            self.generate_reads_count_pie_chart(icebergs_ca_df.loc[:, 'read count m'].sum(),
                                                icebergs_sb_df.loc[:, 'read count m'].sum(),
                                                icebergs_noise_df.loc[:, 'read count m'].sum(),
                                                experiments_overview_output_folder / 'control-reads-pie-chart.html',
                                                'reads')

            self.generate_reads_count_pie_chart(icebergs_ca_df[icebergs_ca_df['read count t'] > 0].shape[0],
                                                icebergs_sb_df[icebergs_sb_df['read count t'] > 0].shape[0],
                                                icebergs_noise_df[icebergs_noise_df['read count t'] > 0].shape[0],
                                                experiments_overview_output_folder / 'treatment-sites-pie-chart.html',
                                                'sites')

            self.generate_reads_count_pie_chart(icebergs_ca_df[icebergs_ca_df['read count m'] > 0].shape[0],
                                                icebergs_sb_df[icebergs_sb_df['read count m'] > 0].shape[0],
                                                icebergs_noise_df[icebergs_noise_df['read count m'] > 0].shape[0],
                                                experiments_overview_output_folder / 'control-sites-pie-chart.html',
                                                'sites')

            self.generate_reads_count_mapq_scatter_plot(icebergs_ca_df, icebergs_sb_df, icebergs_noise_df,
                                                        experiments_overview_output_folder / 'treatment-reads-count-mapq-scatter-plot.html',
                                                        't')

            self.generate_reads_count_mapq_scatter_plot(icebergs_ca_df, icebergs_sb_df, icebergs_noise_df,
                                                        experiments_overview_output_folder / 'control-reads-count-mapq-scatter-plot.html',
                                                        'm')

            self.generate_sites_amount_in_each_class(icebergs_ca_df, icebergs_sb_df, icebergs_noise_df,
                                                     experiments_overview_output_folder / 'sites-amount-in-each-class.html')

        except Exception:
            logging.warning('please check files paths and make sure igv_reports is install correctly '
                            '(by typing create_report in conlsole).',
                            exc_info=False)

    def generate_html_steps_reads_count_funnel_chart(self,
                                                     step_reads: dict,
                                                     out_dir: Path,
                                                     treatment_experiment_name: str,
                                                     control_experiment_name: str) -> None:
        """
        Generate a funnel chart of reads count in each stage of the iceberg pipeline.

        :param step_reads: list of steps and their output files total reads count (for treatment and mock).
         i.e - args['remaining_reads].
        :param out_dir: The directory where the output html file will be saved.
        :param treatment_experiment_name: A name for the treatment experiment.
        :param control_experiment_name: A name for the control experiment.
        """
        hist_t = []
        hist_m = []
        steps = []
        for (step, files_reads) in step_reads.items():
            steps.append(step)
            if len(files_reads) == 6:
                hist_t.append(files_reads[0][1] + files_reads[2][1] + files_reads[4][1])
                hist_m.append(files_reads[1][1] + files_reads[3][1] + files_reads[5][1])
            elif len(files_reads) == 4:
                hist_t.append(files_reads[0][1] + files_reads[1][1])
                hist_m.append(files_reads[2][1] + files_reads[3][1])
            elif len(files_reads) == 2:
                hist_t.append(files_reads[0][1])
                hist_m.append(files_reads[1][1])
            else:
                logging.info(files_reads)

        hist_t.extend(hist_m)
        # print(hist_t)
        logging.info(f'{hist_t}')
        logging.info(f'{steps}')
        data = dict(Quantity=hist_t,
                    Stage=steps * 2,
                    Experiment=[f'Treatment - {treatment_experiment_name}'] * len(steps) +
                               [f'Control - {control_experiment_name}'] * len(steps))

        fig = px.funnel(data, y='Stage', x='Quantity', color='Experiment',
                        color_discrete_map={f"Treatment - {treatment_experiment_name}": "#374B53",
                                            f"Control - {control_experiment_name}": "#617588"},
                        template="simple_white",
                        labels={"Stage": ""})
        fig.write_html(str(out_dir / 'reads_count.html'))

    def generate_html_command_line(self, out_dir: Path, args: dict) -> None:
        """
        Generate the user command and important argument as html file for the report.

        :param out_dir: The directory where the output html file will be saved.
        :param args: A dict of all steps internal arguments (as keys) and their values.

        """
        cmd = 'python '
        for i, word in enumerate(sys.argv):
            if i < len(sys.argv) - 1:
                cmd += word + ' '
            else:
                cmd += word

        args_str = ''
        for i, (arg, val) in enumerate(args.items()):
            if i < len(args.items()) - 1:
                args_str += arg + ': ' + f'{str(val)}<br>'
            else:
                args_str += arg + ': ' + f'{str(val)}'
        args_str = '\n'.join(self.dict_to_html(args))

        text = settings.HTML_ICEBERG_CMD_AND_ARGUMENTS_TEMPLATE.format(cmd=cmd, args=args_str)
        with codecs.open(str(out_dir / "command.html"), "w") as file:
            file.write(text)

    def dict_to_html(self, d, c=0):
        for a, b in d.items():
            yield '{}<li>{}</li>'.format('   ' * c, f'{a}:' if isinstance(b, dict) else f'{a}: {b}')
            if isinstance(b, dict):
                yield '{}<ul>\n{}\n{}</ul>'.format('   ' * c, "\n".join(self.dict_to_html(b, c + 1)), '   ' * c)


def validate_arguments(args: dict) -> None:
    """
    Validates the step arguments from the given arguments dict.

    :param args: The arguments dict (see main function description).
    """
    try:
        validation.validate_path_existence(args['curr_in_dir'])
        step_critical_path_args = ['csv_ca', 'csv_sb', 'csv_noise', 'guide_alignments_dir', 'bam1', 'bam2']
        if any(i not in args for i in step_critical_path_args):
            raise ValueError()

        for arg in step_critical_path_args:
            validation.validate_path_existence(Path(args['OUTPUT_FOLDER_PATH'], args[arg]))

    except ValueError as e:
        logging.error(f'not all necessary arguments for {STEP_NAME} step where giving.')
        raise e
    except FileNotFoundError as e:
        logging.error(f'not all files exists for {STEP_NAME} step.')
        raise e


def prepare_step(args: dict) -> None:
    """
    Prepares the step by - setting the logging configurations, setting the output directory and
    validate the step arguments.

    :param args: A dict of arguments (as keys) and their values (see main function description)
    """
    validate_arguments(args)
    logs.set_logging_config(STEP_NAME, args['OUTPUT_FOLDER_PATH'] / settings.LOGS_FOLDER_RELATIVE_PATH / STEP_LOG_FILE)
    temp = validation.set_curr_out_dir(args['OUTPUT_FOLDER_PATH'], STEP_OUT_FOLDER)
    html_path = os.path.abspath(Path(temp, 'HTML'))
    if os.path.exists(html_path):
        shutil.rmtree(html_path)
    Path(html_path).mkdir(parents=True, exist_ok=True)
    if os.path.exists(os.path.abspath(args['guide_alignments_dir'])):
        shutil.move(os.path.abspath(args['guide_alignments_dir']),
                    os.path.abspath(Path(html_path)))

    if os.path.exists(os.path.abspath(args['sites_profiles_dir'])):
        shutil.move(os.path.abspath(args['sites_profiles_dir']),
                    os.path.abspath(Path(html_path)))

    shutil.copy(settings.PROJECT_ROOT_DIRECTORY / 'images' / 'iceberg-logo.jpg', html_path)
    args['curr_out_dir'] = temp


def main(args: dict) -> None:
    """
    Run the REPORTS step on treatment and control.

    :param args: The pipeline arguments dictionary - dictionary of arguments (as keys) and their values.
    """
    try:
        rm = ReportsManager()
        pbar = logs.get_progress_bar(step_name=STEP_NAME,
                                     total=6,
                                     out_of="report")

        reports_output_dir = validation.set_curr_out_dir(Path(args['curr_out_dir'] / 'HTML'),
                                                         'Reports')

        icebergs_ca_df = csv_to_df.get_merged_experiments_df(args['OUTPUT_FOLDER_PATH'] / args['csv_ca'],
                                                             csv_contains_alignments_fields=True)

        icebergs_sb_df = csv_to_df.get_merged_experiments_df(args['OUTPUT_FOLDER_PATH'] / args['csv_sb'])

        icebergs_noise_df = csv_to_df.get_merged_experiments_df(args['OUTPUT_FOLDER_PATH'] / args['csv_noise'])

        logging.info('generating crispr activity sites report (on / off targets).')
        rm.generate_html_report(args['bam1'],
                                icebergs_ca_df,
                                reports_output_dir,
                                os.path.basename(args['csv_ca']).replace(".csv", "-report.html"),
                                args['EXPERIMENTS']['GENERAL']['REFERENCE_GENOME_PATH'],
                                args['bam2'],
                                html_alignments_dir_name='Alignments')
        pbar.update(1)

        logging.info('generating spontaneous break sites report.')
        rm.generate_html_report(args['bam1'],
                                icebergs_sb_df,
                                reports_output_dir,
                                os.path.basename(args['csv_sb']).replace(".csv", "-report.html"),
                                args['EXPERIMENTS']['GENERAL']['REFERENCE_GENOME_PATH'],
                                args['bam2'])
        pbar.update(1)

        logging.info('generating noise sites report.')
        rm.generate_html_report(args['bam1'],
                                icebergs_noise_df,
                                reports_output_dir,
                                os.path.basename(args['csv_noise']).replace(".csv", "-report.html"),
                                args['EXPERIMENTS']['GENERAL']['REFERENCE_GENOME_PATH'],
                                args['bam2'])
        pbar.update(1)

        experiments_over_view_output_dir = validation.set_curr_out_dir(args['curr_out_dir'] / 'HTML',
                                                                       'Experiment-Overview')
        rm.generate_html_main_report(args['csv_ca'],
                                     args['csv_sb'],
                                     args['csv_noise'],
                                     icebergs_ca_df,
                                     icebergs_sb_df,
                                     icebergs_noise_df,
                                     args['curr_out_dir'] / 'iceberg_report.html',
                                     experiments_over_view_output_dir)

        pbar.update(1)

        args['curr_out_dir'] = experiments_over_view_output_dir

        logging.info('generating user input file html report.')
        rm.generate_html_command_line(args['curr_out_dir'], args)

        pbar.update(1)

        rm.generate_html_steps_reads_count_funnel_chart(args['remaining_reads'],
                                                        args['curr_out_dir'],
                                                        args['EXPERIMENTS']['TX']['NAME'],
                                                        args['EXPERIMENTS']['CONTROL']['NAME'])
        pbar.update(1)

        logging.info(f'{STEP_NAME} step done.')
    except Exception as e:
        logging.error(f'Failed to run {STEP_NAME} step, please read the following for more information.',
                      exc_info=False)
        # sys.exit()
        raise e
