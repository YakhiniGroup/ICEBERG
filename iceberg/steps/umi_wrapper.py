"""
This Python module is used to run the UMI step (Unigue Molecular Index by GuideSeq) on treatment and
control experiments.

Additional option  - the user can choose the amount of bps from r1 and r2 start that will help to build
the umi for the paired end reads (GuideSeq use 6 as fixed number).
"""

from utils import logs, validation
from resources.umi import umitag, consolidate

import logging
from pathlib import Path

STEP_NAME = 'UMI'
STEP_OUT_FOLDER = 'UMI_RESULTS'
STEP_LOG_FILE = 'UMI.log'

CONSOLIDATE_FILE_SUFFIX = '.consolidated'
UMI_FILE_SUFFIX = '.umi'


class UMIWrapper(object):
    def __init__(self,
                 min_quality: int,
                 min_frequency: float,
                 umi_bps_amount_from_reads_start: int,
                 in_folder: Path,
                 out_folder: Path):
        """
         Wrapper that run the umitag step followed by the consolidate step of the UMI package.

        :param min_quality: The minimum base quality of bp in a read for it to be consolidated in the consolidation.
        :param min_frequency: The minimum base frequency of bp in a read for it to be consolidated in the consolidation.
        :param in_folder: The folder where the four files are stored.
        :param out_folder: The folder where the output fastq files will be saved.
        """
        self.min_quality = min_quality
        self.min_frequency = min_frequency
        self.umi_bps_amount_from_reads_start = umi_bps_amount_from_reads_start
        self.in_folder = in_folder
        self.out_folder = out_folder

    def identify_reads(self,
                       r1_file_name: str,
                       r2_file_name: str,
                       i1_file_name: str,
                       i2_file_name: str) -> (str, str):
        """
        Given paired-end reads from experiment, R1 and R2, with their indexes, I1 and I2,
        the function assign to each paired-end reads an Unique Molecular Index (UMI)

        :param r1_file_name: The fastq file name of R1.
        :param r2_file_name: The fastq file name of I1.
        :param i1_file_name: The index fastq file name of R2.
        :param i2_file_name: The index fastq file name of I2.

        :return: the names of the umi files.
        """
        try:
            r1_umi_file_name = r1_file_name.replace('.fastq', UMI_FILE_SUFFIX + '.fastq')
            r2_umi_file_name = r2_file_name.replace('.fastq', UMI_FILE_SUFFIX + '.fastq')
            umitag.umitag(str(self.in_folder / r1_file_name),
                          str(self.in_folder / r2_file_name),
                          str(self.in_folder / i1_file_name),
                          str(self.in_folder / i2_file_name),
                          str(self.out_folder / r1_umi_file_name),
                          str(self.out_folder / r2_umi_file_name),
                          str(self.out_folder),
                          self.umi_bps_amount_from_reads_start)

            return r1_umi_file_name, r2_umi_file_name
        except Exception as e:
            logging.error(f'Failed to run umitag step on files: {r1_file_name} and {r1_file_name},'
                          f' please check files paths.',
                          exc_info=False)
            raise e

    def consolidate_reads(self,
                          umi_file_name: str) -> str:
        """
        Given one of the fastq files from paired-end reads experiment (R1 or R2), the function consolidates
        the reads by their Unique Molecular Index (UMI), which been given to each paired-end reads in the umitag step.

        :param umi_file_name: fastq file name (after umitag).
        """
        try:
            consolidated_file_name = umi_file_name.replace('.fastq', CONSOLIDATE_FILE_SUFFIX + '.fastq')
            consolidate.consolidate(str(self.out_folder / umi_file_name),
                                    str(self.out_folder / consolidated_file_name),
                                    self.min_quality,
                                    self.min_frequency)

            return consolidated_file_name
        except Exception as e:
            logging.error(f'Failed to run umi consolidate step on files: {umi_file_name} please check file path.',
                          exc_info=False)
            # sys.exit()
            raise e

    def collapse(self,
                 r1_file_name: str,
                 r2_file_name: str,
                 i1_file_name: str,
                 i2_file_name: str,
                 experiment_name: str) -> (str, str):
        """
        Given paired-end reads from experiment, R1 and R2, with their indexes, I1 and I2,
        the function first assign to each paired-end reads an Unique Molecular Index (UMI) and then
        consolidate the reads by it.

        :param r1_file_name: The fastq file name of R1.
        :param r2_file_name: The fastq file name of R2.
        :param i1_file_name: The index fastq file name of I1.
        :param i2_file_name: The index fastq file name of I2.
        :param experiment_name: The experiment name.

        :return: The names of the fastq files after collapse (umitag and consolidated).
        """
        pbar = logs.get_progress_bar(step_name=STEP_NAME,
                                     experiment_name=experiment_name,
                                     total=3,
                                     out_of="step")

        logging.info("assign unique molecular identifier (umi) to each read header in R1 and R2 fastq files...")
        r1_umi_file_name, r2_umi_file_name = self.identify_reads(r1_file_name,
                                                                 r2_file_name,
                                                                 i1_file_name,
                                                                 i2_file_name)
        pbar.update(1)
        logging.info(f"consolidate reads by unique molecular identifier (umi collapse) in {r1_umi_file_name}...")
        r1_consolidated_file_name = self.consolidate_reads(r1_umi_file_name)
        pbar.update(1)

        logging.info(f"consolidate reads by unique molecular identifier (umi collapse) in {r2_umi_file_name}...")
        r2_consolidated_file_name = self.consolidate_reads(r2_umi_file_name)
        pbar.update(1)

        return r1_consolidated_file_name, r2_consolidated_file_name


def validate_arguments(args: dict) -> None:
    """
    Validates arguments for the step that may changed/created during the pipeline from the given arguments dict.

    :param args: The arguments dict (see main function description).
    """
    try:
        # validation.validate_path_existence(args['curr_in_dir'])
        files_folder = args['curr_in_dir'] if args['curr_in_dir'] == "" else Path(args['curr_in_dir'])

        for file in [args['EXPERIMENTS']['TX']['R1'], args['EXPERIMENTS']['TX']['I1'],
                     args['EXPERIMENTS']['TX']['R2'], args['EXPERIMENTS']['TX']['I2']]:
            if files_folder == "":
                validation.validate_path_existence(args['EXPERIMENTS']['TX']['EXPERIMENT_FOLDER_PATH'] / file)
            else:
                validation.validate_path_existence(Path(args['curr_in_dir'], file))

        for file in [args['EXPERIMENTS']['CONTROL']['R1'], args['EXPERIMENTS']['CONTROL']['I1'],
                     args['EXPERIMENTS']['CONTROL']['R2'], args['EXPERIMENTS']['CONTROL']['I2']]:

            if files_folder == "":
                validation.validate_path_existence(args['EXPERIMENTS']['CONTROL']['EXPERIMENT_FOLDER_PATH'] / file)
            else:
                validation.validate_path_existence(Path(args['curr_in_dir'], file))

    except ValueError as e:
        logging.error('Failed to validate arguments for umi step.')
        raise e


def prepare_step(args: dict) -> None:
    """
    Prepares the step by - setting the logging configurations, setting the output directory and
    validate the step arguments.

    :param args: A dict of arguments (as keys) and their values (see main function description)
    """
    logs.set_logging_config('umi', Path(args['OUTPUT_FOLDER_PATH'], 'DEBUG', 'LOGS', 'umi.log'))
    args['curr_out_dir'] = validation.set_curr_out_dir(args['OUTPUT_FOLDER_PATH'], STEP_OUT_FOLDER)
    validate_arguments(args)


def main(args) -> None:
    """
    Run the UMI step on treatment and control.

    :param args: The pipeline arguments dictionary - dictionary of arguments (as keys) and their values.
    """
    try:
        umiw = UMIWrapper(args['HYPERPARAMATERS']['MIN_QUALITY'],
                          args['HYPERPARAMATERS']['MIN_FREQUENCY'],
                          args['HYPERPARAMATERS']['UMI_BPS_AMOUNT_FROM_READS_START'],
                          args['EXPERIMENTS']['TX']['EXPERIMENT_FOLDER_PATH'],
                          args['curr_out_dir'])

        treatment = args['EXPERIMENTS']['TX']
        treatment['R1'], treatment['R2'] = umiw.collapse(treatment['R1'],
                                                         treatment['R2'],
                                                         treatment['I1'],
                                                         treatment['I2'],
                                                         treatment['NAME'])


        umiw.in_folder = args['EXPERIMENTS']['CONTROL']['EXPERIMENT_FOLDER_PATH']

        control = args['EXPERIMENTS']['CONTROL']
        control['R1'], control['R2'] = umiw.collapse(control['R1'],
                                                     control['R2'],
                                                     control['I1'],
                                                     control['I2'],
                                                     control['NAME'])

        args['curr_in_dir'] = validation.set_curr_in_dir(args['OUTPUT_FOLDER_PATH'], STEP_OUT_FOLDER)
        logging.info('umi step done.')
    except Exception as e:
        logging.error('error occurred while running umi step, please read the following for more information.',
                      exc_info=False)
        # sys.exit()
        raise e