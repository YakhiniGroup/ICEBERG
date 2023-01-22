"""
This Python module is used to run the BWA step on treatment and control experiments.

The step flow is -

Pre-processing:
1. The genome file indexed if needed (using bwa).
2. The fastq files are indexed (using bwa).

Processing:
1. The fastq files are aligned to the genome and generates a sam file (using bwa).

Post-processing:
4. The sam file are sorted by chromosomes where the aligned reads are sort in each chromosome (using samtools).
5. The sam file is converted to bam file (using samtools).
6. The bam file is indexed for later usages in igv browser (using samtools).
"""

import settings
from utils import validation, logs

import subprocess
import os
import logging
from pathlib import Path

STEP_NAME = 'BWA'
STEP_OUT_FOLDER = 'BWA_RESULTS'
STEP_LOG_FILE = 'BWA.log'

SORTED_SUFFIX = '.sorted'
BAM_SUFFIX = '.bam'
SAM_SUFFIX = '.sam'


class BWAWrapper:

    def __init__(self,
                 in_folder: Path,
                 out_folder: Path):
        """
        Initiate the BWAWrapper object.

        :param in_folder: The folder where the four files are stored.
        :param out_folder: The folder where the output sam and bam files will be saved.
        """
        self.in_folder = in_folder
        self.out_folder = out_folder

    def index_genome(self,
                     genome_ref: Path) -> None:
        """
        Index genome using bwa index command.

        :param genome_ref: The path to the genome file.
        """
        try:
            cmd = settings.INDEX_GENOME_CMD_TEMPLATE.format(genome_ref=str(genome_ref))
            proc = subprocess.Popen(cmd,
                                    shell=True,
                                    env=os.environ.copy(),
                                    stdout=subprocess.PIPE,
                                    universal_newlines=True,
                                    stderr=subprocess.STDOUT)

            logs.log_subprocess(proc)
        except Exception as e:
            logging.error(f'Failed to index genome - {genome_ref}.'
                          'please check genome path and make sure bwa is install correctly '
                          '(by typing bwa in console).',
                          exc_info=False)
            # sys.exit()
            raise e

    def index_fastq(self,
                    file: str) -> None:
        """
        Indexes the fastq file using bwa index command.

        :param file: The path to the directory that contain the fastq file.
        """
        try:
            cmd = settings.INDEX_FASTQ_CMD_TEMPLATE.format(file_path=str(self.in_folder / file))

            proc = subprocess.Popen(cmd,
                                    shell=True,
                                    env=os.environ.copy(),
                                    stdout=subprocess.PIPE,
                                    universal_newlines=True,
                                    stderr=subprocess.STDOUT)
            logs.log_subprocess(proc)
        except Exception as e:
            logging.error(f'Failed to index fastq file - {file}.'
                          'please check fastq file path and make sure bwa is install correctly '
                          '(by typing bwa in console).',
                          exc_info=False)
            # sys.exit()
            raise e

    def align_pair_end_reads_bwa(self,
                                 genome_ref: Path,
                                 r1_file_name: str,
                                 r2_file_name: str,
                                 out_sam_name: str) -> None:
        """
        Aligns the pair end reads to the genome using bwa mem command.

        :param genome_ref: The path to the genome file.
        :param r1_file_name: The fastq file name of r1.
        :param r2_file_name: The fastq file name of r2.
        :param out_sam_name: The wanted name for the output sam file.
        """
        try:
            cmd = settings.BWA_CMD_TEMPLATE.format(genome_ref=str(genome_ref),
                                                   r1_fastq_path=str(self.in_folder / r1_file_name),
                                                   r2_fastq_path=str(self.in_folder / r2_file_name),
                                                   out_sam_path=str(self.out_folder / out_sam_name))
            proc = subprocess.Popen(cmd,
                                    shell=True,
                                    env=os.environ.copy(),
                                    stdout=subprocess.PIPE,
                                    universal_newlines=True,
                                    stderr=subprocess.STDOUT)
            logs.log_subprocess(proc)
        except Exception as e:
            logging.error(f'Failed to align paired end reads to genome - {genome_ref} please check files paths and '
                          'make sure bwa is install correctly (by typing bwa in console) and '
                          'validate your ram has at least 5.5GB free.', exc_info=False)
            # sys.exit()
            raise e

    def sort_sam(self,
                 sam_name: str) -> str:
        """
        Sorts the sam file using samtool sort command.

        :param sam_name: The name of the sam file.

        :return: The sorted sam file name.
        """
        try:
            sorted_out_sam_name = sam_name.replace(SAM_SUFFIX, SORTED_SUFFIX + SAM_SUFFIX)

            cmd = settings.SORT_CMD_TEMPLATE.format(out_format='sam',
                                                    in_file_path=self.out_folder / sam_name,
                                                    out_file_path=self.out_folder / sorted_out_sam_name)

            proc = subprocess.Popen(cmd,
                                    shell=True,
                                    env=os.environ.copy(),
                                    stdout=subprocess.PIPE,
                                    universal_newlines=True,
                                    stderr=subprocess.STDOUT)

            logs.log_subprocess(proc)
            return sorted_out_sam_name
        except Exception as e:
            logging.error(f'Failed to sort igvtools sam file - {sam_name}. '
                          f'please check sam file path and make sure igvtools is install correctly '
                          f'(by typing igvtools in console).',
                          exc_info=False)
            # sys.exit()
            raise e

    def sam_to_bam(self,
                   sam_name: str) -> str:
        """
        Converts the sam file to bam file using samtool view command.

        :param sam_name: The name of the sam file.

        :return: the new bam file name
        """
        try:
            out_bam_name = sam_name.replace(SAM_SUFFIX, BAM_SUFFIX)

            cmd = settings.VIEW_SAM_TO_BAM_CMD_TEMPLATE.format(sam_path=self.out_folder / sam_name,
                                                               bam_path=self.out_folder / out_bam_name)

            proc = subprocess.Popen(cmd,
                                    shell=True,
                                    env=os.environ.copy(),
                                    stdout=subprocess.PIPE,
                                    universal_newlines=True,
                                    stderr=subprocess.STDOUT)

            logs.log_subprocess(proc)
            return out_bam_name
        except Exception as e:
            logging.error(f'Failed to convert sam file to bam file - {sam_name}.'
                          f' please check sam file path and make sure samtools is install correctly '
                          f'(by typing samtools in console).',
                          exc_info=False)
            # sys.exit()
            raise e

    def sort_bam(self,
                 bam_name: str) -> str:
        """
        Sorts the bam file using samtool sort command.

        :param bam_name: The name of the bam file.

        :return: the sorted bam file name.
        """
        try:
            sorted_out_bam_name = bam_name.replace(BAM_SUFFIX, SORTED_SUFFIX + BAM_SUFFIX)

            cmd = settings.SORT_CMD_TEMPLATE.format(out_format='bam',
                                                    in_file_path=self.out_folder / bam_name,
                                                    out_file_path=self.out_folder / sorted_out_bam_name)

            proc = subprocess.Popen(cmd,
                                    shell=True,
                                    env=os.environ.copy(),
                                    stdout=subprocess.PIPE,
                                    universal_newlines=True,
                                    stderr=subprocess.STDOUT)

            logs.log_subprocess(proc)
            return sorted_out_bam_name
        except Exception as e:
            logging.error(f'Failed to sort bam file - {bam_name}. '
                          'please check bam file path and make sure samtools is install correctly '
                          '(by typing samtools in console).',
                          exc_info=False)
            # sys.exit()
            raise e

    def index_bam(self,
                  bam_name: str) -> None:
        """
        Indexes the bam file using samtools index command.

        :param bam_name: The name of the bam file.
        """
        try:
            cmd = settings.INDEX_BAM_CMD_TEMPLATE.format(bam_path=self.out_folder / bam_name)

            proc = subprocess.Popen(cmd,
                                    shell=True,
                                    env=os.environ.copy(),
                                    stdout=subprocess.PIPE,
                                    universal_newlines=True,
                                    stderr=subprocess.STDOUT)
            logs.log_subprocess(proc)
        except Exception as e:
            logging.error(f'Failed to index bam file - {bam_name}. '
                          f'Please check bam file path and make sure samtools is install correctly '
                          f'(by typing samtools in console).',
                          exc_info=False)
            # sys.exit()
            raise e

    def align_pair_end_reads_to_genome(self,
                                       r1_file_name: str,
                                       r2_file_name: str,
                                       genome_ref: Path,
                                       out_sam_name: str,
                                       experiment_name: str) -> (str, Path):
        """
        Align pair end fastq file to a given genome file using the following flow:

        Pre-processing:
        indexes the genome (if needed) and indexes the experiment fastq files.

        Processing:
        Aligns the pre-processed paired-end fastq files to the genome.

        Post-processing:
        sorts the Aligned sam file, generate bam file (from the sorted sam file) and
        index the sorted bam file.

        :param genome_ref: The path to the genome file.
        :param r1_file_name: The fastq file name of r1.
        :param r2_file_name: The fastq file name of r2.
        :param out_sam_name: The wanted name for the output sam file.
        :param experiment_name: The experiment name.

        """
        try:
            pbar = logs.get_progress_bar(step_name=STEP_NAME,
                                         experiment_name=experiment_name,
                                         total=8,
                                         out_of='step')

            genome_index_files_extensions = ['.pac', '.amb', '.ann', '.bwt', '.sa', '.fai']
            if not all(Path(str(genome_ref) + extension).exists() for extension in genome_index_files_extensions):
                logging.info(f'index genome reference {genome_ref} for bwa...')
                self.index_genome(genome_ref)
            # for extension in genome_index_files_extensions:
            #     if not Path(str(genome_ref) + extension).exists():
            #         logging.info(f'index genome reference {genome_ref} for bwa...')
            #         self.index_genome(genome_ref)
            #         break
            pbar.update(1)

            logging.info(f'index fastq file {r1_file_name} for bwa...')
            self.index_fastq(r1_file_name)
            pbar.update(1)

            logging.info(f'index fastq file {r2_file_name} for bwa...')
            self.index_fastq(r2_file_name)
            pbar.update(1)

            logging.info(f'align reads in {r1_file_name} and {r2_file_name} to genome reference '
                         f'and save as sam file...')
            self.align_pair_end_reads_bwa(genome_ref,
                                          r1_file_name,
                                          r2_file_name,
                                          out_sam_name)
            pbar.update(1)

            logging.info(f'sort the aligned reads in sam file {out_sam_name}...')
            sorted_out_sam_name = self.sort_sam(out_sam_name)
            pbar.update(1)

            logging.info(f'generate bam file from sam file {sorted_out_sam_name}...')
            out_bam_name = self.sam_to_bam(out_sam_name)
            pbar.update(1)

            logging.info(f'sorting bam file {out_bam_name}...')
            sorted_out_bam_name = self.sort_bam(out_bam_name)
            pbar.update(1)

            logging.info(f'index the sorted bam file {sorted_out_bam_name}...')
            self.index_bam(sorted_out_bam_name)
            pbar.update(1)

            return sorted_out_sam_name, self.out_folder / sorted_out_bam_name
        except Exception as e:
            logging.exception("Exception occurred, please try to sort and index the pair-end files and "
                              "the genome file again.",
                              exc_info=False)
            raise e


def validate_arguments(args: dict) -> None:
    """
    Validates the step arguments from the given arguments dict.

    :param args: The arguments dict (see main function description).
    """
    try:
        files_folder = args['curr_in_dir'] if args['curr_in_dir'] == "" else Path(args['curr_in_dir'])

        for file in [args['EXPERIMENTS']['TX']['R1'], args['EXPERIMENTS']['TX']['R2']]:
            if files_folder == "":
                validation.validate_path_existence(args['EXPERIMENTS']['TX']['EXPERIMENT_FOLDER_PATH'] / file)
            else:
                validation.validate_path_existence(Path(args['curr_in_dir'], file))

        for file in [args['EXPERIMENTS']['CONTROL']['R1'], args['EXPERIMENTS']['CONTROL']['R2']]:
            if files_folder == "":
                validation.validate_path_existence(args['EXPERIMENTS']['CONTROL']['EXPERIMENT_FOLDER_PATH'] / file)
            else:
                validation.validate_path_existence(Path(args['curr_in_dir'], file))
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


def main(args):
    """
     Run the BWA step on treatment and control.

    :param args: The pipeline arguments dictionary - dictionary of arguments (as keys) and their values.
    """
    try:
        bwaw = BWAWrapper(args['curr_in_dir'],
                          args['curr_out_dir'])

        experiments = args['EXPERIMENTS']
        treatment, control = experiments['TX'], experiments['CONTROL']

        args['sam1'], args['bam1'] = bwaw.align_pair_end_reads_to_genome(treatment['R1'],
                                                                         treatment['R2'],
                                                                         experiments['GENERAL']['REFERENCE_GENOME_PATH'],
                                                                         '{0}.sam'.format(treatment['NAME']),
                                                                         treatment['NAME'])

        args['sam2'], args['bam2'] = bwaw.align_pair_end_reads_to_genome(control['R1'],
                                                                         control['R2'],
                                                                         experiments['GENERAL']['REFERENCE_GENOME_PATH'],
                                                                         '{0}.sam'.format(control['NAME']),
                                                                         control['NAME'])

        args['curr_in_dir'] = validation.set_curr_in_dir(args['OUTPUT_FOLDER_PATH'], STEP_OUT_FOLDER)
        logging.info(f'{STEP_NAME} step done.')
    except Exception as e:
        logging.error(f'Failed to run {STEP_NAME} step, please read the following for more information.',
                      exc_info=False)
        # sys.exit()
        raise e
