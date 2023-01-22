"""
This Python module is used to run the EXPERIMENT_LIBRARIES_DETECTION step on treatment and control experiments.

Each paired-end reads, r1 and r2, from R1.fastq and R2.fastq respectively, is from one of the following libraries:
Forward library: paired-end reads where part of the experiment tag reverse complement is at the start of r2.
Reverse library: paired-end reads where part of the experiment tag is at the start of r2.

This happens due the PCR rounds for amplify sites containing the experiment injected tag, The PCR primers are -
Forward primer: tag_reverse_complement[12:] (22 nt).
Reverse primer: tag[6:] (28 nt).

The library of each paired-end reads is detected by aligning the libraries PCR-primers
to the start of r2, r2[:22] or r2[:28] respectively.

The library of the primer with the maximum alignment-score is the one detected for both r1 and r2.
The detected PCR primer with the alignment Hamming distance are appended to the (identical) umi section
in r1 and r2 fastq headers.
"""

from utils import validation, generators, aligners, logs
import settings

import re
import logging
from pathlib import Path

STEP_NAME = 'EXPERIMENT_LIBRARIES_DETECTION'
STEP_OUT_FOLDER = 'EXPERIMENT_LIBRARIES_DETECTION_RESULTS'
STEP_LOG_FILE = 'EXPERIMENT_LIBRARIES_DETECTION.log'

DETECTED_FILE_SUFFIX = '.detected'


class ExperimentLibrariesDetector(object):
    def __init__(self,
                 tag: str,
                 tag_rc: str,
                 in_folder: Path,
                 out_folder: Path) -> None:
        """
        Initiate the ExperimentLibrariesDetector object.

        :param tag (str): The tag which inject in the cut sites.
        :param tag_rc: The reverse complement of the tag which inject in the cut sites.
        :param in_folder: The folder where the four files are stored.
        :param out_folder: The folder where the output fastq files will be saved.

        :return: a list contains four lines which represent one read.
        """
        self.forward_primer = tag_rc[12:]  # forward primer (22nt)
        self.reverse_primer = tag[6:]  # reverse primer (28nt)
        self.in_folder = in_folder
        self.out_folder = out_folder

    def align_primer_to_read_sequence_start(self,
                                            sequence: str,
                                            read_name: str) -> tuple:
        """
        Align the libraries PCR-primers to the start of r2, r2[:22] or r2[:28] respectively.
        The library of the primer with the maximum alignment-score is the one detected for both r1 and r2.

        :param sequence: The read sequence.
        :param read_name: The read name.

        :return: A tuple contains the primer with the maximum alignment score and the alignment Hamming distance,
         returns tuple contains none and len of primer if the alignment fails.
        """
        aligner = aligners.get_aligner()

        alignments_forward_primer = aligner.align(sequence[:len(self.forward_primer)],
                                                  self.forward_primer)

        alignments_reverse_primer = aligner.align(sequence[:len(self.reverse_primer)],
                                                  self.reverse_primer)

        logging.info(settings.TAG_AND_TAG_RC_COMPARISON_TEMPLATE.
                     format(read_name=read_name,
                            alignments_forward_score=alignments_forward_primer.score,
                            alignments_forward=alignments_forward_primer[0],
                            alignments_reverse_score=alignments_reverse_primer.score,
                            alignments_reverse=alignments_reverse_primer[0]))

        if alignments_reverse_primer.score < alignments_forward_primer.score:
            alignments = alignments_forward_primer
            primer = self.forward_primer

        else:
            alignments = alignments_reverse_primer
            primer = self.reverse_primer

        if len(alignments) == 0:
            return 'none', len(primer)

        alignment = alignments[0]

        alignment_lines = str(alignment).split('\n')
        site_sequence_aligned = alignment_lines[0]
        alignment_symbols = alignment_lines[1]
        grna_aligned = alignment_lines[2]
        hamming_distance = sum(1 for c1, c2 in zip(grna_aligned.upper(), site_sequence_aligned.upper()) if c1 != c2)
        return primer, hamming_distance

    def detect_reads_experiment_libraries(self,
                                          r1_file_name: str,
                                          r2_file_name: str,
                                          experiment_name: str) -> tuple:
        """
        Detect the experiment-library (forward-library or reverse-library) of each paired-end reads, r1 and r2.
        The primer of the library with the highest score is added to the read header in the new fastq file.

        :param r1_file_name: The fastq file name of r1.
        :param r2_file_name: The fastq file name of r2.
        :param experiment_name: The experiment name.

        :return: Two files names - file for R1.fastq and file for R1.fastq after the experiment libraries detection.
        """
        logging.info(f'detect the library of each paired-end read in {r1_file_name} and {r2_file_name}...')

        pbar = logs.get_progress_bar(step_name=STEP_NAME,
                                     experiment_name=experiment_name,
                                     total=sum(1 for _ in generators.get_fastq_read_lines(self.in_folder / r1_file_name)),
                                     out_of="paired-end reads")

        r1_detected = open(self.out_folder / r1_file_name.replace('.fastq', DETECTED_FILE_SUFFIX + '.fastq'), 'w')
        r2_detected = open(self.out_folder / r2_file_name.replace('.fastq', DETECTED_FILE_SUFFIX + '.fastq'), 'w')

        count_reverse_primer, count_forward_primer, count_both, total_count = 0, 0, 0, 0

        experiment_reads = zip(generators.get_fastq_read_lines(self.in_folder / r1_file_name),
                               generators.get_fastq_read_lines(self.in_folder / r2_file_name))

        for read_r1, read_r2 in experiment_reads:

            primer, hamming_distance = self.align_primer_to_read_sequence_start(read_r2[1],
                                                                                read_r2[0])
            if primer == self.forward_primer:
                count_forward_primer += 1
            elif primer == self.reverse_primer:
                count_forward_primer += 1

            for i, (line_r1, line_r2) in enumerate(zip(read_r1, read_r2)):
                if i == 0:
                    umi_regex = '(@[a-zA-Z0-9_]*)'
                    r1_detected.write(re.sub(umi_regex, r'\1_' + f'{primer}_{hamming_distance}', line_r1) + '\n')
                    r2_detected.write(re.sub(umi_regex, r'\1_' + f'{primer}_{hamming_distance}', line_r2) + '\n')
                else:
                    r1_detected.write(line_r1 + '\n')
                    r2_detected.write(line_r2 + '\n')
            pbar.update(1)

        r1_detected.close()
        r2_detected.close()

        logging.info(f'{self.forward_primer}: {count_forward_primer}'
                     f'{self.reverse_primer}: {count_reverse_primer}'
                     f'total:{total_count}')
        pbar.close()
        return (r1_file_name.replace('.fastq', DETECTED_FILE_SUFFIX + '.fastq'),
                r2_file_name.replace('.fastq', DETECTED_FILE_SUFFIX + '.fastq'))


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


def main(args: dict) -> None:
    """
    Run the EXPERIMENT_LIBRARIES_DETECTION step on treatment and control.

    :param args: The pipeline arguments dictionary = dictionary of arguments (as keys) and their values.
    """
    try:
        experiments = args['EXPERIMENTS']
        eld = ExperimentLibrariesDetector(experiments['GENERAL']['EXPERIMENTS_TAG'],
                                          experiments['GENERAL']['EXPERIMENTS_TAG_RC'],
                                          args['curr_in_dir'],
                                          args['curr_out_dir'])

        treatment, control = experiments['TX'], experiments['CONTROL']
        treatment['R1'], treatment['R2'] = eld.detect_reads_experiment_libraries(treatment['R1'],
                                                                                 treatment['R2'],
                                                                                 treatment['NAME'])

        control['R1'], control['R2'] = eld.detect_reads_experiment_libraries(control['R1'],
                                                                             control['R2'],
                                                                             control['NAME'])

        args['curr_in_dir'] = validation.set_curr_in_dir(args['OUTPUT_FOLDER_PATH'], STEP_OUT_FOLDER)

        logging.info(f'{STEP_NAME} done.')
    except Exception as e:
        logging.error(f'Failed to run {STEP_NAME}, please read the following for more information.',
                      exc_info=False)
        # sys.exit()
        raise e
