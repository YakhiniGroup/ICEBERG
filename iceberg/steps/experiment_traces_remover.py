"""
This Python module is used to run the EXPERIMENT_TRACES_REMOVAL step on treatment and control experiments.

For each paired-end reads, r1 and r2, from R1.fastq and R2.fastq respectively, the following steps are done -

First, the library primer that detected at the EXPERIMENT_LIBRARIES_DETECTION step is trimmed from the start of r2.
r2[:22] is trimmed if Forward primer detected and r2[:28] is trimmed if Reverse primer detected.

Second, the tag and the tag reverse complement are aligned to the start of r2 and r1, r[:32], and the
maximum alignment score is chosen, while the alignment hamming distance is less
then MAX_ALIGNMENTS_HAMMING_DISTANCE then r[:32] is trimmed, else, the process stop.

Finally, the tag and the tag reverse complement are trimmed if them appear in r2 or r1.

The second and the last steps are for the cases where the tag injected more then once and in different orientation.
"""

from utils import validation, generators, aligners, logs

import logging
from pathlib import Path


STEP_NAME = 'EXPERIMENT_TRACES_REMOVAL'
STEP_OUT_FOLDER = 'EXPERIMENT_TRACES_REMOVAL_RESULTS'
STEP_LOG_FILE = 'EXPERIMENT_TRACES_REMOVAL.log'

TRIMMED_FILE_SUFFIX = '.trimmed'


class ExperimentTracesRemover(object):
    def __init__(self,
                 tag: str,
                 tag_rc: str,
                 maximum_hamming_distance_allowed: int,
                 in_folder: Path,
                 out_folder: Path) -> None:
        """
        Initiate the ExperimentTracesRemover object.

        :param tag: The tag which inject in the cut sites.
        :param tag_rc: The reverse complement of the tag which inject in the cut sites.
        :param maximum_hamming_distance_allowed: The maximum Hamming distance allowed for alignment of the 
         tag or tag reverse complement to r2 to be valid.
        :param in_folder: The folder where the four files are stored.
        :param out_folder: The folder where the output fastq files will be saved.
        """

        self.tag = tag
        self.tag_rc = tag_rc
        self.maximum_hamming_distance_allowed = maximum_hamming_distance_allowed
        self.files_dir = in_folder
        self.out_dir = out_folder

    def align_tag_to_read_sequence_start(self,
                                         sequence: str) -> (str, int):
        """
        Align the tag and the tag reverse complement to a read sequence and return the hemming distance and the tag
        in the correct orientation for the alignment with the maximum alignment score.

        :param sequence:

        :return: A tuple contains the tag orientation with the maximum alignment score and
         the alignment Hamming distance.
        """
        aligner = aligners.get_aligner()
        logging.info(sequence)
        alignments_tag = aligner.align(sequence[:len(self.tag)], self.tag)
        alignments_tag_rc = aligner.align(sequence[:len(self.tag_rc)], self.tag_rc)

        logging.info(f'\n'
                     f'tag1 alignment score: {alignments_tag.score}\n'
                     f'{alignments_tag[0]}\n'
                     f'tag2 alignment score: {alignments_tag_rc.score}\n'
                     f'{alignments_tag_rc[0]}')

        if alignments_tag_rc.score < alignments_tag.score:
            alignments = alignments_tag
            tag = self.tag
        else:
            alignments = alignments_tag_rc
            tag = self.tag_rc

        alignment = alignments[0]
        alignment_lines = str(alignment).split('\n')
        site_sequence_aligned = alignment_lines[0]
        alignment_symbols = alignment_lines[1]
        grna_aligned = alignment_lines[2]
        hamming_distance = sum(1 for c1, c2 in zip(grna_aligned.upper(), site_sequence_aligned.upper()) if c1 != c2)
        return tag, hamming_distance

    def trim_libraries_primer(self,
                              read_r1: list,
                              read_r2: list) -> (list, list):
        """
        Trim the primer that detected at the EXPERIMENT_LIBRARIES_DETECTION step from the start of r2.
        r2[:22] is trimmed if Forward primer detected and r2[:28] is trimmed if Reverse primer detected.

        :param read_r1: read from R1.fastq.
        :param read_r2: read from R2.fastq (the pair of read_r1).

        :return: A tuple contains read_r1 and read_r2 after trimming.
        """
        library_primer = read_r2[0].split('_')[4]
        read_r2[1] = read_r2[1][len(library_primer):]
        read_r2[3] = read_r2[3][len(library_primer):]
        return read_r1, read_r2

    def trim_tag_duplicates_from_read_start(self,
                                            read: list) -> list:
        """
        The tag and the tag reverse complement are aligned to the start of r2 and r1, r[:32], and the maximum
        alignment score is chosen, while the alignment hamming distance is less than MAX_ALIGNMENTS_HAMMING_DISTANCE
        then r[:32] is trimmed, else, the process stop.

        :param read: A read from fastq file.

        :return: The read after trimming.
        """
        while len(read[1]) > 0:
            tag, hamming_distance = self.align_tag_to_read_sequence_start(read[1].upper())
            if hamming_distance > self.maximum_hamming_distance_allowed:
                return read
            read[1], read[3] = read[1][len(tag):], read[3][len(tag):]
        return read

    def naive_trim(self,
                   read: list,
                   sequence_to_trim: str) -> list:
        """
        Trims given sequence from fastq read.

        :param read: A read from fastq file.
        :param sequence_to_trim: A sequenct to trim from the read.

        return: The given read after trimming the sequence to trim.
        """
        final_read_sequence = ''
        final_read_base_qual = ''
        if sequence_to_trim in read[1]:
            read_sequence_parts = read[1].split(sequence_to_trim)
            logging.info(f'read contains {sequence_to_trim} {len(read_sequence_parts)-1} times, read header:\n'
                         f'{read[0]}')
            i = 0
            for part in read_sequence_parts:
                final_read_sequence += part
                final_read_base_qual += read[3][i:i+len(part)]
                i += len(part) + len(sequence_to_trim)
            read[1], read[3] = final_read_sequence, final_read_base_qual
        return read

    def naive_trim_tag_duplicates(self,
                                  read: list) -> list:
        """
        Trim the tag and the tag reverse complement if they appear in r1 or r2.
        For trimming tags that are not in the beginning of the read

        :param read: A read from fastq file.

        return: The given read after the tag and the tag reverse complement are trimmed.
        """
        read = self.naive_trim(read, self.tag)
        read = self.naive_trim(read, self.tag_rc)
        return read

    def trim_experiment_traces_from_pair_end_reads(self,
                                                   read_r1: list,
                                                   read_r2: list) -> (list, list):
        """
        Trim the experiments traces from r1 and r2, in R1.fastq and R2.fastq respectively

        First, the library primer that detected at the EXPERIMENT_LIBRARIES_DETECTION step is trimmed from the start
        of r2. r2[:22] is trimmed if Forward primer detected and r2[:28] is trimmed if Reverse primer detected.

        Second, the tag and the tag reverse complement are aligned to the start of r2 and r1, r[:32], and the maximum
        alignment score is chosen, while the alignment hamming distance is less than MAX_ALIGNMENTS_HAMMING_DISTANCE
        then r[:32] is trimmed, else, we stop.

        Finally, the tag and the tag reverse complement are trimmed if them appear in r2 or r1.

        :param read_r1: read from R1.fastq.
        :param read_r2: read from R2.fastq (the pair of read_r1).

        return: A tuple contains read_r1 and read_r2 after trimming.

        """
        read_r1_trimmed, read_r2_trimmed = self.trim_libraries_primer(read_r1,
                                                                          read_r2)

        read_r1_trimmed = self.trim_tag_duplicates_from_read_start(read_r1_trimmed)
        read_r2_trimmed = self.trim_tag_duplicates_from_read_start(read_r2_trimmed)

        return self.naive_trim_tag_duplicates(read_r1_trimmed), self.naive_trim_tag_duplicates(read_r2_trimmed)

    def trim_experiment_traces_from_pair_end_fastqs(self,
                                                    r1_file_name: str,
                                                    r2_file_name: str,
                                                    experiment_name: str) -> (str, str):
        """
        Trim the experiments traces from pair end fastq files.

        :param r1_file_name: The fastq file name of R1.
        :param r2_file_name: The fastq file name of R2.
        :param experiment_name: The experiment name.

        :return: Two fastq files names - file for R1 and R2 after trimming.
        """
        logging.info(f'trims experiment traces from {r1_file_name} and {r2_file_name}, '
                     'maximum hamming distance allowed for alignment to be valid is: '
                     f'{self.maximum_hamming_distance_allowed})...')
        r1_tag = open(self.out_dir / r1_file_name.replace('.fastq', TRIMMED_FILE_SUFFIX + '.fastq'), 'w')
        r2_tag = open(self.out_dir / r2_file_name.replace('.fastq', TRIMMED_FILE_SUFFIX + '.fastq'), 'w')

        experiment_reads = zip(generators.get_fastq_read_lines(self.files_dir / r1_file_name),
                               generators.get_fastq_read_lines(self.files_dir / r2_file_name))

        pbar = logs.get_progress_bar(step_name=STEP_NAME,
                                     experiment_name=experiment_name,
                                     total=sum(1 for _ in generators.get_fastq_read_lines(self.files_dir / r1_file_name)),
                                     out_of='paired-end reads')

        for read_r1, read_r2 in experiment_reads:
            read_r1_trimmed, read_r2_trimmed = self.trim_experiment_traces_from_pair_end_reads(read_r1, read_r2)
            if len(read_r1_trimmed[1]) > 0 and len(read_r2_trimmed[1]) > 0:
                for line_r1, line_r2 in zip(read_r1_trimmed, read_r2_trimmed):
                    r1_tag.write(line_r1 + '\n')
                    r2_tag.write(line_r2 + '\n')
            pbar.update(1)
        r1_tag.close()
        r2_tag.close()
        pbar.close()
        return (r1_file_name.replace('.fastq', TRIMMED_FILE_SUFFIX + '.fastq'),
                r2_file_name.replace('.fastq', TRIMMED_FILE_SUFFIX + '.fastq'))


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
    logs.set_logging_config(STEP_NAME, Path(args['OUTPUT_FOLDER_PATH'], 'DEBUG', 'LOGS', STEP_LOG_FILE))
    args['curr_out_dir'] = validation.set_curr_out_dir(args['OUTPUT_FOLDER_PATH'], STEP_OUT_FOLDER)
    validate_arguments(args)


def main(args: dict) -> None:
    """
     Run the EXPERIMENT_TRACES_REMOVAL step on treatment and control.

    :param args: The pipeline arguments dictionary - dictionary of arguments (as keys) and their values.
    """
    try:
        experiments = args['EXPERIMENTS']
        th = ExperimentTracesRemover(experiments['GENERAL']['EXPERIMENTS_TAG'],
                                     experiments['GENERAL']['EXPERIMENTS_TAG_RC'],
                                     args['HYPERPARAMATERS']['MAX_ALIGNMENTS_HAMMING_DISTANCE'],
                                     args['curr_in_dir'],
                                     args['curr_out_dir'])

        treatment, control = experiments['TX'], experiments['CONTROL']
        treatment['R1'], treatment['R2'] = th.trim_experiment_traces_from_pair_end_fastqs(treatment['R1'],
                                                                                          treatment['R2'],
                                                                                          treatment['NAME'])

        control['R1'], control['R2'] = th.trim_experiment_traces_from_pair_end_fastqs(control['R1'],
                                                                                      control['R2'],
                                                                                      control['NAME'])

        args['curr_in_dir'] = validation.set_curr_in_dir(args['OUTPUT_FOLDER_PATH'], STEP_OUT_FOLDER)
        logging.info(f'{STEP_NAME} done.')
    except Exception as e:
        logging.error(f'Failed to run {STEP_NAME} step, please read the following for more information.',
                      exc_info=False)
        # sys.exit()
        raise e
