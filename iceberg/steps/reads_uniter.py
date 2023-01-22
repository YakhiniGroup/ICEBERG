"""
This Python module is used to run the UNITE_READS_TO_ICEBERGS step on treatment and control experiments.

Iterates over the sorted aligned bam file reads and unite them to icebergs, where the output
is a csv file where each row represent a single iceberg.

If read distance from current iceberg (where iceberg initiate as one read) is less then MAX_READS_DISTANCE
than the read is added to the current iceberg.

Each iceberg contains information about:
1. The iceberg chromosome.
2. The iceberg start index.
3. The iceberg end index.
4. The iceberg read count.
4. The iceberg valid pair-end reads depth.
5. The iceberg length.
6. The iceberg mapq - The median mapq of the iceberg reads.
7. The paired end reads orientation - Amount of pair end reads from the right of the tag and from the left of the tag.
8. The cut position candidate with information about the experiment libraries.
9. The iceberg depth (amount of reads in highest candidate from both forward and reverse libraries).

!Most of the fields allowed only paired-end reads with reasonable distance, non supplementary and non secondary reads!
"""

from utils import generators, logs, validation

from statistics import median
import pandas as pd
import sys
import logging
from pathlib import Path
import numpy as np
import pysam

np.set_printoptions(threshold=sys.maxsize)

STEP_NAME = 'UNITE_READS_TO_ICEBERGS'
STEP_OUT_FOLDER = 'UNITE_READS_TO_ICEBERGS_RESULTS'
STEP_LOG_FILE = 'UNITE_READS_TO_ICEBERGS.log'


class ReadsUniter(object):

    def __init__(self,
                 tag: str,
                 tag_rc: str,
                 max_reads_distance: int,
                 in_folder: Path,
                 out_folder: Path) -> None:
        """
        Initiate the ReadsUniter object.

        :param tag: The tag which inject in the cut sites.
        :param tag_rc: The reverse complement of the tag which inject in the cut sites.
        :param max_reads_distance: The maximum distance allowed between read and iceberg (in nts)
         for the read to be joined to the iceberg.
        :param in_folder: The folder where the four files are stored.
        :param out_folder: The folder where the output sam and bam files will be saved.

        :return: A list contains four lines which represent one read.
        """
        self.forward_primer = tag_rc[12:]  # forward primer (22nt)
        self.reverse_primer = tag[6:]  # reverse primer (28nt)
        self.files_dir = in_folder
        self.out_dir = out_folder
        self.max_reads_distance = max_reads_distance

    def get_paired_end_reads_library_primer(self,
                                            read):
        """
        Get the paired end reads library primer from the read header.

        :param read: Aligned read.

        :return: The primer from the read header.
        """
        library_primer = read.query_name.split('_')[4]
        if library_primer == self.forward_primer:
            return "forward-primer"
        elif library_primer == self.reverse_primer:
            return "reverse-primer"
        else:
            return "nomatch"

    def get_paired_end_reads_limit_indexes(self,
                                           read):
        """
        Given aligned read, the function returns the start and end of the pair end reads if their
        template_length is less then 4 * max_reads_distance, else return the start and the end of the current
        read.

        :param read: Aligned read.
        """
        infer_read_length = read.infer_read_length()
        infer_read_length = infer_read_length if infer_read_length is not None else 0

        start_index = read.reference_start
        end_index = read.reference_start + infer_read_length

        if read.flag & 2:  # read mapped in proper pair.
            template_length = read.template_length
            if abs(template_length) > 4 * self.max_reads_distance:
                return start_index, end_index
            if read.flag & 128 and not read.flag & 2048:  # non supplementary read of R2 reads.
                if template_length > 0:  # R2 is before R1 (in the genome).
                    start_index = read.reference_start
                    end_index = start_index + template_length

                elif template_length < 0:  # R1 is before R2 (in the genome).
                    start_index = read.next_reference_start
                    end_index = start_index + abs(template_length)

            elif read.flag & 64 and not read.flag & 2048:  # non supplementary read of R1 reads.
                if template_length > 0:  # R1 is before R2 (in the genome).
                    start_index = read.reference_start
                    end_index = start_index + template_length

                elif template_length < 0:  # R2 is before R1 (in the genome).
                    start_index = read.next_reference_start
                    end_index = start_index + abs(template_length)

        logging.info(f"start_index: {start_index}, end_index: {end_index}")
        return start_index, end_index

    def update_iceberg_orientation(self,
                                   orientations: dict,
                                   read):
        """
        Given the iceberg paired end reads orientations and aligned read, the function update the iceberg orientations
        counter dictionary with the orientation of the read.

        :param orientations: The paired end reads orientations counter dictionary.
        :param read: Aligned read.
        """
        template_length = read.template_length
        if abs(template_length) <= 4 * self.max_reads_distance:
            if read.flag & 128 and not read.flag & 2048:  # non supplementary read of R2 reads.
                if read.template_length > 0:  # R2 is before R1 (in the genome).
                    orientation = 'tag_right_side'
                    if orientation not in orientations:
                        orientations[orientation] = 0
                    orientations[orientation] += 1
                elif read.template_length < 0:  # R1 is before R2 (in the genome).
                    orientation = 'tag_left_side'
                    if orientation not in orientations:
                        orientations[orientation] = 0
                    orientations[orientation] += 1

    def update_iceberg_cut_position_candidates(self,
                                               cut_position_candidates: dict,
                                               read):
        """
        Given the iceberg cut position candidates counter dictionary and aligned read, the function update the
        dictionary at the read inferred cut position if exist, else create new cut position candidate.

        :param cut_position_candidates: The paired end reads orientations counter dictionary.
        :param read: Aligned read.
        """
        if read.flag & 128 and not read.flag & 2048:  # non supplementary read of R2 reads.
            cut_position_candidate = None
            if read.template_length > 0:  # R2 is before R1 (in the genome).
                cut_position_candidate = read.reference_start
            elif read.template_length < 0:  # R1 is before R2 (in the genome).
                cut_position_candidate = read.next_reference_start + abs(read.template_length) - 1

            primer = self.get_paired_end_reads_library_primer(read)
            if cut_position_candidate:
                if cut_position_candidate not in cut_position_candidates:
                    cut_position_candidates[cut_position_candidate] = {}
                if primer not in cut_position_candidates[cut_position_candidate]:
                    cut_position_candidates[cut_position_candidate][primer] = 0
                cut_position_candidates[cut_position_candidate][primer] += 1
            return cut_position_candidate

    def initiate_iceberg(self,
                         read) -> dict:
        """
        Given the first read fields initiate the iceberg boundaries, cut position candidates and
        the iceberg orientation fields.

        :param read: The chromosome of the read.

        :return: The iceberg dictionary.
        """
        cut_position_candidates = {}
        self.update_iceberg_cut_position_candidates(cut_position_candidates, read)

        orientations = {}
        self.update_iceberg_orientation(orientations, read)
        start_index, end_index = self.get_paired_end_reads_limit_indexes(read)

        iceberg = {'chromosome': read.reference_name,
                   'start index': start_index,
                   'end index': end_index,
                   'mapq': [read.mapping_quality],
                   'cut-position-candidates': cut_position_candidates,
                   'read count': 1,
                   'orientations': orientations}
        return iceberg

    def join_read_to_iceberg(self,
                             iceberg: dict,
                             read) -> None:
        """
        Given the read fields, update the iceberg boundaries, cut position candidates and iceberg orientation fields.

        :param iceberg: The iceberg dictionary.
        :param read: The chromosome of the read.
        """
        self.update_iceberg_cut_position_candidates(iceberg['cut-position-candidates'], read)
        self.update_iceberg_orientation(iceberg["orientations"], read)
        paired_end_reads_start, paired_end_reads_end = self.get_paired_end_reads_limit_indexes(read)

        # if cut_position_candidate is not None and iceberg['start index'] > cut_position_candidate:
        #     iceberg['start index'] = cut_position_candidate
        # read_length = read.infer_read_length()
        # end_idx = int(read.reference_start) + (read_length if read_length is not None else 0)
        #
        # if end_idx > iceberg['end index']:
        #     iceberg['end index'] = end_idx
        iceberg['start index'] = min(iceberg['start index'], paired_end_reads_start)
        iceberg['end index'] = max(iceberg['end index'], paired_end_reads_end)
        iceberg['mapq'].append(read.mapping_quality)
        iceberg['read count'] += 1

    def finalize_iceberg(self,
                         iceberg: dict) -> None:
        """
        Given the ready iceberg dictionary the function calculate and update the iceberg values.
        (for calculation that can be calculate only when all iceberg reads were added to the iceberg)

        :param iceberg: The iceberg fields.
        """
        iceberg['mapq'] = round(median(iceberg['mapq']))
        iceberg['iceberg length'] = iceberg['end index'] - iceberg['start index'] + 1
        iceberg_cut_position_candidates = iceberg['cut-position-candidates']
        if iceberg_cut_position_candidates == {}:
            depth = 0
            iceberg['iceberg depth'] = 0
        else:
            highest_candidate = max(iceberg_cut_position_candidates,
                                    key=lambda candidate: sum(iceberg_cut_position_candidates[candidate][primer]
                                                              for primer in iceberg_cut_position_candidates[candidate]))

            depth = sum(iceberg_cut_position_candidates[highest_candidate][primer]
                        for primer in iceberg_cut_position_candidates[highest_candidate])

        iceberg['iceberg depth'] = depth

    def unite_reads_by_distance(self,
                                sorted_bam_file_path: Path,
                                experiment_name: str) -> str:
        """
        Unites the reads into icebergs/clusters (by their chromosome and their index).

        :param sorted_bam_file_path: The sorted bam file name.
        :param experiment_name: The experiment name.

        :return: The iceberg csv file name.
        """
        logging.info(f'unite reads to icebergs based on bps distance (={self.max_reads_distance}) '
                     f'in file: {sorted_bam_file_path.name}')

        pbar = logs.get_progress_bar(step_name='UNITE_READS_TO_ICEBERGS',
                                     experiment_name=experiment_name,
                                     total=sum(1 for _ in generators.get_bam_read(self.files_dir / sorted_bam_file_path.name)),
                                     out_of="read")

        icebergs = []
        # sorted_bam_file_coverage = pysam.AlignmentFile(sorted_bam_file_path, "rb")
        experiment_reads = generators.get_bam_read(self.files_dir / sorted_bam_file_path.name)
        for read in experiment_reads:
            if len(icebergs) == 0:
                icebergs.append(self.initiate_iceberg(read))
                continue
            read_distance_from_curr_iceberg = read.reference_start - icebergs[-1]['end index']
            # if read is from current iceberg - join the read to the current iceberg.
            read_from_same_chromosome_as_iceberg = read.reference_name == icebergs[-1]['chromosome']
            if read_from_same_chromosome_as_iceberg and read_distance_from_curr_iceberg <= self.max_reads_distance:
                self.join_read_to_iceberg(icebergs[-1], read)
            # else - finalize parameters of last iceberg and save the new read as new iceberg.
            else:
                self.finalize_iceberg(icebergs[-1])
                icebergs.append(self.initiate_iceberg(read))
            pbar.update(1)

        if len(icebergs) > 0 and 'iceberg length' not in icebergs[-1]:
            self.finalize_iceberg(icebergs[-1])
        x = pd.DataFrame(icebergs)
        x.to_csv(path_or_buf=self.out_dir / str(sorted_bam_file_path.name).replace('.bam', '.icebergs.csv'),
                 index=False)

        return str(sorted_bam_file_path.name).replace('.bam', '.icebergs.csv')


def validate_arguments(args: dict) -> None:
    """
    Validates the step arguments from the given arguments dict.

    :param args: The iceberg pipeline arguments dict.
    """
    try:
        validation.validate_path_existence(Path(args['curr_in_dir'], args['bam1']))
        validation.validate_path_existence(Path(args['curr_in_dir'], args['bam2']))
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
    Run the UNITE_READS_TO_ICEBERGS step on treatment and control.

    :param args: The pipeline arguments dictionary - dictionary of arguments (as keys) and their values.
    """
    try:
        ru = ReadsUniter(args['EXPERIMENTS']['GENERAL']['EXPERIMENTS_TAG'],
                         args['EXPERIMENTS']['GENERAL']['EXPERIMENTS_TAG_RC'],
                         args['HYPERPARAMATERS']['MAX_READS_DISTANCE'],
                         args['curr_in_dir'],
                         args['curr_out_dir'])

        args['csv1'] = ru.unite_reads_by_distance(args['bam1'], args['EXPERIMENTS']['TX']['NAME'])

        args['csv2'] = ru.unite_reads_by_distance(args['bam2'], args['EXPERIMENTS']['CONTROL']['NAME'])

        args['curr_in_dir'] = validation.set_curr_in_dir(args['OUTPUT_FOLDER_PATH'], STEP_OUT_FOLDER)
        logging.info(f'{STEP_NAME} step done.')
    except Exception as e:
        logging.error(f'Failed to run {STEP_NAME} step, please read the following for more information.',
                      exc_info=False)
        # sys.exit()
        raise e
