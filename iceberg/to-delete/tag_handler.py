import os
import re
import gzip
import subprocess
import argparse
import logging
from pathlib import Path

import utils.validation
from utils import directories_operations
from utils.logs import set_logging_config, log_subprocess
import pandas as pd
#import sys
#import khmer
#import collections

tag_suffix = '.tag'
no_tag_suffix = '.no-tag'
trimmed_suffix = '.trimmed'
untrimmed_suffix = '.untrimmed'
tag_dir = 'TAG-RESULTS'


class tag_handler:
    def __init__(self, tag1: str, tag2: str, k: int, jaccard_threshold: float) -> None:
        """
        __init__:
        initiate the tag_handler object.

        :param tag1: The tag which inject in the cut sites.
        :param tag2: The reverse complement of the tag which inject in the cut sites.
        :param k: The length of the kmers.
        :param jaccard_threshold: The threshold for the jaccard similarity of the read kmers and the tag kmers.

        :return: a list contains four lines which represent one read.
        """
        self.tag1 = tag1
        self.tag2 = tag2
        self.k = k
        self.jaccard_threshold = jaccard_threshold

        #for khmer (not used)
        #alphabet_length = 5  # ['A', 'C', 'G', 'T', 'N']
        #self.kmers_count_graph = khmer.Countgraph(self.k, alphabet_length ** self.k, 1)

        self.tag1_kmer_dict = self.naive_kmer_counts_calculation(tag1)
        self.tag2_kmer_dict = self.naive_kmer_counts_calculation(tag2)


    def get_fastq_lines(self, file: Path) -> list:
        """
        get_fastq_lines:
        Yield one read from the given fastq file (for iterate the reads in the file).

        :param file: The fastq file path.

        :return: A list contains four lines which represent one read.
        """
        if re.search('.gz$', str(file)):
            fastq = gzip.open(str(file), 'rb')
        else:
            fastq = open(str(file), 'r')
        with fastq as f:
            while True:
                l1 = f.readline()
                if not l1:
                    break
                l2 = f.readline()
                l3 = f.readline()
                l4 = f.readline()
                yield [l1, l2, l3, l4]

    def naive_kmer_counts_calculation(self, sequence: str) -> dict:
        """
        naiveKmerCountsCalculation:
        Creates a kmers dictionary for the sequence with kmers as keys and the
        amount of their appearance in the sequence as their values).

        :param sequence: The read sequence.

        :return: The kmers dictionary.
        """
        k_mer_counts_dict = {}
        # Loop over the string (read/tag)
        for i in range(len(sequence) - self.k):
            kmer = sequence[i:i + self.k]
            # Add the kmer to the dictionary if it's not there or sum +1 to its value
            if kmer not in k_mer_counts_dict:
                k_mer_counts_dict[kmer] = 0
            k_mer_counts_dict[kmer] += 1
        # sort the kmer dictionary by amount (not used yet)
        #k_mer_counts_dict = collections.OrderedDict(sorted(k_mer_counts_dict.items()))
        return k_mer_counts_dict

    # kmer with khmer (shay, time not improved much)
    # def naive_kmer_counts_calculation(self, sequence):
    #     """
    #     naiveKmerCountsCalculation:
    #     Creates a kmers dictionary for the sequence
    #     (with kmers as keys and amount of appearance in the sequence as the values).
    #
    #     return: The kmers dictionary.
    #     """
    #     kmers = self.kmers_count_graph.get_kmers(sequence)
    #
    #     k_mer_counts_dict = {}
    #     # Loop over the string (read/tag)
    #     for kmer in kmers:
    #         # Add the kmer to the dictionary if it's not there or sum +1 to its value
    #         if kmer not in k_mer_counts_dict:
    #             k_mer_counts_dict[kmer] = 0
    #         k_mer_counts_dict[kmer] += 1
    #     # sort the kmer dictionary by amount
    #     #k_mer_counts_dict = collections.OrderedDict(sorted(k_mer_counts_dict.items()))
    #     return k_mer_counts_dict

    def jaccard_similarity(self, kmer_dict: dict) -> (float, str):
        """
        jaccard_similarity:
        Calculate the Jaccard similarity of the read kmers with tag1 and tag2 kmers,
        then chooses the best tag form for the read.

        :param kmer_dict: The read dict where each key is a kmer and it value is the kmer count in the read.

        :return: A float which represent the highest Jaccard similarity result, and the tag ID (tag1 or tag2).
        """
        a = set(kmer_dict.keys())
        t1 = set(self.tag1_kmer_dict.keys())
        t2 = set(self.tag2_kmer_dict.keys())
        jaccard_t1 = len(a.intersection(t1)) / len(a.union(t1))
        jaccard_t2 = len(a.intersection(t2)) / len(a.union(t2))
        if jaccard_t1 >= jaccard_t2:
            tag_id = 'tag1'
            jaccard = jaccard_t1
        else:
            tag_id = 'tag2'
            jaccard = jaccard_t2
        return jaccard, tag_id

    def add_paired_end_reads_to_kmers_df(self, read_r1: list, read_r2: list, r1_res: float, r2_res: float, r1_tag_id: str, r2_tag_id: str, kmers_memory: list) -> None:
        """
        add_paired_end_reads_to_kmers_df:
        given the kmer memory list and the next paired-end reads data, adds the new paired-end
        reads data to the memory.

        :param read_r1: The read from r1 that need to be add to memory (paired with read_r2).
        :param read_r2: The read from r2 that need to be add to memory (paired with read_r1).
        :param r1_res: The Jaccard similarity result of read_r1 and r1_tag_id
        :param r2_res: The Jaccard similarity result of read_r2 and r2_tag_id
        :param r1_tag_id: The tag that gives the highest Jaccard similarity result with read_r1.
        :param r2_tag_id: The tag that gives the highest Jaccard similarity result with read_r2.
        :param kmers_memory: The list with the previous paired-end reads results.
        """
        kmers_memory.append({'read r1': read_r1[1], 'read r1 id': read_r1[0], 'jaccard results r1': r1_res, 'r1 tag id': r1_tag_id,
                             'read r2': read_r2[1], 'read r2 id': read_r2[0], 'jaccard results r2': r2_res, 'r2 tag id': r2_tag_id})

    def classify_tag_pair_end_reads(self, r1_file_name: str, r2_file_name: str, files_dir: Path, out_dir: Path) -> tuple:
        """
        classify_tag_pair_end_reads:
        split the reads in to two files - reads with tag and reads without tag based on the Jaccard similarity on
        the kmers of the read and the tags.

        :param r1_file_name: The fastq file name of r1.
        :param r2_file_name: The fastq file name of r2.
        :param files_dir: The directory where the four files are stored.
        :param out_dir: The directory where the output fastq files will be saved.

        :return: Four fastq files names - two files for r1 (tag and without tag) and two files for r2 (tag and without tag).
        """
        logging.info(f'classify reads in {r1_file_name} and {r2_file_name} to- reads with tag and reads without tag, '
                     f'using Jaccard similarity on each read with the tag kmers (k={self.k})...')
        r1_tag = open(out_dir / r1_file_name.replace('.fastq', tag_suffix + '.fastq'), 'w')
        r1_no_tag = open(out_dir / r1_file_name.replace('.fastq', no_tag_suffix + '.fastq'), 'w')
        r2_tag = open(out_dir / r2_file_name.replace('.fastq', tag_suffix + '.fastq'), 'w')
        r2_no_tag = open(out_dir / r2_file_name.replace('.fastq', no_tag_suffix + '.fastq'), 'w')

        kmers_memory = []
        for read_r1, read_r2 in zip(self.get_fastq_lines(files_dir / r1_file_name), self.get_fastq_lines(files_dir / r2_file_name)):
            # generating the kmer dict for the read
            r1_curr_read_kmer_dict = self.naive_kmer_counts_calculation(read_r1[1])
            r1_res, r1_tag_id = self.jaccard_similarity(r1_curr_read_kmer_dict)

            r2_curr_read_kmer_dict = self.naive_kmer_counts_calculation(read_r2[1])
            r2_res, r2_tag_id = self.jaccard_similarity(r2_curr_read_kmer_dict)

            self.add_paired_end_reads_to_kmers_df(read_r1, read_r2, r1_res, r2_res, r1_tag_id, r2_tag_id, kmers_memory)
            if r1_res > self.jaccard_threshold or r2_res > self.jaccard_threshold:
                for line_r1, line_r2 in zip(read_r1, read_r2):
                    r1_tag.write(line_r1)
                    r2_tag.write(line_r2)
            else:
                for line_r1, line_r2 in zip(read_r1, read_r2):
                    r1_no_tag.write(line_r1)
                    r2_no_tag.write(line_r2)
        r1_tag.close()
        r1_no_tag.close()
        r2_tag.close()
        r2_no_tag.close()
        kmers_memory = pd.DataFrame(kmers_memory)
        kmers_memory.to_csv(out_dir / r1_file_name.replace('.fastq', '-kmer-validation.csv'))
        return (r1_file_name.replace('.fastq', tag_suffix + '.fastq'),
                r1_file_name.replace('.fastq', no_tag_suffix + '.fastq'),
                r2_file_name.replace('.fastq', tag_suffix + '.fastq'),
                r2_file_name.replace('.fastq', no_tag_suffix + '.fastq'))

    #omri tool command (not so well)
    # def trim_tag_pair_end_reads(self, r1, r2, files_dir, out_dir):
    #         '''
    #         trim_tag_pair_end_reads:
    #         Trim the tags from the reads (in a pair-end reads experiment).
    #         '''
    #         try:
    #             logging.info('trim tag only in reads that classified as reads with tag...')
    #             out_r = r1.replace('.fastq', trimmed_suffix + '+-' + '.fastq')
    #             out_l = r2.replace('.fastq', trimmed_suffix + '+-' + '.fastq')
    #             r1_out_untrimmed = r1.replace('.fastq', untrimmed_suffix + '.fastq')
    #             r2_out_untrimmed = r2.replace('.fastq', untrimmed_suffix + '.fastq')
    #             logging.info('cutadapt summary:')
    #             # here cutadapt split the outputs to files.
    #             # cmd = f'~/.local/bin/cutadapt -b {self.tag1} -B {self.tag2} --untrimmed-output {str(out_dir / r1_out_untrimmed)} --untrimmed-paired-output {str(out_dir / r2_out_untrimmed)} ' \
    #             #       f'-o {str(out_dir / out_r)} -p {str(out_dir / out_l)} ' \
    #             #       f'{str(files_dir / r1)} {str(files_dir / r2)}'
    #
    #             cmd = f'~/.local/bin/cutadapt -b "{self.tag1};min_overlap={self.k}" -B "{self.tag2};min_overlap=8" ' \
    #                   f'-o {str(out_dir / out_r)} -p {str(out_dir / out_l)} ' \
    #                   f'{str(files_dir / r1)} {str(files_dir / r2)}'
    #
    #             cmd = f'~/.local/bin/cutadapt -a "{self.tag1};min_overlap={self.k}" -a "{self.tag2};min_overlap={self.k}" ' \
    #                   f'-a "{self.tag1rev};min_overlap={self.k}" -a "{self.tag2rev};min_overlap={self.k}" -n 4 ' \
    #                   f'-o {str(out_dir / out_r)} -p {str(out_dir / out_l)} ' \
    #                   f'{str(files_dir / r1)} {str(files_dir / r2)}'
    #
    #             proc = subprocess.Popen(cmd, shell=True, env=os.environ.copy(), stdout=subprocess.PIPE,
    #                                     universal_newlines=True, stderr=subprocess.STDOUT)
    #             self.log_subprocess(proc)
    #
    #
    #             return out_r, out_l
    #         except Exception as e:
    #             logging.error(
    #                 f'couldnt run cutadapt step (read trimming) on files: {r1} and {r2} please check files paths and make sure cutadapt is install correctly (by typing  ~/.local/bin/cutadapt --version in conlsole).',
    #                 exc_info=False)
    #             # sys.exit()
    #             raise e

    def trim_tag_pair_end_reads(self, r1_file_name: str, r2_file_name: str, files_dir: Path, out_dir: Path) -> (str, str):
        """
        trim_tag_pair_end_reads:
        Trims the tags from the given experiment paired-end reads.

        :param r1_file_name: The fastq file name of r1.
        :param r2_file_name: The fastq file name of r2.
        :param files_dir: The directory where the four files are stored.
        :param out_dir: The directory where the output fastq files will be saved.

        :return: The names of the fastq files of the trimmed reads.
        """
        try:
            logging.info('trim tag only in reads that classified as reads with tag...')
            out_r = r1_file_name.replace('.fastq', trimmed_suffix + 'R1:tag1-R2:tag2' + '.fastq')
            out_l = r2_file_name.replace('.fastq', trimmed_suffix + 'R1:tag1-R2:tag2' + '.fastq')
            logging.info('cutadapt summary:')

            #here cutadapt split the outputs to files if classify_tag_pair_end_reads not used (kmer guy & shay).
            #r1_out_untrimmed = r1.replace('.fastq', untrimmed_suffix + '.fastq')
            #r2_out_untrimmed = r2.replace('.fastq', untrimmed_suffix + '.fastq')
            # cmd = f'~/.local/bin/cutadapt -b {self.tag1} -B {self.tag2} --untrimmed-output {str(out_dir / r1_out_untrimmed)} --untrimmed-paired-output {str(out_dir / r2_out_untrimmed)} ' \
            #       f'-o {str(out_dir / out_r)} -p {str(out_dir / out_l)} ' \
            #       f'{str(files_dir / r1)} {str(files_dir / r2)}'


            # tag1-r1,tag2-r2 paired-end reads.
            cmd = f'~/.local/bin/cutadapt -b "{self.tag1};min_overlap={self.k}" -B "{self.tag2[12:]};min_overlap=15" ' \
                  f'-o {str(out_dir / out_r)} -p {str(out_dir / out_l)} ' \
                  f'{str(files_dir / r1_file_name)} {str(files_dir / r2_file_name)}'
            proc = subprocess.Popen(cmd, shell=True, env=os.environ.copy(), stdout=subprocess.PIPE,
                                    universal_newlines=True, stderr=subprocess.STDOUT)
            log_subprocess(proc)

            # tag1-r2,tag2-r1 paired-end reads.
            out_r_final = r1_file_name.replace('.fastq', trimmed_suffix + 'R1:tag2-R2:tag1' + '.fastq')
            out_l_final = r2_file_name.replace('.fastq', trimmed_suffix + 'R1:tag2-R2:tag1' + '.fastq')

            cmd = f'~/.local/bin/cutadapt -b "{self.tag2};min_overlap={self.k}" -B "{self.tag1};min_overlap={self.k}" ' \
                  f'-o {str(out_dir / out_r_final)} -p {str(out_dir / out_l_final)} ' \
                  f'{str(out_dir / out_r)} {str(out_dir / out_l)}'
            proc = subprocess.Popen(cmd, shell=True, env=os.environ.copy(), stdout=subprocess.PIPE,
                                    universal_newlines=True, stderr=subprocess.STDOUT)
            log_subprocess(proc)
            return out_r_final, out_l_final

        except Exception as e:
            logging.error(f'couldnt run cutadapt step (read trimming) on files: {r1_file_name} and {r2_file_name}'
                          f' please check files paths and make sure cutadapt is install correctly '
                          f'(by typing  ~/.local/bin/cutadapt --version in conlsole).', exc_info=False)
            # sys.exit()
            raise e


def validate_arguments(args: dict) -> None:
    """
    validate_arguments:
    Validates the step arguments from the given arguments dict.

    :param args: The arguments dict (see main function description).
    """
    try:
        utils.validation.validate_path_existence(args['curr_in_dir'])
        if args['tag1'] == '' or args['tag2'] == '':
            raise ValueError()
        if args['read1_1'] != '':
            for file in [args['read1_1'], args['read2_1']]:
                if file == '':
                    raise ValueError()
                utils.validation.validate_path_existence(Path(args['curr_in_dir'], file))
        if args['read1_2'] != '':
            for file in [args['read1_2'], args['read2_2']]:
                if file == '':
                    raise ValueError()
                utils.validation.validate_path_existence(Path(args['curr_in_dir'], file))
        if args['read1_1'] == '' or args['read1_1'] == '' and args['read1_2'] == '':
            raise ValueError()

    except ValueError as e:
        logging.error('not all necessary arguments for tag step where giving.')
        raise e


def prepare_step(args: dict) -> None:
    """
    prepare_step:
    Prepares the step by - setting the logging configurations, setting the output directory and
    validate the step arguments.

    :param args: A dict of arguments (as keys) and their values (see main function description)
    """
    set_logging_config('tag', Path(args['out_dir'], 'DEBUG', 'LOGS', 'tag.log'))
    args['curr_out_dir'] = utils.validation.set_curr_out_dir(args['out_dir'], tag_dir)
    validate_arguments(args)


def main(args: dict) -> None:
    """
    main:
    finds all the reads that contains the tag in the treatment fastq files using
    kmers approach, then trims the tag from those reads.
    If mock files are provided, apply the same for mock.

    :param args: A dict of arguments (as keys) and their values.
     Can be one of the following:
     1.The pipeline arguments dict - if you run analyzer_old.py.
     2.The step arguments dict - if you run this class directly.
    """
    try:
        th = tag_handler(args['tag1'], args['tag2'], args['kmer'], args['jaccard_threshold'])
        #th = tag_handler('TTGAGTTGTCATATGTTAATAACGGTAT', 'ACATATGACAACTCAATTAAAC', args['kmer'], args['jaccard_threshold'])

        r1_tag_1, r1_no_tag_1, r2_tag_1, r2_no_tag_1 = th.classify_tag_pair_end_reads(args['read1_1'], args['read2_1'],
                                                                                      # args['index1_1'], args['index2_1'],
                                                                                      args['curr_in_dir'], args['curr_out_dir'])

        r1_tag_trimmed_1, r2_tag_trimmed_1 = th.trim_tag_pair_end_reads(r1_tag_1, r2_tag_1,
                                                                        #args['index1_1'], args['index2_1'],
                                                                        args['curr_out_dir'], args['curr_out_dir'])
        # if classify_tag_pair_end_reads (kmer guy & shay) is not used
        # r1_tag_trimmed_1, r2_tag_trimmed_1 = th.trim_tag_pair_end_reads(args['read1_1'], args['read2_1'],
        #                                                                 args['index1_1'], args['index2_1'],
        #                                                                 args['curr_in_dir'], args['curr_out_dir'])
        args['read1_1'] = r1_tag_trimmed_1
        args['read2_1'] = r2_tag_trimmed_1

        if args['read1_2'] != '':
            r1_tag_2, r1_no_tag_2, r2_tag_2, r2_no_tag_2 = th.classify_tag_pair_end_reads(args['read1_2'], args['read2_2'],
                                                                                          #args['index1_2'], args['index2_2'],
                                                                                          args['curr_in_dir'], args['curr_out_dir'])

            r1_tag_trimmed_2, r2_tag_trimmed_2 = th.trim_tag_pair_end_reads(r1_tag_2, r2_tag_2,
                                                                            #args['index1_2'], args['index2_2'],
                                                                            args['curr_out_dir'], args['curr_out_dir'])

            # r1_tag_trimmed_2, r2_tag_trimmed_2 = th.trim_tag_pair_end_reads(args['read1_2'], args['read2_2'],
            #                                                                 args['index1_2'], args['index2_2'],
            #                                                                 args['curr_in_dir'], args['curr_out_dir'])
            args['read1_2'] = r1_tag_trimmed_2
            args['read2_2'] = r2_tag_trimmed_2
        args['curr_in_dir'] = utils.validation.set_curr_in_dir(args['out_dir'], tag_dir)
        logging.info('tag step done.')
    except Exception as e:
        logging.error('error occurred while running tag step, please read the following for more information.',
                      exc_info=False)
        # sys.exit()
        raise e


def parse_args() -> dict:
    """
    parse_args:
    Parses step arguments (used only in the --help command or if user run this class alone).

    :return: A dict of the step arguments (as keys) and their values.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--files_dir', help="ff")
    parser.add_argument('--read1_1', help="ff")
    parser.add_argument('--read2_1', help="ff")
    #parser.add_argument('index1_1', help="ff")
    #parser.add_argument('index2_1', help="ff")
    parser.add_argument('--out_dir', help="ff")
    parser.add_argument('--tag1', default='', help="ff")
    parser.add_argument('--tag2', default='', help="ff")
    parser.add_argument('--read1_2', default='', help="ff", metavar='')
    parser.add_argument('--read2_2', default='', help="ff", metavar='')
    #parser.add_argument('--index1_2', default='', help="ff", metavar='')
    #parser.add_argument('--index2_2', default='', help="ff", metavar='')
    parser.add_argument('--kmer', default=8, help="ff", type=int, metavar='')
    parser.add_argument('--jaccard_threshold', default=0.05, type=float, help="ff", metavar='')

    args, rest = parser.parse_known_args()
    return dict([(k, v) for k, v in vars(args).items()])


if __name__ == '__main__':
    args = parse_args()
    args['curr_in_dir'] = Path(args['files_dir'])
    utils.validation.set_curr_out_dir(Path(args['out_dir'], 'DEBUG'), 'LOGS')
    prepare_step(args)
    # print(args['curr_in_dir'])
    main(args)
