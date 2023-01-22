import os
import warnings
import argparse
import logging
from pathlib import Path
import cProfile, pstats
import subprocess
import sys

import utils.validation
from utils.logs import set_logging_config
from iceberg import settings
from utils import directories_operations, reads_counter

warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

def print_step_help_massage(steps: list) -> None:
    """
    print_step_help_massage:
    prints the help massage of each step in the given list.

    :param steps: list of steps names.
    """
    for step in steps:
        if step in [step.STEP_NAME for step in settings.STEPS_METADATA]:
            cmd = '{} -h'.format(settings.STEPS_METADATA[step]['cmd'])
            subprocess.check_call(cmd, shell=True, env=os.environ.copy())


def parse_steps_arguments() -> tuple:
    """
    parse_steps_arguments:
    Parses all steps arguments (step_name as key, the value is true if user specify --step_name, else the value is false).

    :return: A dict of the steps and their values.
    """
    global parser
    parser = argparse.ArgumentParser(description="preform analysis by steps for treatment and mock.\n"
                                                 "if only treatment or mock was given, preform single expirement analysis (analysis without expirements comparison).\n"
                                                 "if no one of the following arguments was given, preform all iceberg steps, else only preform the specified steps.\n ",
                                     add_help=True)
    parser.add_argument('-u', '--umi', action="store_true", help="ff")
    parser.add_argument('-l', '--libraries_detection', action="store_true", help="ff")
    parser.add_argument('-t', '--tag', action="store_true", help="ff")
    parser.add_argument('-b', '--bwa', action="store_true", help="ff")
    parser.add_argument('-n', '--unite', action="store_true", help="ff")
    parser.add_argument('-c', '--classify', action="store_true", help="ff")
    parser.add_argument('-g', '--guide_alignment', action="store_true", help="ff")
    parser.add_argument('-r', '--report', action="store_true", help="ff")
    return parser.parse_known_args()


def parse_steps_internal_arguments(steps_internal_args: dict) -> tuple:
    """
    parse_steps_internal_arguments:
    Parses all steps internal arguments.

    :return: A dict of the steps internal arguments (as keys) and their values.
    """
    global parser
    subparsers = parser.add_subparsers()
    rest_parser = subparsers.add_parser('rest', help='a help')
    rest_parser.add_argument('--files_dir', default='', help="ff")
    rest_parser.add_argument('--read1_1', default='', help="ff")
    rest_parser.add_argument('--read2_1', default='', help="ff")
    rest_parser.add_argument('--index1_1', default='', help="ff")
    rest_parser.add_argument('--index2_1', default='', help="ff")
    rest_parser.add_argument('--out_dir', help="ff")
    rest_parser.add_argument('--read1_2', default='', help="ff", metavar='')
    rest_parser.add_argument('--read2_2', default='', help="ff", metavar='')
    rest_parser.add_argument('--index1_2', default='', help="ff", metavar='')
    rest_parser.add_argument('--index2_2', default='', help="ff", metavar='')
    rest_parser.add_argument('-q', '--min_qual', default=15, type=int, help="ff",  metavar='')
    rest_parser.add_argument('-f', '--min_freq', default=0.9, type=float, help="ff", metavar='')
    rest_parser.add_argument('--tag1', default='GTTTAATTGAGTTGTCATATGTTAATAACGGTAT', help="ff")
    rest_parser.add_argument('--tag2', default='ATACCGTTATTAACATATGACAACTCAATTAAAC', help="ff")
    rest_parser.add_argument('--jaccard_threshold', default=0.05, type=float, help="ff", metavar='')
    rest_parser.add_argument('--kmer', default=8, help="ff", type=int, metavar='')

    rest_parser.add_argument('--genome',default='', help="ff")
    rest_parser.add_argument('--sam1', default='', help="ff", metavar='')
    rest_parser.add_argument('--sam2', default='', help="ff", metavar='')
    rest_parser.add_argument('--reads_distance', default=2000, type=int, help="ff")
    rest_parser.add_argument('--icebergs_distance', default=4000, type=int, help="ff")
    rest_parser.add_argument('--csv1', default='', help="ff")
    rest_parser.add_argument('--csv2', default='', help="ff")
    rest_parser.add_argument('--noise_threshold', default=0.999, type=float, help="ff", metavar='')
    rest_parser.add_argument('--crispr_activity_threshold', default=0.75, type=float, help="ff", metavar='')

    rest_parser.add_argument('--guide_rna', help="ff")

    parser.add_argument('--guide_alignment_txt', help="ff")
    rest_parser.add_argument('--bam1', help="ff")
    rest_parser.add_argument('--csv_ca', help="ff")
    rest_parser.add_argument('--csv_sb', help="ff")
    rest_parser.add_argument('--csv_noise', help="ff")
    rest_parser.add_argument('--bam2', help="ff")
    return rest_parser.parse_known_args(steps_internal_args)

# def check_if_dir_exist(files_dir_path):
#     try:
#         if not os.path.exists(files_dir_path):
#             raise FilesDirectoryNotExistError
#     except FilesDirectoryNotExistError:
#         logging.error('files directory: {} not exist'.format(files_dir_path), exc_info=True)
#         sys.exit()
#
#
# def check_dir_files_amount(files_dir_path):
#     try:
#         if len(fnmatch.filter(os.listdir(files_dir_path), '*.fastq')) < 4 and len(fnmatch.filter(os.listdir(files_dir_path), '*.fasta')) < 4:
#             raise FilesNumberToSmallError
#         if len(fnmatch.filter(os.listdir(files_dir_path), '*.fastq')) == 4 or len(fnmatch.filter(os.listdir(files_dir_path), '*.fasta')) == 4:
#             logging.warning('files directory: {} contain enough files only for single paired end sample analysis.'.format(files_dir_path))
#     except FilesNumberToSmallError:
#         logging.error('files directory: {} not contain enough files for single or dual paired end samples analysis.'.format(files_dir_path), exc_info=True)
#         sys.exit()

# def validate_files_dir(files_dir_path):
#         check_if_dir_exist(files_dir_path)
#         check_dir_files_amount(files_dir_path)

# def prepare_out_dir(out_dir_path):
#     if not os.path.exists(out_dir_path):
#         os.makedirs(out_dir_path)
#         logging.info('out directory: {} created in order to store the iceberg results.'.format(out_dir_path))
#     else:
#         logging.warning('out directory: {} already exists, if you repeats on analysis which you all ready done the files will be overwritten.'.format(out_dir_path))
#     profiler_dir = Path(out_dir_path, 'DEBUG', 'LOGS')
#     if not os.path.exists(profiler_dir):
#         os.makedirs(profiler_dir)


#def validate_ram_capacity(): pass


def enable_profiler() -> cProfile.Profile:
    """
    enable_profiler:
    enable profiler (using cProfile).

    :return: the profiler generated using cProfile
    """
    profiler = cProfile.Profile()
    profiler.enable()
    return profiler


def disable_profiler(profiler: cProfile.Profile, save_profiler_to_file: bool) -> None:
    """
    disable_profiler:
    disable profiler and export sneakviz file for visualization if needed.

    :param profiler: The profiler (generated by cProfile).
    :param save_profiler_to_file: If sneakviz file is wanted set to true  else set to false.
    """
    profiler.disable()
    if save_profiler_to_file:
        stats = pstats.Stats(profiler).sort_stats('cumtime')
        stats.strip_dirs()
        stats.dump_stats(Path(steps_internal_args['out_dir'], 'DEBUG', 'profiler'))


def main(steps: dict, args: dict) -> None:
    """
    main:
    Running the iceberg pipeline by steps.

    :param steps: A dict of the steps as keys and boolean as value (true if user enter --step_name, else false)
    :param args: A dict of all steps internal arguments (as keys) and their values.
    """
    args['curr_in_dir'] = args['files_dir']
    args['remaining_reads'] = {}

    treatment_in_file_name = args['read1_1']
    mock_in_file_name = args['read1_2']
    start = True
    for step in ['umi', 'libraries_detection', 'tag', 'bwa', 'unite', 'classify', 'guide_alignment', 'report']:
        if steps[step]:
            if start:
                reads_counter.count_experiments_reads_by_step(args, step, start=True)
                start = False
            logging.debug(f'running {step} step....')
            settings.STEPS_METADATA[step]['class'].prepare_step(args)
            settings.STEPS_METADATA[step]['class'].main(args)
            set_logging_config('analyzer', Path(args['out_dir'], 'DEBUG', 'LOGS', 'analyzer.log'), writing_mode='a')
            reads_counter.count_experiments_reads_by_step(args, step)

    reads_counter.generate_reads_count_funnel_chart_html(args['remaining_reads'], args['curr_out_dir'],
                                                         treatment_in_file_name, mock_in_file_name)
    reads_counter.generate_command_line_html(args['curr_out_dir'], steps, args)

    logging.info('analyzer done.')


if __name__ == '__main__':
    print(sys.argv)
    try:
        profiler = enable_profiler()
        save_profiler_to_file = False
        steps_args, steps_internal_args = parse_steps_arguments()
        steps_args = dict([(k, v) for k, v in vars(steps_args).items()])
        if len(steps_internal_args) == 0:
            print_step_help_massage([k for k in steps_args.keys() if steps_args[k]])
        else:
            save_profiler_to_file = True
            steps_internal_args, _ = parse_steps_internal_arguments(steps_internal_args)
            steps_internal_args = dict([(k, v) for k, v in vars(steps_internal_args).items()])
            if steps_internal_args['files_dir'] != '' and steps_internal_args['out_dir'] != '':
                set_logging_config('analyzer', None)
                steps_internal_args['files_dir'] = utils.validation.set_curr_in_dir(Path(steps_internal_args['files_dir']), '')
                utils.validation.set_curr_out_dir(Path(steps_internal_args['out_dir'], 'DEBUG'), 'LOGS')
                set_logging_config('analyzer', Path(steps_internal_args['out_dir'], 'DEBUG', 'LOGS', 'analyzer.log'), writing_mode='w')

                # validate_files_dir(Path(steps_internal_args['files_dir']))
                # prepare_out_dir(Path(steps_internal_args['out_dir']))
                # validate_ram_capacity()
                if all(not value for key, value in steps_args.items()):
                    main(settings.ICEBERG_ALL_STEPS_DICT, steps_internal_args)
                elif len(steps_internal_args) != 0:
                    main(steps_args, steps_internal_args)
    except SystemExit:  # for help message.
        pass

    except:
        logging.error('error occurred while running the iceberg pipeline, please read the following for more information.',
                      exc_info=True)
        print('sssss')
        sys.exit()
    finally:
        print('ffff')
        disable_profiler(profiler, save_profiler_to_file)
        logging.shutdown()



