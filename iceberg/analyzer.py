"""
This Python module is used to run the iceberg pipeline by steps for GuideSeq treatment and conrol experiments.
"""

import settings
from utils import arguments_parsing, reads_counter, logs, validation

from collections import OrderedDict
import warnings
import logging
from pathlib import Path
import cProfile
import pstats

warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

global profiler
global args


def enable_profiler() -> None:
    """
    Enable profiler (using cProfile).

    :return: the profiler generated using cProfile.
    """
    global profiler
    profiler = cProfile.Profile()
    profiler.enable()


def disable_profiler(save_to_file=True) -> None:
    """
    Disable profiler and generate sneakviz file for visualization if needed.

    :param save_to_file: If sneakviz file is wanted set to true, else, set to false.
    """
    global profiler
    profiler.disable()
    if save_to_file:
        stats = pstats.Stats(profiler).sort_stats('cumtime')
        stats.strip_dirs()
        stats.dump_stats(Path(args['OUTPUT_FOLDER_PATH'], 'DEBUG', 'profiler'))


def run_step(args,
             step_module) -> None:
    """
    Run single step of the iceberg pipeline.

    :param args: A dict of all steps internal arguments (as keys) and their values.
    :param step_module: the step module.

    """
    logging.debug(f'Running {step_module.STEP_NAME} step....')
    step_module.prepare_step(args)
    step_module.main(args)
    logs.set_logging_config('analyzer',
                            args['OUTPUT_FOLDER_PATH'] / settings.LOGS_FOLDER_RELATIVE_PATH / 'analyzer.log',
                            writing_mode='a')


def main(args) -> None:
    """
    Running the iceberg pipeline by steps.

    :param args: A dict of all steps internal arguments (as keys) and their values.
    """
    logs.set_logging_config('ICEBERG_ANALYZER', None)
    _ = validation.set_curr_out_dir(args['OUTPUT_FOLDER_PATH'], settings.LOGS_FOLDER_RELATIVE_PATH)
    logs.set_logging_config('ICEBERG_ANALYZER',
                            args['OUTPUT_FOLDER_PATH'] / settings.LOGS_FOLDER_RELATIVE_PATH / 'analyzer.log',
                            writing_mode='w')

    logging.info(settings.START_MESSAGE)
    args['curr_in_dir'] = ""
    args['remaining_reads'] = OrderedDict()
    pbar = logs.get_progress_bar(step_name='ICEBERG_ANALYZER',
                                 experiment_name="",
                                 total=len(args["ANALYZER STEPS"]),
                                 out_of="step",
                                 bar_length=40)

    start = True
    for step in settings.STEPS_METADATA:
        if args["ANALYZER STEPS"][step.STEP_NAME]:
            if start:
                reads_counter.count_experiments_reads_by_step(args, step.STEP_NAME, start=True)
                start = False
            run_step(args, step)
            reads_counter.count_experiments_reads_by_step(args, step.STEP_NAME)
            pbar.update(1)
    logging.info('analyzer done.')


if __name__ == '__main__':
    try:
        print(settings.START_MESSAGE)
        args = arguments_parsing.parse()
        if args['enable_profiler']:
            enable_profiler()

        main(args)

        if args['enable_profiler']:
            disable_profiler()

    except Exception as e:
        logging.error('Failed to run the iceberg pipeline, please read the following for more information.',
                      exc_info=True)
        if args is not None and args['enable_profiler']:
            disable_profiler(save_to_file=False)
    finally:
        logging.shutdown()



