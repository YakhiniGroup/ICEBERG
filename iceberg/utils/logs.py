"""
This Python module is used to handle the logging configurations and the pipeline progress bars.
"""

import logging
from tqdm import tqdm
from pathlib import Path
import subprocess


class TqdmLoggingHandler(logging.Handler):
    def __init__(self, level=logging.NOTSET):
        super().__init__(level)

    def emit(self, record):
        try:
            msg = self.format(record)
            tqdm.write(msg)
            self.flush()
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            self.handleError(record)


logging.addLevelName(logging.DEBUG, 'debug')
logging.addLevelName(logging.INFO, 'info')
logging.addLevelName(logging.WARNING, 'warning')
logging.addLevelName(logging.ERROR, 'error')
logging.addLevelName(logging.CRITICAL, 'critical')

import sys
# def set_logging_config(step_name: str,
#                        log_file_path: Path,
#                        writing_mode: str = 'w',
#                        ) -> None:
#     """
#     set_logging_config:
#     Sets the logging configurations for step in the iceberg pipeline.
#
#     :param step_name: The step name.
#     :param log_file_path: The path to the logging file.
#     :param writing_mode: The log_file_path writing mode.
#     """
#     try:
#         if step_name in loggers:
#             return loggers[step_name]
#
#         logger = logging.getLogger(step_name)
#         logger.setLevel(logging.DEBUG)
#
#         if log_file_path is None:
#             formatter = logging.Formatter(fmt='%(asctime)s - analyzer - %(levelname)s: %(message)s',
#                                           datefmt='%m-%d-%Y %H:%M:%S')
#
#             stdout_handler = logging.StreamHandler(sys.stdout)
#             stdout_handler.setLevel(logging.DEBUG)
#             stdout_handler.setFormatter(formatter)
#
#             #logger.addHandler(stdout_handler)
#             logger.addHandler(TqdmLoggingHandler())
#         elif step_name == 'analyzer':
#             formatter = logging.Formatter(fmt='%(asctime)s - analyzer - %(levelname)s: %(message)s',
#                                           datefmt='%m-%d-%Y %H:%M:%S')
#
#             stdout_handler = logging.StreamHandler(sys.stdout)
#             stdout_handler.setLevel(logging.DEBUG)
#             stdout_handler.setFormatter(formatter)
#
#             file_handler = logging.FileHandler(log_file_path, mode=writing_mode)
#             file_handler.setLevel(logging.DEBUG)
#             file_handler.setFormatter(formatter)
#
#             #logger.addHandler(stdout_handler)
#             logger.addHandler(file_handler)
#             logger.addHandler(TqdmLoggingHandler())
#
#         else:
#             formatter = logging.Formatter(fmt=f'%(asctime)s - analyzer - {step_name} - %(levelname)s: %(message)s',
#                                           datefmt='%m-%d-%Y %H:%M:%S')
#
#             stdout_handler = logging.StreamHandler(sys.stdout)
#             stdout_handler.setLevel(logging.DEBUG)
#             stdout_handler.setFormatter(formatter)
#
#             file_handler = logging.FileHandler(log_file_path, mode=writing_mode)
#             file_handler.setLevel(logging.DEBUG)
#             file_handler.setFormatter(formatter)
#
#             #logger.addHandler(stdout_handler)
#             logger.addHandler(file_handler)
#             logger.addHandler(TqdmLoggingHandler())
#         loggers[step_name] = logger
#         return logger
#
#                 logging.basicConfig(format=f'%(asctime)s - analyzer - %(levelname)s: %(message)s',
#                                     datefmt='%Y-%m-%d %H:%M:%S')
#
#     except Exception as e:
#         logging.error('couldnt config logging basic configuration.')
#         raise e


def set_logging_config(step_name: str,
                       log_file_path: Path,
                       writing_mode: str = 'w') -> None:
    """
    set_logging_config:
    Sets the logging configurations for step in the iceberg pipeline.

    :param step_name: The step name.
    :param log_file_path: The path to the logging file.
    :param writing_mode: The log_file_path writing mode.
    """
    try:
        if log_file_path is None:
            logging.basicConfig(format=f'%(asctime)s - ICEBERG_ANALYZER - %(levelname)s: %(message)s',
                                datefmt='%Y-%m-%d %H:%M:%S')

        elif step_name == 'ICEBERG_ANALYZER':
            logging.basicConfig(format=f'%(asctime)s - ICEBERG_ANALYZER - %(levelname)s: %(message)s',
                                datefmt='%Y-%m-%d %H:%M:%S',
                                level=logging.DEBUG, force=True,
                                handlers=[logging.FileHandler(log_file_path, mode=writing_mode)])
                                          # , logging.StreamHandler()])
        else:
            logging.basicConfig(format=f'%(asctime)s - ICEBERG_ANALYZER - {step_name} - %(levelname)s: %(message)s',
                                datefmt='%Y-%m-%d %H:%M:%S',
                                level=logging.DEBUG, force=True,
                                handlers=[logging.FileHandler(log_file_path, mode=writing_mode)])
                                          #, logging.StreamHandler()])
    except Exception as e:
        logging.error('couldnt config logging basic configuration.')
        raise e


def log_subprocess(proc: subprocess.Popen) -> None:
    """
    log_subprocess:
    Adds the process console outputs to the log.

    :param proc: The process which created using subprocess.Popen.
    """
    output = '\n'
    for line in proc.stdout:
        output += line
    proc.wait()
    if output != '\n':
        logging.info(output)


def get_progress_bar(step_name, experiment_name="", total=0, colour='#DEFFCD', out_of="", bar_length=25):
    spaces = (50 - len(step_name) - len(experiment_name)) * " "
    # bar_format = l_bar + bar + r_bar (with our styling)
    bar_format = step_name + ": " + experiment_name + spaces + "{percentage:3.0f}%" + \
                 "|{bar:" + str(bar_length) + "}| " + \
                 "running " + out_of + " {n_fmt} out of {total_fmt} [{elapsed}<{remaining}]" #, {rate_fmt}{postfix}# for iteration per second

    if step_name == "ICEBERG_ANALYZER":
        return tqdm(total=total,
                    colour='#D3F5FE',
                    bar_format=bar_format
                    )
                    #smoothing=0.5)

    return tqdm(total=total,
                colour=colour,
                bar_format=bar_format + 25 * " ")
                #smoothing=0.5,)