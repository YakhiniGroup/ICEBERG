"""
This Python module is used to validate paths existence, inputs and etc.
"""

import settings

import logging
import os
from pathlib import Path
from cerberus import Validator
from Bio.Seq import Seq


def validate_input_file_fields(input_file_fields: dict):
    schema = settings.ICEBERG_ARGUMENTS_SCHEMA
    v = Validator(schema)
    v.validate(input_file_fields, schema)
    if len(v.errors) > 0:
        errors_str = "Failed to validate input.yaml file: \n" + \
                     f"{v.errors}"
        raise ValueError(errors_str)
    input_file_fields['OUTPUT_FOLDER_PATH'] = Path(os.path.abspath(input_file_fields['OUTPUT_FOLDER_PATH']))
    _validate_experiments_data(input_file_fields['EXPERIMENTS'])


def _validate_experiments_data(experiments):
    experiments['GENERAL']['REFERENCE_GENOME_PATH'] = Path(experiments['GENERAL']['REFERENCE_GENOME_PATH'])
    validate_path_existence(experiments['GENERAL']['REFERENCE_GENOME_PATH'])

    experiments['GENERAL']['EXPERIMENTS_TAG_RC'] = str(Seq(experiments['GENERAL']['EXPERIMENTS_TAG']).reverse_complement())

    for experiment in experiments.keys() - {'GENERAL'}:
        experiments[experiment]['EXPERIMENT_FOLDER_PATH'] = Path(experiments[experiment]['EXPERIMENT_FOLDER_PATH'])
        validate_path_existence(experiments[experiment]['EXPERIMENT_FOLDER_PATH'])
        for key in experiments[experiment].keys():
            if key in ['R1', 'R2', 'I1', 'I2']:
                path = experiments[experiment]['EXPERIMENT_FOLDER_PATH'] / experiments[experiment][key]
                validate_path_existence(path)


def validate_path_existence(dir_path: Path) -> None:
    """
    Validate the given path existence, throw FileNotFoundError if path not found.

    :param dir_path: A path.
    """
    try:
        if not os.path.exists(dir_path):
            raise FileNotFoundError(f'{dir_path}')

    except FileNotFoundError as e:
        logging.error('files directory not exist: {}'.format(dir_path), exc_info=False)
         # sys.exit()
        raise e


def set_curr_out_dir(base_directory: Path,
                     dir_to_set: str) -> Path:
    """
    Sets a path (and create parts of it if needed).

    :param base_directory: The directory where to create the dir_to_set directory.
    :param dir_to_set: The directory to create.
    """
    curr_out_dir = Path(base_directory, dir_to_set)
    if not os.path.exists(curr_out_dir):
        curr_out_dir.parent.mkdir(parents=True, exist_ok=True)
        curr_out_dir.mkdir(parents=True, exist_ok=True)
        logging.info('directory: {} created in order to store the iceberg results.'.format(curr_out_dir))
    else:
        logging.warning('directory: {} already exists, if you repeats on analysis which you all ready '
                        'done the files will be overwritten.'.format(curr_out_dir))
    return curr_out_dir


def set_curr_in_dir(base_directory: Path,
                    dir_to_set: str):
    """
    Validate the given path existence and return it.

    :param base_directory: The path to directory that contains  the dir_to_set.
    :param dir_to_set: The end directory.
    """
    curr_in_dir = Path(os.path.abspath(base_directory), dir_to_set)
    validate_path_existence(curr_in_dir)
    return curr_in_dir
