"""
This Python module is used to parse the user command line arguments and the input.yaml file.
"""

from utils import validation

import argparse
import yaml
from pathlib import Path

args = None


def parse():
    global args
    if args is None:
        parse_input_yaml_file()

    return args


def parse_input_yaml_file() -> tuple:
    """
    Parses all steps arguments (step_name as key, the value is true if user specify --step_name, else the value is false).

    :return: A dict of the steps and their values.
    """
    global parser
    global args
    parser = argparse.ArgumentParser(description="preform analysis by steps for treatment and control.\n",
                                     add_help=True)

    parser.add_argument('--input_file_path',
                        help="Absolute path to iceberg input.yaml file. See README.md file for more information.")

    parser.add_argument('--enable_profiler',  action="store_true",
                        help="If True, performance profiler is enable using cProfiler.")

    cmd_args, _ = parser.parse_known_args()

    with open(Path(cmd_args.input_file_path)) as input_file_path:
        input_file_fields = yaml.safe_load(input_file_path)

    validation.validate_input_file_fields(input_file_fields)
    args = input_file_fields

    args['enable_profiler'] = cmd_args.enable_profiler
    return args


