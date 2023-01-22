from iceberg import analyzer, settings
from utils import validation

import unittest
import os
import yaml
import filecmp
from pathlib import Path
import logging
from utils import logs

test_control_dir = Path(os.getcwd()) / "test" / "TIGIT-2_18-4-21_OUT_IB_NEW_TEST"


def are_files_sizes_alike(file1, file2):
    size1 = os.path.getsize(file1)
    size2 = os.path.getsize(file2)

    if abs(size1 - size2) < 0.1 * size1:
        return True
    return False


def are_dir_trees_files_sizes_alike(dir1, dir2):
    """
    Compare two directories recursively. Files in each directory are
    assumed to be equal if their names and contents are equal.

    :param dir1: First directory path
    :param dir2: Second directory path

    :return: True if the directory trees are the same and
        there were no errors while accessing the directories or files,
        False otherwise.
   """
    for subdir, dirs, files in os.walk(dir1):
        for file in files:
            file_path = subdir + os.sep + file
            if not are_files_sizes_alike(file1=file_path,
                                         file2=file_path.replace(str(dir1), str(dir2))):
                return False
    return True


def are_dir_trees_equal(dir1, dir2):
    """
    Compare two directories recursively. Files in each directory are
    assumed to be equal if their names and contents are equal.

    @param dir1: First directory path
    @param dir2: Second directory path

    @return: True if the directory trees are the same and
        there were no errors while accessing the directories or files,
        False otherwise.
   """

    dirs_cmp = filecmp.dircmp(dir1, dir2)
    if len(dirs_cmp.left_only) > 0 or len(dirs_cmp.right_only) > 0 or len(dirs_cmp.funny_files) > 0:
        return False
    (_, mismatch, errors) = filecmp.cmpfiles(dir1,
                                             dir2,
                                             dirs_cmp.common_files,
                                             shallow=False)
    if len(mismatch) > 0 or len(errors) > 0:
        return False
    for common_dir in dirs_cmp.common_dirs:
        new_dir1 = os.path.join(dir1, common_dir)
        new_dir2 = os.path.join(dir2, common_dir)
        if not are_dir_trees_equal(new_dir1, new_dir2):
            return False
    return True


class TestIceberg(unittest.TestCase):

    def test_analyzer(self):
        with open(Path(os.getcwd()) / "test" / "test_input.yaml") as file:
            args = yaml.safe_load(file)
        validation.validate_input_file_fields(args)
        analyzer.main(args)

        pbar = logs.get_progress_bar(step_name='ICEBERG_ANALYZER_TESTS',
                                     total=len(args["ANALYZER STEPS"]),
                                     out_of='folders-comparison')

        are_dirs_equal = are_dir_trees_equal(test_control_dir / settings.umi.STEP_OUT_FOLDER,
                                             args["OUTPUT_FOLDER_PATH"] / settings.umi.STEP_OUT_FOLDER)
        self.assertTrue(are_dirs_equal, f'{settings.umi.STEP_OUT_FOLDER} not as expected')
        pbar.update(1)

        are_dirs_equal = are_dir_trees_equal(test_control_dir / settings.bwa.STEP_OUT_FOLDER,
                                             args["OUTPUT_FOLDER_PATH"] / settings.bwa.STEP_OUT_FOLDER)
        self.assertTrue(are_dirs_equal, f'{settings.bwa.STEP_OUT_FOLDER} not as expected')
        pbar.update(1)

        are_dirs_equal = are_dir_trees_equal(test_control_dir / settings.experiment_traces_remover.STEP_OUT_FOLDER,
                                             args["OUTPUT_FOLDER_PATH"] / settings.experiment_traces_remover.STEP_OUT_FOLDER)
        self.assertTrue(are_dirs_equal, f'{settings.experiment_traces_remover.STEP_OUT_FOLDER} not as expected')
        pbar.update(1)

        are_dirs_alike = are_dir_trees_files_sizes_alike(test_control_dir / settings.bwa.STEP_OUT_FOLDER,
                                                         args["OUTPUT_FOLDER_PATH"] / settings.bwa.STEP_OUT_FOLDER)
        self.assertTrue(are_dirs_alike, f'{settings.bwa.STEP_OUT_FOLDER} not as expected')
        pbar.update(1)

        are_dirs_alike = are_dir_trees_files_sizes_alike(test_control_dir / settings.reads_uniter.STEP_OUT_FOLDER,
                                                         args["OUTPUT_FOLDER_PATH"] / settings.reads_uniter.STEP_OUT_FOLDER)
        self.assertTrue(are_dirs_alike, f'{settings.reads_uniter.STEP_OUT_FOLDER} not as expected')
        pbar.update(1)

        are_dirs_alike = are_dir_trees_files_sizes_alike(test_control_dir / settings.experiments_merge_manager.STEP_OUT_FOLDER,
                                                         args["OUTPUT_FOLDER_PATH"] / settings.experiments_merge_manager.STEP_OUT_FOLDER)
        self.assertTrue(are_dirs_alike, f'{settings.experiments_merge_manager.STEP_OUT_FOLDER} not as expected')
        pbar.update(1)

        are_dirs_alike = are_dir_trees_files_sizes_alike(test_control_dir / settings.icebergs_uniter.STEP_OUT_FOLDER,
                                                         args["OUTPUT_FOLDER_PATH"] / settings.icebergs_uniter.STEP_OUT_FOLDER)
        self.assertTrue(are_dirs_alike, f'{settings.icebergs_uniter.STEP_OUT_FOLDER} not as expected')
        pbar.update(1)

        are_dirs_alike = are_dir_trees_files_sizes_alike(test_control_dir / settings.classify.STEP_OUT_FOLDER,
                                                         args["OUTPUT_FOLDER_PATH"] / settings.classify.STEP_OUT_FOLDER)
        self.assertTrue(are_dirs_alike, f'{settings.classify.STEP_OUT_FOLDER} not as expected')
        pbar.update(1)

        are_dirs_alike = are_dir_trees_files_sizes_alike(test_control_dir / settings.sites_profile_calculator.STEP_OUT_FOLDER,
                                                         args["OUTPUT_FOLDER_PATH"] / settings.sites_profile_calculator.STEP_OUT_FOLDER)
        self.assertTrue(are_dirs_alike, f'{settings.sites_profile_calculator.STEP_OUT_FOLDER} not as expected')
        pbar.update(1)

        are_dirs_alike = are_dir_trees_files_sizes_alike(test_control_dir / settings.guide_alignment.STEP_OUT_FOLDER,
                                                         args["OUTPUT_FOLDER_PATH"] / settings.guide_alignment.STEP_OUT_FOLDER)
        self.assertTrue(are_dirs_alike, f'{settings.guide_alignment.STEP_OUT_FOLDER} not as expected')
        pbar.update(1)

        are_dirs_alike = are_dir_trees_files_sizes_alike(test_control_dir / settings.reports.STEP_OUT_FOLDER,
                                                         args["OUTPUT_FOLDER_PATH"] / settings.reports.STEP_OUT_FOLDER)

        self.assertTrue(are_dirs_alike, f'{settings.reports.STEP_OUT_FOLDER} not as expected')
        pbar.update(1)


if __name__ == '__main__':
    try:
        logs.set_logging_config('ICEBERG_ANALYZER_TESTS', None)
        print(settings.START_MESSAGE)
        unittest.main()
    except Exception:
        logging.error(f'Failed to run the iceberg test! :( ', exc_info=False)
