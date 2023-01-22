from abc import ABC, abstractmethod
from pathlib import Path

@property
@abstractmethod
def STEP_NAME1(self):
    pass


@STEP_NAME1.getter
def STEP_NAME(self):
    pass


@STEP_NAME1.setter
def STEP_NAME(self):
    pass


# class Step(ABC):
#     pass
    #
    # @abstractmethod
    # def prepare_step(args: dict) -> None:
    #     pass
    #
    # @abstractmethod
    # def validate_arguments(args: dict) -> None:
    #     pass
    #
    # @abstractmethod
    # def run(args: dict) -> None:
    #     pass

class Step(object):

    def __init__(self,
                 in_folder: Path,
                 out_folder: Path):
        self.in_folder = in_folder
        self.out_folder = out_folder
