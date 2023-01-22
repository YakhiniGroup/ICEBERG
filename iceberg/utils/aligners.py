"""
This Python module is used to generate aligners with defaults configurations.
"""

from Bio import Align


def get_aligner(mode: str = 'global'):
    """
    generated aligner using Bio.Align.PairwiseAligner.

    :parm mode:'global' for global alignments or 'local' for local alignments.

    :return: Aligner generated using Bio.Align.PairwiseAligner
    """
    aligner = Align.PairwiseAligner(match=30,
                                    mismatch=-6,
                                    mode=mode,
                                    query_left_open_gap_score=0,
                                    query_left_extend_gap_score=0,
                                    target_left_open_gap_score=-20,
                                    target_left_extend_gap_score=-18,
                                    query_internal_open_gap_score=-20,
                                    query_internal_extend_gap_score=-18,
                                    target_internal_open_gap_score=-30,
                                    target_internal_extend_gap_score=-27,
                                    query_right_open_gap_score=0,
                                    query_right_extend_gap_score=0,
                                    target_right_open_gap_score=-20,
                                    target_right_extend_gap_score=-18)
    return aligner
