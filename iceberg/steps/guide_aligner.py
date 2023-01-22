"""
This Python module is used to run the GUIDERNA_ALIGNMENT step on treatment and control experiments.

The alignment of the GuideRNA with each site is done by opening window around the highest cut position candidate.
"""

import settings
from utils import aligners, csv_to_df, logs, validation

import pandas as pd
from pyfaidx import Fasta
from pathlib import Path
import logging
import numpy as np
from Bio import Align
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Range1d
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot
from bokeh.resources import CDN
from bokeh.embed import file_html
from bokeh.models import NumeralTickFormatter
import panel as pn

pn.extension()


STEP_NAME = 'GUIDERNA_ALIGNMENT'
STEP_OUT_FOLDER = 'GUIDERNA_ALIGNMENT_RESULTS'
STEP_LOG_FILE = 'GUIDERNA_ALIGNMENT.log'


class GuideAligner(object):

    def __init__(self,
                 guide: str,
                 genome_path: Path,
                 in_folder: Path,
                 out_folder: Path) -> None:
        """
        Initiate the GuideAligner object with the guideRNA and the genome reference file (opened with pyfaidx.Fasta).

        :param guide: The experiment guideRNA.
        :param genome_path: The path to the (indexed) genome reference file.
        :param in_folder: The folder where the input csv file is stored.
        :param out_folder: The folder where the output files directory will be saved.
        """
        self.guide = guide
        self.genome = Fasta(str(genome_path))
        self.in_folder = in_folder
        self.out_folder = out_folder

    def get_colors(self,
                   seqs: list) -> list:
        """
        make colors for bases in sequence.

        :param seqs: the sequence to color.

        :return: list of colors (as strings) for each bp in the given sequence.
        """
        text = [i for s in list(seqs) for i in s]
        clrs = {'A': 'green', 'T': 'red', 'G': 'orange', 'C': 'blue',
                '-': 'white', 'N': 'white', ' ': 'white', '✄': 'white'}
        colors = [clrs[i if i == '✄' else i.upper()] for i in text]
        return colors

    def get_cut_str(self,
                    strand: str,
                    grna_aligned: str) -> str:
        """
        Return the cut site by the alignment.
        (if strand is +: 3 bps before the end of the guideRNA, if strand is -: 3 bps after the start of the guideRNA) .

        :param strand: The chromosome of the iceberg site.
        :param grna_aligned: A string represent the aligned guideRNA (with proceeding and preceding spaces).

        :return: the estimated cut position of the CRISPR activity site (at the specific chromosome).
        """
        if strand == '+':
            cut_str = ''.join(
                np.hstack((len(grna_aligned.rstrip(' ')[:-4]) * [' '],
                          ['✄', ' ', ' ', ' '],
                          (len(grna_aligned) - len(grna_aligned.rstrip())) * [' ']))
            )
        else:
            cut_str = ''.join(
                np.hstack(([c for c in grna_aligned.rstrip(' ') if c == ' '],
                           [' ', ' ', '✄'],
                           (len(grna_aligned) - len([c for c in grna_aligned.rstrip(' ') if c == ' ']) - 3) * [' ']))
            )

        return cut_str

    def view_alignment(self,
                       alignments: list,
                       strand: str,
                       fontsize: str = "9pt",
                       plot_width: int = 1200,
                       start: int = 0) -> gridplot:
        """
        Generate Bokeh sequence alignment view as an html file.

        :param alignments: The alignments of the guideRNA to the iceberg site sequence (using Align.PairwiseAligner).
        :param strand: The chromosome of the iceberg site.
        :param fontsize: The start index (in the chromosome) of the iceberg site.
        :param plot_width: The end index (in the chromosome) of the iceberg site.
        :param start: The start index (in the chromosome) of the iceberg site.

        :return: the visualization as Bokeh gridplot.
        """
        try:
            # make sequence and id lists from the aln object
            seqs = []
            ids = []
            test_str = ''
            for i, a in enumerate(alignments):
                if i == 0:
                    seqs.append(''.join([' '] * len(a[0])))
                    ids.append(f'alignments score:{a[2]}')
                    alignment_start = sum([1 if c == ' ' else 0 for c in a[1].rstrip(' ')]) - 1
                alignment_start = min(alignment_start,
                                      sum([1 if c == ' ' else 0 for c in a[1].rstrip(' ')]) - 1)

                cut_str = self.get_cut_str(strand, a[1])
                seqs.extend([cut_str, a[1], a[0]])
                test_str += f'{i}, {len(a[0])}, {len(a[1])}, {len(cut_str)};'
                ids.extend([f'alignment{i+1}', 'guideRNA sequence' + str(i+1), 'site sequence' + str(i+1)])
            logging.info(test_str)
            seqs.reverse()
            ids.reverse()
            text = [i for s in list(seqs) for i in s]
            colors = self.get_colors(seqs)
            N = len(seqs[0])
            S = len(seqs)
            width = .4
            x = np.arange(start, start + N)
            y = np.arange(0, S, 1)
            # creates a 2D grid of coords from the 1D arrays
            xx, yy = np.meshgrid(x, y)
            # flattens the arrays
            gx = xx.ravel()
            gy = yy.flatten()
            # use recty for rect coords with an offset
            recty = gy + .5
            h = 1 / S
            # now we can create the ColumnDataSource with all the arrays
            source = ColumnDataSource(dict(x=gx,
                                           y=gy,
                                           recty=recty,
                                           text=text,
                                           colors=colors))

            plot_height = len(seqs) * 15 + 50
            x_range = Range1d(start, start+N, bounds='auto')

            if N > 100:
                viewlen = 100
            else:
                viewlen = N
            # view_range is for the close up view
            view_range = (start, start + viewlen+1)

            tools = "xpan, xwheel_zoom, box_zoom, reset, save"

            # entire sequence view (no text, with zoom)
            p = figure(title=None,
                       plot_width=plot_width,
                       plot_height=50,
                       x_range=x_range,
                       y_range=(0, S), tools=tools,
                       min_border=0,
                       toolbar_location='below')

            rects = Rect(x="x",
                         y="recty",
                         width=1,
                         height=1,
                         fill_color="colors",
                         line_color=None,
                         fill_alpha=0.6)

            p.add_glyph(source, rects)
            p.yaxis.visible = False
            p.grid.visible = False
            p.xaxis.formatter = NumeralTickFormatter(format="00")

            # sequence text view with ability to scroll along x axis
            p1 = figure(title=None,
                        plot_width=plot_width,
                        plot_height=plot_height,
                        x_range=view_range,
                        y_range=ids,
                        tools="xpan,reset",
                        min_border=0,
                        toolbar_location='below')  # , lod_factor=1)

            glyph = Text(x="x",
                         y="y",
                         text="text",
                         text_align='center',
                         text_color="black",
                         text_font_size=fontsize)

            rects = Rect(x="x",
                         y="recty",
                         width=1,
                         height=1,
                         fill_color="colors",
                         line_color=None,
                         fill_alpha=0.4)

            p1.add_glyph(source, glyph)
            p1.add_glyph(source, rects)

            p1.grid.visible = False
            p1.xaxis.formatter = NumeralTickFormatter(format="00")
            p1.xaxis.major_label_text_font_style = "bold"
            p1.yaxis.minor_tick_line_width = 0
            p1.yaxis.major_tick_line_width = 0

            p = gridplot([[p], [p1]],
                         toolbar_location='below')
            return p
        except Exception as e:
            logging.error(f'Failed to generate guideRNA site sequence Alignment html file', exc_info=False)
            raise e

    def generate_site_guide_aligment_html(self,
                                          alignments: Align.PairwiseAlignments,
                                          chrom: str,
                                          start: int,
                                          visualization_dir: Path,
                                          strand: str,
                                          out_file_name: str) -> None:
        """
        Generate the alignment visualization as a separate html file.

        :param alignments: The alignments of the guideRNA to the iceberg site sequence (using Align.PairwiseAligner).
        :param chrom: The chromosome of the iceberg site.
        :param start: The start index (in the chromosome) of the iceberg site.
        :param visualization_dir: The directory for writing the html visualization files.
        :param strand: The genome strand (+ or -).
        :param out_file_name: The result file name.

        """
        out_dir = validation.set_curr_out_dir(visualization_dir, '')
        alignment_sequences = []
        row_length = 0
        for i, alignment in enumerate(alignments):
            if i < 5:
                alignment_lines = str(alignment).split('\n')
                aligned_grna_seq = alignment_lines[0]
                aligned_sequence = alignment_lines[2]
                alignment_longest_row = max(len(aligned_sequence), len(aligned_grna_seq))
                if alignment_longest_row > row_length:
                    row_length = alignment_longest_row
                alignment_sequences.append([aligned_grna_seq, aligned_sequence, alignment.score])

            else:
                break
        for alignment_sequence in alignment_sequences:
            for i in range(len(alignment_sequence) - 1):
                row = alignment_sequence[i]
                if len(row) < row_length:
                    alignment_sequence[i] = row + ''.join([' '] * (row_length - len(row)))
        p = self.view_alignment(alignment_sequences, strand=strand, start=start)
        html = file_html(p, CDN, "my plot")
        with open(out_dir / out_file_name, 'w') as f:
            f.write(html)

    def levenshtein(self,
                    seq1: str,
                    seq2: str) -> int:
        """
        Calculate the levenshtein distance between two sequences.

        :param seq1:  The first sequence from the alignment.
        :param seq2: The second sequence from the alignment.

        :return: The Levenshtein distance.
        """
        size_x = len(seq1) + 1
        size_y = len(seq2) + 1
        matrix = np.zeros((size_x, size_y))
        for x in range(size_x):
            matrix[x, 0] = x
        for y in range(size_y):
            matrix[0, y] = y

        for x in range(1, size_x):
            for y in range(1, size_y):
                if seq1[x - 1] == seq2[y - 1]:
                    matrix[x, y] = min(matrix[x - 1, y] + 1,
                                       matrix[x - 1, y - 1],
                                       matrix[x, y - 1] + 1)
                else:
                    matrix[x, y] = min(matrix[x - 1, y] + 1,
                                       matrix[x - 1, y - 1] + 1,
                                       matrix[x, y - 1] + 1)

        return int(matrix[size_x - 1, size_y - 1])

    def get_alignment_fields(self,
                             alignments: Align.PairwiseAlignments,
                             start: int,
                             strand: str) -> (str, str, str, int, int):
        """
        Given the alignments results of one CRISPR activity site (using Align.PairwiseAligner), processing the important
        information about the site.

        :param alignments: The alignments of the guideRNA to the iceberg site sequence (using Align.PairwiseAligner).
        :param start: The start index (in the chromosome) of the iceberg site.
        :param strand: The genome strand (+ or -).

        :return: An Important information about the alignment and the CRISPR activity site.
         (cut_position, short_gRNA_alignment, short_site_alignment, Hamming_distance, Levenshtein_distance).
        """
        try:
            for i, alignment in enumerate(alignments):
                if i == 0:
                    alignment_lines = str(alignment).split('\n')
                    # if strand == '-':
                    #     for i in range(len(alignment_lines)):
                    #         alignment_lines[i] = alignment_lines[i][::-1]
                    site_sequence_aligned = alignment_lines[0]
                    alignment_symbols = alignment_lines[1]
                    grna_aligned = alignment_lines[2]
                    cut_position = start + sum(1 for i in self.get_cut_str(strand, grna_aligned).rstrip() if i == ' ')
                    short_grna_alignment = grna_aligned.strip(' ')
                    short_site_alignment = site_sequence_aligned[len([c for c in grna_aligned.rstrip(' ') if c == ' ']):len([c for c in grna_aligned.rstrip(' ') if c == ' ']) + len(grna_aligned.strip(' '))]
                    hamming_distance = sum(1 for c1, c2 in zip(short_grna_alignment.upper(), short_site_alignment.upper()) if c1 != c2)
                    levenshtein_distance = self.levenshtein(short_grna_alignment.replace('-', '').upper(),
                                                            short_site_alignment.replace('-', '').upper())
                else:
                    break
            return cut_position, short_grna_alignment, short_site_alignment, hamming_distance, levenshtein_distance
        except Exception as e:
            logging.error('Failed to align guideRNA to site sequence', exc_info=False)
            raise e

    def align(self,
              site_sequence: str,
              chrom: str,
              start: int,
              highest_candidate: list,
              visualization_dir: Path,
              out_file_name: str) -> (str, str, str, str, int, int):
        """
        align the guideRNA with the site sequence.

        :param site_sequence: The iceberg site sequence (from the reference genome).
        :param chrom: The chromosome of the iceberg site.
        :param start: The start index (in the chromosome) of the iceberg site.
        :param highest_candidate: The highest cut position candidate.
        :param visualization_dir: The directory for writing the html visualization files.
        :param out_file_name: The result file name.

        :return: An important information about the alignment and the CRISPR activity site.
         (cut_position, strand, short_gRNA_alignment, short_site_alignment, Hamming_distance, Levenshtein_distance).
        """
        try:
            aligner = aligners.get_aligner(mode='local')
            alignments = aligner.align(str(site_sequence).upper(), self.guide)
            alignments_complement = aligner.align(str(site_sequence.complement).upper(), self.guide[::-1])

            if len(alignments) == 0:
                alignments_score = 0
            else:
                alignments_score = alignments[0].score

            if len(alignments_complement) == 0:
                alignments_complement_score = 0
            else:
                alignments_complement_score = alignments_complement[0].score

            if alignments_score == alignments_complement_score == 0:
                return highest_candidate, -1, 'none', 'none', 'none', -1, -1

            if alignments_score >= alignments_complement_score:
                strand = '+'
                best_alignments = alignments
            else:
                strand = '-'
                best_alignments = alignments_complement

            res = self.get_alignment_fields(best_alignments,
                                            start,
                                            strand)

            cut_position, short_grna_alignment, short_site_alignment, hamming_distance, levenshtein_distance = res
            self.generate_site_guide_aligment_html(best_alignments,
                                                   chrom,
                                                   start,
                                                   visualization_dir,
                                                   strand,
                                                   out_file_name)
            return (highest_candidate,
                    cut_position,
                    strand,
                    short_grna_alignment,
                    short_site_alignment,
                    hamming_distance,
                    levenshtein_distance)

        except Exception as e:
            logging.error('Failed to align guideRNA to site sequence', exc_info=False)
            raise e

    def align_guide_to_site(self,
                            icebergs_site: pd.Series,
                            visualization_dir: Path) -> (str, str, str, str, int, int):
        """
        Reads the iceberg site sequence from the genome reference and align the guideRNA with the site sequence.

        :param icebergs_site: Iceberg from the crispr activities DataFrame.
        :param visualization_dir: The directory for writing the html visualization files.

        :return: An Important information about the alignment and the CRISPR activity site.
         (cut_position, strand, short_gRNA_alignment, short_site_alignment, Hamming_distance, Levenshtein_distance).
        """
        chromosome = icebergs_site['chromosome']
        start = icebergs_site['start index']
        end = icebergs_site['end index']
        # cut_position_candidates = ast.literal_eval(icebergs_site["cut-position-candidates t"])
        cut_position_candidates = icebergs_site["cut-position-candidates t"]

        if len(cut_position_candidates) == 0:
            return -1, -1, -1, 'none', 'none', -1, -1

        highest_candidate = max(cut_position_candidates,
                                key=lambda candidate: sum(cut_position_candidates[candidate][primer]
                                                          for primer in cut_position_candidates[candidate]))

        window_size = 20
        window_start = highest_candidate - window_size
        window_end = highest_candidate + window_size
        site_sequence = self.genome[chromosome][window_start:window_end]
        logging.info(f'align guide to {chromosome}:{window_start}-{window_end}')
        out_file_name = f'{chromosome}-{start}-{end}.html'
        logging.info(f'site window sequence: {str(site_sequence)}')
        logging.info(f'experiment GuideRNA: {self.guide}')

        return self.align(site_sequence,
                          chromosome,
                          window_start,
                          highest_candidate,
                          visualization_dir,
                          out_file_name)

    def align_guide_to_sites(self,
                             files_dir: Path,
                             csv_ca: Path) -> Path:
        """
        Iterate over the icebergs sites within the CRISPR activities dataframe and align the guideRNA with the
        sites sequence, and add to the icebergs sites additional important information about the alignments.

        :param files_dir: The directory where the classified icebergs csv files from treatment and control are stored.
        :param csv_ca: The crispr activity csv file name.

        :return: The Path to the new alignments txt file.
        """
        logging.info(f'guideRNA: {self.guide}')

        cuts_position = []
        short_grna_alignments = []
        short_site_alignments = []
        hamming_distances = []
        levenshtein_distances = []
        strands = []
        peaks_t = []

        icebergs_sites_df = csv_to_df.get_merged_experiments_df(files_dir / csv_ca)

        pbar = logs.get_progress_bar(step_name=STEP_NAME,
                                     total=icebergs_sites_df.shape[0],
                                     out_of="icebergs-site")

        for i in icebergs_sites_df.index:
            res = self.align_guide_to_site(icebergs_sites_df.iloc[i, :],
                                           self.out_folder / 'Alignments')
            max_score_peak, cut_position, strand, short_grna_alignment, \
            short_site_alignment, hamming_distance, levenshtein_distance = res

            cuts_position.append(cut_position)
            strands.append(strand)
            short_grna_alignments.append(short_grna_alignment)
            short_site_alignments.append(short_site_alignment)
            hamming_distances.append(hamming_distance)
            levenshtein_distances.append(levenshtein_distance)
            peaks_t.append(max_score_peak)
            pbar.update(1)
        icebergs_sites_df['cut position'] = cuts_position
        icebergs_sites_df['strand'] = strands
        icebergs_sites_df['peak t'] = peaks_t
        icebergs_sites_df['gRNA alignment'] = short_grna_alignments
        icebergs_sites_df['short site alignment'] = short_site_alignments
        icebergs_sites_df['Hamming distance'] = hamming_distances
        icebergs_sites_df['Levenshtein distances'] = levenshtein_distances
        icebergs_sites_df.to_csv(files_dir / csv_ca, index=False)
        return self.out_folder / 'Alignments'


def validate_arguments(args: dict) -> None:
    """
    Validates the step arguments from the given arguments dict.

    :param args: The arguments dict (see main function description).
    """
    try:
        validation.validate_path_existence(args['curr_in_dir'])
        if args['csv_ca'] == '':
            raise ValueError()
        else:
            validation.validate_path_existence(Path(args['OUTPUT_FOLDER_PATH'], args['csv_ca']))
    except ValueError as e:
        logging.error(f'not all necessary arguments for {STEP_NAME} step where giving.')
        raise e


def prepare_step(args: dict) -> None:
    """
    Prepares the step by - setting the logging configurations, setting the output directory and
    validate the step arguments.

    :param args: A dict of arguments (as keys) and their values (see main function description)
    """
    logs.set_logging_config(STEP_NAME, args['OUTPUT_FOLDER_PATH'] / settings.LOGS_FOLDER_RELATIVE_PATH / STEP_LOG_FILE)
    args['curr_out_dir'] = validation.set_curr_out_dir(args['OUTPUT_FOLDER_PATH'], STEP_OUT_FOLDER)
    validate_arguments(args)


def main(args: dict) -> None:
    """
    Run the GUIDERNA_ALIGNMENT step on treatment and control.

    :param args: The pipeline arguments dictionary - dictionary of arguments (as keys) and their values.
    """
    try:
        ga = GuideAligner(args['EXPERIMENTS']['TX']['GUIDERNA'],
                          args['EXPERIMENTS']['GENERAL']['REFERENCE_GENOME_PATH'],
                          args['OUTPUT_FOLDER_PATH'],
                          args['curr_out_dir'])
        args['guide_alignments_dir'] = ga.align_guide_to_sites(args['OUTPUT_FOLDER_PATH'],
                                                               args['csv_ca'])
        # args['curr_in_dir'] = directories_operations.set_curr_in_dir(args['out_dir'], aligner_dir)
        logging.info(f'{STEP_NAME} step done.')
    except Exception as e:
        logging.error(f'Failed to run {STEP_NAME}, please read the following for more information.',
                      exc_info=False)
        # sys.exit()
        raise e
