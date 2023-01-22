"""
This Python module contains the iceberg pipeline settings -validation scheme, paths, commands templates,
strings templates, html templates and etc.
"""

from steps import umi_wrapper as umi
from steps import experiment_libraries_detector
from steps import experiment_traces_remover
from steps import bwa_wrapper as bwa
from steps import reads_uniter
from steps import experiments_merge_manager
from steps import icebergs_uniter
from steps import break_sites_classifier as classify
from steps import sites_profile_calculator
from steps import guide_aligner as guide_alignment
from steps import reports_manager as reports

import os
from pathlib import Path

PROJECT_ROOT_DIRECTORY = Path(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

STEPS_METADATA = [umi, experiment_libraries_detector, experiment_traces_remover, bwa, reads_uniter,
                  experiments_merge_manager, icebergs_uniter, classify, sites_profile_calculator,
                  guide_alignment,
                  reports]

LOGS_FOLDER_RELATIVE_PATH = Path('DEBUG', 'LOGS')


########################################################################################################################
#                                                                                                                      #
#                                                  validation                                                          #
#                                                                                                                      #
########################################################################################################################

ICEBERG_ARGUMENTS_SCHEMA = {
                         'OUTPUT_FOLDER_PATH': {
                             'required': True,
                             'type': 'string'
                            },

                         'ANALYZER STEPS': {
                            'type': 'dict',
                            'schema': {
                                'UMI': {
                                    'required': True,
                                    'type': 'boolean',
                                    'allowed': [True]
                                },
                                'EXPERIMENT_LIBRARIES_DETECTION': {
                                    'required': True,
                                    'type': 'boolean',
                                    'allowed': [True]
                                },
                                'EXPERIMENT_TRACES_REMOVAL': {
                                    'required': True,
                                    'type': 'boolean',
                                    'allowed': [True]
                                },
                                'BWA': {
                                    'required': True,
                                    'type': 'boolean',
                                    'allowed': [True]
                                },
                                'UNITE_READS_TO_ICEBERGS': {
                                    'required': True,
                                    'type': 'boolean',
                                    'allowed': [True]
                                },
                                'MERGE_TREATMENT_AND_CONTROL_ICEBERGS_BY_LOCUS': {
                                    'required': True,
                                    'type': 'boolean',
                                    'allowed': [True]
                                },
                                'UNITE_CLOSE_ICEBERGS_SITES': {
                                    'required': True,
                                    'type': 'boolean',
                                    'allowed': [True]
                                },
                                'BREAKS_CLASSIFY': {
                                    'required': True,
                                    'type': 'boolean',
                                    'allowed': [True]
                                },
                                'CALCULATE_ICEBERGS_SITES_PROFILE': {
                                    'required': True,
                                    'type': 'boolean',
                                    'allowed': [True]
                                },
                                'GUIDERNA_ALIGNMENT': {
                                    'required': True,
                                    'type': 'boolean',
                                    'allowed': [True]
                                },
                                'REPORTS': {
                                    'required': True,
                                    'type': 'boolean',
                                    'allowed': [True]
                                }
                            }
                         },

                         'EXPERIMENTS': {
                            'type': 'dict',
                            'schema': {
                                'GENERAL': {
                                    'type': 'dict',
                                    'schema': {
                                        'REFERENCE_GENOME_PATH': {
                                            'required': True,
                                            'type': 'string',
                                        },
                                        'EXPERIMENTS_TAG': {
                                            'required': True,
                                            'type': 'string',
                                            'minlength': 10,
                                            'maxlength': 100
                                        }
                                    }
                                },
                                'TX': {
                                    'type': 'dict',
                                    'schema': {
                                        'NAME': {
                                            'required': True,
                                            'type': 'string',
                                            'regex': r'.*[a-zA-Z0-9].*'
                                        },
                                        'EXPERIMENT_FOLDER_PATH': {
                                            'required': True,
                                            'type': 'string'
                                        },
                                        'R1': {
                                            'required': True,
                                            'type': 'string'
                                        },
                                        'R2': {
                                            'required': True,
                                            'type': 'string'
                                        },
                                        'I1': {
                                            'required': True,
                                            'type': 'string'
                                        },
                                        'I2': {
                                            'required': True,
                                            'type': 'string'

                                        },
                                        'GUIDERNA': {
                                            'required': True,
                                            'type': 'string'
                                        }
                                    }
                                },
                                'CONTROL': {
                                    'type': 'dict',
                                    'schema': {
                                        'NAME': {
                                            'required': True,
                                            'type': 'string',
                                            'regex': r'.*[a-zA-Z0-9].*'

                                        },
                                        'EXPERIMENT_FOLDER_PATH': {
                                            'required': True,
                                             'type': 'string'
                                        },
                                        'R1': {
                                            'required': True,
                                            'type': 'string'
                                        },
                                        'R2': {
                                            'required': True,
                                            'type': 'string'
                                        },
                                        'I1': {
                                            'required': True,
                                            'type': 'string'
                                        },
                                        'I2': {
                                            'required': True,
                                            'type': 'string'
                                        }
                                    }
                                }
                            }
                        },

                        'HYPERPARAMATERS': {
                            'type': 'dict',
                            'schema': {
                                'UMI_BPS_AMOUNT_FROM_READS_START': {
                                    'required': True,
                                    'type': 'integer',
                                    'min': 0,
                                    'max': 30
                                },
                                'MIN_QUALITY': {
                                    'required': True,
                                    'type': 'integer',
                                    'min': 0,
                                    'max': 60
                                },
                                'MIN_FREQUENCY': {
                                    'required': True,
                                    'type': 'float',
                                    'min': 0.01,
                                    'max': 0.999
                                },
                                'MAX_READS_DISTANCE': {
                                    'required': True,
                                    'type': 'integer',
                                    'min': 0,
                                    'max': 10000
                                },
                                'MAX_ICEBERG_DISTANCE': {
                                    'required': True,
                                    'type': 'integer',
                                    'min': 0,
                                    'max': 20000
                                },
                                'MAX_ALIGNMENTS_HAMMING_DISTANCE': {
                                    'required': True,
                                    'type': 'integer',
                                    'min': 0,
                                    'max': 34
                                },
                                'NOISE_BINS_AND_CONTROL_MAPQ_PERCENTILE': {
                                    'type': 'list',
                                    'schema': {
                                        'type': 'list',
                                        'items': [{'type': 'list',
                                                   'items': [{'type': 'integer', 'required': True, 'min': 0, 'max': 60},
                                                             {'type': 'integer', 'required': True, 'min': 0, 'max': 61}]},
                                                  {'type': 'float', 'required': True, 'min': 0.01, 'max': 0.999},
                                                  {'type': 'integer', 'required': True, 'min': 0, 'max': 1000}]
                                    }
                                },
                                'CRISPR_ACTIVITY_THRESHOLD': {
                                    'required': True,
                                    'type': 'float',
                                    'min': 0.01,
                                    'max': 1
                                }
                            }
                        }
                    }


########################################################################################################################
#                                                                                                                      #
#                                                   ANALYZER                                                           #
#                                                                                                                      #
########################################################################################################################
STARTING_READ_AMOUNT_TEMPLATE = \
    '''
    starting reads amount (before {step}):
                                     {description}
    '''

READ_AMOUNT_TEMPLATE = \
    '''
    remaining reads amount after {step}:
                                     {description}
    '''


########################################################################################################################
#                                                                                                                      #
#                                        EXPERIMENT_LIBRARIES_DETECTION                                                #
#                                                                                                                      #
########################################################################################################################
TAG_AND_TAG_RC_COMPARISON_TEMPLATE = \
    '''
read name r2:{read_name}
forward-library primer alignment score: {alignments_forward_score}            
{alignments_forward}
    
reverse-library primer alignment score: {alignments_reverse_score}
{alignments_reverse}
    '''

########################################################################################################################
#                                                                                                                      #
#                                            EXPERIMENT_TRAILS_REMOVAL                                                 #
#                                                                                                                      #
########################################################################################################################
TAG_AND_TAG_RC_COMPARISON_TEMPLATE = \
    '''
read name r2:{read_name}
forward-library primer alignment score: {alignments_forward_score}            
{alignments_forward}

reverse-library primer alignment score: {alignments_reverse_score}
{alignments_reverse}
    '''


########################################################################################################################
#                                                                                                                      #
#                                                     BWA                                                              #
#                                                                                                                      #
########################################################################################################################
INDEX_GENOME_CMD_TEMPLATE = 'bwa index -a bwtsw "{genome_ref}"'

INDEX_FASTQ_CMD_TEMPLATE = 'bwa index "{file_path}"'

BWA_CMD_TEMPLATE = 'bwa mem -M -t 2 "{genome_ref}" "{r1_fastq_path}" "{r2_fastq_path}" > "{out_sam_path}"'

SORT_CMD_TEMPLATE = 'samtools sort -O "{out_format}" "{in_file_path}" -o "{out_file_path}"'

VIEW_SAM_TO_BAM_CMD_TEMPLATE = 'samtools view -S -b "{sam_path}" > "{bam_path}"'

INDEX_BAM_CMD_TEMPLATE = 'samtools index "{bam_path}"'

########################################################################################################################
#                                                                                                                      #
#                                              UNITE_READS_TO_ICEBERGS                                                 #
#                                                                                                                      #
########################################################################################################################

########################################################################################################################
#                                                                                                                      #
#                                   MERGE_TREATMENT_AND_CONTROL_ICEBERGS_BY_LOCUS                                      #
#                                                                                                                      #
########################################################################################################################

########################################################################################################################
#                                                                                                                      #
#                                              UNITE_CLOSE_ICEBERGS_SITES                                              #
#                                                                                                                      #
########################################################################################################################

########################################################################################################################
#                                                                                                                      #
#                                                     CLASSIFY                                                         #
#                                                                                                                      #
########################################################################################################################


EXPERIMENT_CHROMOSOME_STAT_TEMPLATE = \
    '''    
        ------------------------------------------------------------------------------------
                                        {experiment_name}:
            
            icebergs count: {icebergs_count:.0f}.
            reads count: {reads_count:.0f}.
            
            avg icebergs -
                length: {avg_length:.3f}.     reads count: {avg_reads_count:.3f}.
            
            longest iceberg - 
                length: {longest_length:.0f}.    reads count: {longest_reads_count:.0f}.    index: {longest_start:.0f} - {longest_end:.0f}.
            
            largest read amount iceberg -
                length: {largest_length:.0f}.    reads count: {largest_reads_count:.0f}.    index: {largest_start:.0f} - {largest_end:.0f}.
    '''

########################################################################################################################
#                                                                                                                      #
#                                            SITES_PROFILE_CALCULATOR                                                  #
#                                                                                                                      #
########################################################################################################################
html_string = '''
<html>
  <link rel="stylesheet" type="text/css" href="df_style.css"/>
  <body>
    {table}
  </body>
</html>.
'''

########################################################################################################################
#                                                                                                                      #
#                                                GUIDERNA_ALIGNMENT                                                    #
#                                                                                                                      #
########################################################################################################################


########################################################################################################################
#                                                                                                                      #
#                                                   REPORTS                                                            #
#                                                                                                                      #
########################################################################################################################
CREATE_REPORT_CMD_TEMPLATE = 'create_report "{bedfile}" "{genome}" --tracks {tracks} --output "{out_file}" --type junction ' \
                             '--info-columns Reads_Count_Treatment Reads_Count_Control Depth_Treatment Depth_Control Iceberg_Length_Treatment ' \
                             'Iceberg_Length_Control Mapq_Treatment Mapq_Control {guideRNA_alignment} --title "{html_title}"'

HTML_MAIN_REPORT_TEMPLATE = \
    '''
    <!DOCTYPE html>
    <html>
    <head>
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <style>
    * {{
      box-sizing: border-box;
    }}

    /* Create two equal columns that floats next to each other */
    .column {{
      float: left;
      width: 50%;
      padding: 10px;
    }}

    /* Clear floats after the columns */
    .row:after {{
      content: "";
      display: table;
      clear: both;
    }}
    </style>

    <style>
      body {{
        font-size: 18px;
      }}

      ul {{
        list-style-type: none;
        margin: 0;
        padding: 0;
        overflow: hidden;
        background-color: rgb(211, 211, 211);
        position: -webkit-sticky; /* Safari */
        position: sticky;
        top: 0;
      }}

      li {{
        float: right;

      }}
      li {{
        float: left;
        border-right:1px solid rgb(255, 252, 252);

      }}

      li a {{
        display: block;
        color: rgb(58, 58, 58);
        text-align: center;
        padding: 30px 16px;
        text-decoration: none;
      }}

      li a:hover {{
        background-color: rgb(243, 243, 243);
      }}

      .active {{
        background-color: #ffffff;
      }}
      </style>
    </head>

    <body>
    <ul>
      <li><a class="active" href="#home"><img src="HTML/iceberg-logo.jpg" alt="Trulli" width="35px" height="35px"></a></a></li>
      <li><a href="HTML/Reports/{csv_ca:}" target="_blank">CRISPR Activities Sites Report</a></li>
      <li><a href="HTML/Reports/{csv_sb:}" target="_blank">Spontaneous Breaks Sites Report</a></li>
      <li><a href="HTML/Reports/{csv_noise:}" target="_blank">Noise Sites Report</a></li>
      <li style="float:right"><a href="https://github.com/ilaigenish" target="_blank">About</a></li>
    </ul>

    <div align="center">
      <h1>Iceberg Report</h1>
    </div>

    <div align="center">
      <h2>Reads Funnel Chart</h2>
      <iframe id="reads_count" src = "HTML/Experiment-Overview/reads_count.html" height="500px" width="1000px" frameBorder="0" ></iframe>
    </div>

    <iframe id="reads_count" src = "HTML/Experiment-Overview/sites-amount-in-each-class.html" height="200px" width="1000px" frameBorder="0" ></iframe>

    <div align="center">
      <div class="row">
        <div class="column" style="background-color:rgb(255, 255, 255); border-style: solid; border-color: rgb(230, 229, 229);">

          <h2>Treatment</h2>

          <p>Reads Count Per Sites Types</p>
          <iframe id="pie_reads_t" src = "HTML/Experiment-Overview/treatment-reads-pie-chart.html" height="500px" width="500px" frameBorder="0" ></iframe>
          <p>Sites Count Per Sites Types</p>
          <iframe id="pie_sites_t" src = "HTML/Experiment-Overview/treatment-sites-pie-chart.html" height="500px" width="500px" frameBorder="0" ></iframe>
          <p>Sites Reads Count vs Sites MAPQ</p>
          <iframe id="scatter_read_count_mapq_t" src = "HTML/Experiment-Overview/treatment-reads-count-mapq-scatter-plot.html" height="500px" width="500px" frameBorder="0" ></iframe>
        </div>

        <div class="column" style="background-color:rgb(255, 255, 255); border-style: solid; border-color: rgb(230, 229, 229);">

          <h2>Control</h2>

          <p>Reads Count Per Sites Types</p>
          <iframe id="pie_reads_m" src = "HTML/Experiment-Overview/Control-reads-pie-chart.html" height="500px" width="500px" frameBorder="0" ></iframe>
          <p>Sites Count Per Sites Types</p>
          <iframe id="pie_sites_m" src = "HTML/Experiment-Overview/Control-sites-pie-chart.html" height="500px" width="500px" frameBorder="0" ></iframe>
          <p>Sites Reads Count vs Sites MAPQ</p>
          <iframe id="scatter_read_count_mapq_m" src = "HTML/Experiment-Overview/Control-reads-count-mapq-scatter-plot.html" height="500px" width="500px" frameBorder="0" ></iframe>
        </div>
      </div>
    </div>

    <iframe id="command" src = "HTML/Experiment-Overview/command.html" height="700px" width="1000px" frameBorder="0" ></iframe>

    </body>
    </html>
    '''

HTML_LOGO = \
    '''
    <img src="../iceberg-logo.jpg" alt="Iceberg Logo" width="40" height="40">
    '''

HTML_CENTER_HEADER_TEMPLATE = \
    '''
    <div align="center">
        {line}
        <h2>Sites</h2>
    </div>
    '''

HTML_ALIGNMENTS_IFRAME = \
    '''
    <div align="center">
        <h2>GuideRNA And Site Sequence Alignments</h2>
        <p><b>Note:</b> Coordinates indicate genome positions of reference genome without gaps.</p>
        <iframe id="alignment_html" src = "" height="250px" width="1300px" frameBorder="0" ></iframe>
    </div>
    <div align="center">
        <h2> Treatment Site Profile Histogram</h2>
        <iframe id="profile_histogram_html" src = "" height="550px" width="1300px" frameBorder="0" ></iframe>
    </div>
    <div align="center">
        <h2> Treatment Site Profile Table</h2>
        <iframe id="profile_table_html" src = "" height="550px" width="1300px" frameBorder="0" ></iframe>
    </div>
    '''
HTML_CHANGE_ALIGNMENTS_IFRAME_ON_CLICK_TEMPLATE = \
    '''
        document.getElementById("alignment_html").src = "..{os_path_sep}{html_alignments_dir_name}{os_path_sep}".concat(locus.replace(":", "-")).concat(".html")
        document.getElementById("profile_histogram_html").src = "../sites-profiles/".concat(locus.replace(":", "-")).concat("-histogram.html")
        document.getElementById("profile_table_html").src = "../sites-profiles/".concat(locus.replace(":", "-")).concat("-table.html")

    '''

HTML_SITES_AMOUNTS_TEMPLATE = \
    '''
    <html>
        <body>
            <h2>Number of CRISPR activities sites (on/off targets): {ca_sites_count}</h2>
            <h2>Number of spontaniuos breaks sites: {sb_sites_count}</h2>
            <h2>Number of noise sites: {n_sites_count}</h2>
        </body>
    </html>
    '''

HTML_ICEBERG_CMD_AND_ARGUMENTS_TEMPLATE = \
    '''
    <html>
        <body>
            <h1>User Command For Generate This Report:</h1>
            <p>{cmd}<br><br><p>
            <h1>Pipeline Arguments:</h1>
            <p>{args}<p>
        </body>
    </html>
    '''


########################################################################################################################
#                                                                                                                      #
#                                                   Terminal Styling                                                   #
#                                                                                                                      #
########################################################################################################################

START_MESSAGE = \
    '''
    ========================================================================================================================                                                                                                                                
    
        NNNNNNN    NNNNNNNNNN   NNNNNNNNNNNNNNNNN NNNNNNNNNNNNNNNN   NNNNNNNNNNNNNNNN NNNNNNNNNNNNNNN         NNNNNNNNNN   
        NCAATCN  NATAGACGTGAGCN NTGAGATCATCGAGCAN NCGATAGNNANNCAGCTN NATGCTACAGCTACTN NTCTGAGNTCNAGCGN      NATCAACAGTCN   
        NNGCTNN NCTCANNNNNCGTCN NNTCACGNNNNNTTGCN NNTACGGNNNNAACAGTN NNGTCGANNNNTAGCN NNTGCTNNNNNGCATCN    NAAGANNNNTCTN   
         NTCAN NGCATN    NNNNNN  NTTACTN   NNNNNN  NATGCN    NCTGTGN  NAGATN   NNNNNN  NTAGCN    NTAGTAN  NCCTAN   NNNNN   
         NCTGN NACAGN            NTAGCAN           NGCTTN   NATCACN   NTCGAN           NAGATN   NAGTATN  NATCCN             
         NACTN NCAGAN            NACGATNNNNNN      NTAGCNNNNTGAGTN    NCGATCNNNNN      NCATGNNNNNTGTN    NCCACN             
         NTCTN NAGTAN            NCGGTTNNNNNN      NTATGNNNNNNTACACN  NAGTATNNNNN      NGTATNNNNNACGATCN NTACAN   NNNNCTAN
         NCGAN NAACTN            NATCGTN           NCTAGN     NTAGCTN NGCAGN           NCATAN    NGATCGN NAGTACN    NCCN   
         NTACN  NTCTGN   NNNNNN  NTAAGCN   NNNNNN  NAGTAN    NTACATCN NATAAN   NNNNNN  NACGTN    NAGTGAN  NCGTCGN   NGGN   
        NTAGANN NACTANNNNNGTAAN NNAGAGTNNNNNGCGAN NNCGTACNNNNGACACTN NNGATGTNNNACGTCN NNACAGTN   NGTGACN   NGTACNNNNATTN    
        CAGTTCN  NGACGCTGCGTGTN NAGTCAGTGCAGTAGCN NGTCGAGATCACTGAGN  NCGCATGAATCTGCAN NGACACTN   NCTACTN    NGTATATAGTAN   
        NNNNNNN    NNNNNNNNNN   NNNNNNNNNNNNNNNNN NNNNNNNNNNNNNNN    NNNNNNNNNNNNNNNN NNNNNNNN   NNNNNNN     NNNNNN NNNN   
    
                        Analyzing Tool For GUIDE-seq ExperimentsAnalyzing Tool For GUIDE-seq Experiments ©                                    
    
    ========================================================================================================================
    '''


UPPER_FRAME = \
    '''
    ============================================================================================
    ||                                                                                        ||
    ||                                                                                        ||
    '''


CHROMOSOME_TITLE_TEMPLATE = '\t\t\t\t\t\t\t\tchromosome name: {chromosome_name}'


LOWER_FRAME = \
    '''
    ||                                                                                        ||
    ||                                                                                        ||
    ============================================================================================
    '''
