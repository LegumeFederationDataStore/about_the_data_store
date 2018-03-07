#!/usr/bin/env python

import os
import sys
import argparse
import gzip
import logging
from incongruency_detector.Detector import Detector

parser = argparse.ArgumentParser(description='''
    Detect and optionally Normalize Incongruencies with LIS Data Store Standard

    https://github.com/LegumeFederation/datastore/issues/23

''', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--genome', metavar = '</path/to/my/genome.fna>',
help='''Genome to check''')

parser.add_argument('--annotation', metavar = '</path/to/my/annotation.gff>',
help='''Annotation to check''')

parser.add_argument('--directory', metavar = '</path/to/my/data_dir>',
help='''Directory to check.  Will auto detect type.''')

parser.add_argument('--gt_path', metavar = '</path/to/my/genome_tools/>',
help='''Path to genome tools''', required=True)

parser.add_argument('--log_file', metavar = '<FILE>', 
default='./detect_incongruencies.log',
help='''File to write log to.  (default:./detect_incongruencies.log)''')

parser.add_argument('--log_level', metavar = '<LOGLEVEL>', default='INFO',
help='''Log level: DEBUG, INFO, WARNING, ERROR, CRITICAL (default:INFO)''')

parser.add_argument('--normalize', action='store_true',
help='''Normalizes provided files.

    Incongruencies in FASTA will be corrected if the provided genome name
    passes checks.
    
    The gff file will be tidied if it fails gff3 validation in gt:

        gt gff3 -sort -tidy -retainids input.gff3 > out.gff3

''')

parser._optionals.title = "Program Options"
args = parser.parse_args()

log_file = args.log_file

log_level = getattr(logging, args.log_level.upper(), logging.INFO)
msg_format = '%(asctime)s|%(name)s|[%(levelname)s]: %(message)s'
logging.basicConfig(format=msg_format, datefmt='%m-%d %H:%M',
                    level=log_level)
log_handler = logging.FileHandler(
                       log_file, mode='w')
formatter = logging.Formatter(msg_format)
log_handler.setFormatter(formatter)
logger = logging.getLogger('detect_incongruencies')
logger.addHandler(log_handler)


if __name__ == '__main__':
    genome = args.genome
    annotation = args.annotation
    directory = args.directory
    gt_path = args.gt_path
    normalize = args.normalize
    initializers = {'genome': genome, 'annotation': annotation,
                    'directory': directory,
                    'logger': logger, 'gt_path': gt_path,
                    'normalize': normalize}
    detector = Detector(**initializers)
    detector.detect_incongruencies()

