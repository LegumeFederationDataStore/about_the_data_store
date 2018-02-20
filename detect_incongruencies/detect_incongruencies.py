#!/usr/bin/env python

import os
import sys
import argparse
import gzip
import logging
from incongruency_detector.Detector import Detector

parser = argparse.ArgumentParser(description='''
    Detect Incongruencies with LIS Data Store Standard

    https://github.com/LegumeFederation/datastore/issues/23

''', formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('--genome', metavar = '</path/to/my/genome.fna>',
help='''Genome to check''')

parser.add_argument('--annotation', metavar = '</path/to/my/annotation.gff>',
help='''Annotation to check''')

parser.add_argument('--gt_path', metavar = '</path/to/my/genome_tools/>',
help='''Path to genome tools''')

parser.add_argument('--log_level', metavar = '<LOGLEVEL>', default='INFO',
help='''Log level: DEBUG, INFO, WARNING, ERROR, CRITICAL (default:INFO)''')

parser.add_argument('--tidy', action='store_true',
help='''Tidy gff file using:

    gt gff3 -sort -tidy -retainids input.gff3 > out.gff3

''')

parser._optionals.title = "Program Options"
args = parser.parse_args()

log_level = getattr(logging, args.log_level.upper(), logging.INFO)
msg_format = '%(asctime)s|%(name)s|[%(levelname)s]: %(message)s'
logging.basicConfig(format=msg_format, datefmt='%m-%d %H:%M',
                    level=log_level)
log_handler = logging.FileHandler(
                       './detect_incongruencies.log', mode='w')
formatter = logging.Formatter(msg_format)
log_handler.setFormatter(formatter)
logger = logging.getLogger('detect_incongruencies')
logger.addHandler(log_handler)


if __name__ == '__main__':
    genome = args.genome
    annotation = args.annotation
    gt_path = args.gt_path
    tidy = args.tidy
    targets = {'genome' : genome, 'annotation' : annotation,
               'logger' : logger, 'gt_path' : gt_path}
    optionals = {'tidy' : tidy}
    detector = Detector(**targets)
    detector.detect_incongruencies(**optionals)

