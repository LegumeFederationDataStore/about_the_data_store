#!/usr/bin/env python

import os
import sys
import gzip
import click
import logging
from incongruency_detector.Detector import Detector


@click.command()
@click.option('--target', multiple=True,
              help='''TARGETS can be files or directories or both''')
@click.option('--no_busco', is_flag=True,
              help='''TARGETS can be files or directories or both''')
@click.option('--log_file', metavar = '<FILE>', 
              default='./detect_incongruencies.log',
              help='''File to write log to. (default:./detect_incongruencies.log)''')
@click.option('--log_level', metavar = '<LOGLEVEL>', default='INFO',
              help='''Log level: DEBUG, INFO, WARNING, ERROR, CRITICAL (default:INFO)''')
#@click.option('--normalize', is_flag=True,
#             help='''Normalizes provided files.

#Incongruencies in FASTA will be corrected if the provided genome name
#passes checks.
    
#The gff file will be tidied if it fails gff3 validation in gt:

#    gt gff3 -sort -tidy -retainids input.gff3 > out.gff3

#''')
def main(target, no_busco, log_file, log_level):
    '''incongruency_detector: 
    
        Detect Incongruencies with LIS Data Store Standard

        https://github.com/LegumeFederation/datastore/PUT_THING_HERE
    '''
    if not target:
        print('Must specify at least one --target.  Run with --help for usage')
        sys.exit(1)
    options = {'log_level': log_level, 'log_file': log_file,
               'no_busco': no_busco}
    for t in target:
        detector = Detector(t, **options)  # initializers
        detector.detect_incongruencies()  # class method runs all deteciton methods


if __name__ == '__main__':
    main()
