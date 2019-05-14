#!/usr/bin/env python

import os
import sys
import gzip
import click
import logging
from incongruency_detector.Detector import Detector
from incongruency_detector.Normalizer import Normalizer


@click.command()
@click.option('--target', multiple=True,
              help='''TARGETS can be files or directories or both''')
@click.option('--no_busco', is_flag=True,
              help='''Disable BUSCO''')
@click.option('--normalize', is_flag=True,
              help='''Converts input files into data store standard''')
@click.option('--gnm', metavar = '<gnm#>',
              help='''Genome Version for Normalizer. ex. gnm1''')
@click.option('--ann', metavar = '<ann#>',
              help='''Annotation Version for Normalizer. ex. ann1''')
@click.option('--genus', metavar = '<STRING>',
              help='''Genus of organism input file''')
@click.option('--species', metavar = '<STRING>',
              help='''Species of organism input file''')
@click.option('--infra_id', metavar = '<STRING>',
              help='''Line or infraspecific identifier.  ex. A17_HM341''')
@click.option('--unique_key', metavar = '<STRING, len=4>',
              help='''4 Character unique idenfier (get from spreadsheet)''')
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
def main(target, no_busco, normalize, gnm, ann, genus, species, infra_id, unique_key, log_file, log_level):
    '''incongruency_detector:
    
        Detect Incongruencies with LIS Data Store Standard

        https://github.com/LegumeFederation/datastore/PUT_THING_HERE
    '''
    if not target:
        print('Must specify at least one --target.  Run with --help for usage')
        sys.exit(1)
    data_store_home = os.path.dirname(os.path.dirname(
                                                   os.path.abspath(__file__)))
    print(data_store_home)
    sys.exit(1)
    options = {'log_level': log_level, 'log_file': log_file,
               'no_busco': no_busco, 'gnm': gnm, 'ann': ann, 'genus': genus,
               'species': species, 'infra_id': infra_id, 
               'unique_key': unique_key}
    for t in target:
        if normalize:
            normalizer = Normalizer(t, **options)
            t = normalizer.normalize()
        detector = Detector(t, **options)  # initializers
        detector.detect_incongruencies()  # class method runs all deteciton methods


if __name__ == '__main__':
    main()
