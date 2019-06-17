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
@click.option('--gnm',
              help='''Genome Version for Normalizer. ex. 1 = gnm1''')
@click.option('--ann',
              help='''Annotation Version for Normalizer. ex. 1 = ann1''')
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
    templates_dir = '{}/templates'.format(data_store_home)
    readme_template = '{}/template__README.KEY.yml'.format(templates_dir)
    options = {'log_level': log_level, 'log_file': log_file,
               'no_busco': no_busco, 'gnm': gnm, 'ann': ann, 'genus': genus,
               'species': species, 'infra_id': infra_id, 
               'unique_key': unique_key, 'readme_template': readme_template}
    for t in target:
        if normalize:
            normalizer = Normalizer(t, **options)
            t = normalizer.normalize()
        detector = Detector(t, **options)  # initializers
        detector.detect_incongruencies()  # class method runs all deteciton methods

def cli():
    ''']
       Current Tools:
         basic_fasta_stats
         format_fasta
         chunk_fasta
         chunk_fastq
         fastx_converter
         filter_fasta_by_length
         get_fasta_by_id
         get_fastq_by_id
         subset_fastq
       Please run `sequencetools <TOOL> --help` for individual usage
    '''
    if not len(sys.argv) > 1:
        print(cli.__doc__)
        sys.exit(1)
    tool = sys.argv[1].lower()  # get tool
    sys.argv = [tool] + sys.argv[2:]  # make usage for applicaiton correct
    if tool == '--help' or tool =='-h':
        print(cli.__doc__)
        sys.exit(1)
    if tool == 'basic_fasta_stats':
        from .tools import basic_fasta_stats
        basic_fasta_stats.main()
    if tool == 'format_fasta':
        from .tools import format_fasta
        format_fasta.main()
    if tool == 'chunk_fastq':
        from .tools import chunk_fastq
        chunk_fastq.main()
    if tool == 'chunk_fasta':
        from .tools import chunk_fasta
        chunk_fasta.main()
    if tool == 'fastx_converter':
        from .tools import fastx_converter
        fastx_converter.main()
    if tool == 'filter_fasta_by_length':
        from .tools import filter_fasta_by_length
        filter_fasta_by_length.main()
    if tool == 'get_fasta_by_id':
        from .tools import get_fasta_by_id
        get_fasta_by_id.main()
    if tool == 'get_fastq_by_id':
        from .tools import get_fastq_by_id
        get_fastq_by_id.main()
    if tool == 'subset_fastq':
        from .tools import subset_fastq
        subset_fastq.main()

if __name__ == '__main__':
    cli()
#    main()
