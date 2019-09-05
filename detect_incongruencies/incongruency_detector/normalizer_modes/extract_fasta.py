#!/usr/bin/env python

import os
import sys
import logging
import re
import subprocess
import gzip
import click
from ruamel.yaml import YAML, RoundTripDumper
from sequencetools.helpers import sequence_helpers
from ..helper import check_file, return_filehandle, create_directories,\
                     check_subprocess_dependencies, setup_logging


def primary_transcript_check(peptides, logger):
    '''Decides if protein file form run_gffread is primary of not

       if the file is not primary, it contains isoforms, the longest

       is selected for each gene and a protein_primaryTranscript file

       is produced.
    '''
    seq_handle = open(peptides, 'rt')
    longest = {}
    count = 0
    primary = '{}.protein_primaryTranscript.faa'.format('.'.join(peptides.split('.')[:-2]))
    for record in sequence_helpers.get_seqio_fasta_record(seq_handle):
        count += 1
        my_record = '{}\n{}'.format(record.id, record.seq)
        my_obj = {'record': my_record, 'length': len(record.seq)}
        gene_id = '.'.join(record.id.split('.')[:-1])
        if gene_id not in longest:
            longest[gene_id] = my_obj
        else:
            if len(record.seq) > longest[gene_id]['length']:
                longest[gene_id] = my_obj
    if len(longest) == count:
        logger.info('Protein file is already primary transcripts')
        return  # protein file is already primary
    primary_handle = open(primary, 'w')
    for r in longest:
        primary_handle.write(longest[r]['record'] + '\n')
    primary_handle.close()
    return primary


def run_gffread(gff, reference, logger):
    '''Reads gff3 file and writes mRNA, CDS and Peptides'''
    gff_dir = os.path.dirname(gff)
    gff_attributes = os.path.basename(gff).split('.')
    mrna = '{}/{}.mRNA.fna'.format(gff_dir, '.'.join(gff_attributes[:5]))
    cds = '{}/{}.cds.fna'.format(gff_dir, '.'.join(gff_attributes[:5]))
    pep = '{}/{}.protein.faa'.format(gff_dir, '.'.join(gff_attributes[:5]))
    cmd = 'gffread {} -g {} -w {} -x {} -y {} -W'.format(gff, reference, mrna,
                                                         cds, pep)
#    subprocess.check_call(cmd, shell=True)
    primary_transcript_check(pep, logger)


@click.command()
@click.option('--target', type=str,
              help='''Formatted file from normalizer prefix MUST BE UNCOMPRESSED''')
@click.option('--reference', type=str,
              help='''FASTA file of sequences to read gff intervals normalizer prefix is a good start here too MUST BE UNCOMPRESSED''')
@click.option('--log_file', metavar = '<FILE>', 
              default='./normalizer_index.log',
              help='''File to write log to. (default:./normalizer_index.log)''')
@click.option('--log_level', metavar = '<LOGLEVEL>', default='INFO',
              help='''Log level: DEBUG, INFO, WARNING, ERROR, CRITICAL (default:INFO)''')
def cli(target, reference, log_file, log_level):
    '''Determines what typae of index to apply to input target'''
    logger = setup_logging(log_file, log_level)
    check_subprocess_dependencies()
    if not (target and reference):
        logger.error('--target and --reference arguments are required')
        sys.exit(1)
    target = os.path.abspath(target)  # get full path
    target_attributes = os.path.basename(target).split('.')
    if os.path.basename(reference).split('.')[-1] == 'gz':
        logger.error('GFFREAD cannot process compressed fasta as reference')
        sys.exit(1)
    if len(target_attributes) < 7:
        logger.error('Target file {} is not delimited correctly'.format(target))
        sys.exit(1)
    run_gffread(target, reference, logger)  # write readme template
