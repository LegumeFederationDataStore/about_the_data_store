#!/usr/bin/env python

import os
import sys
import logging
import re
import subprocess
import gzip
import click
from ruamel.yaml import YAML, RoundTripDumper
from ..helper import check_file, return_filehandle, create_directories,\
                     check_subprocess_dependencies, setup_logging


def read_yaml(template):
    '''Reads yaml file and checks for valid yaml'''
    yaml = YAML()
    target_yaml = yaml.load(open(template, 'rt'))
    return yaml, target_yaml


def generate_readme(template, attributes, logger):
    '''Writes target YAML file'''
    yaml, my_yaml = read_yaml(template)
    my_readme = '{}/README.{}.yaml'.format(attributes['target_dir'],
                                           attributes['key'])
    identifier = attributes['key']
    sci_name = '{} {}'.format(attributes['genus'].capitalize(), 
                                     attributes['species'].lower())
    sci_name_abr = attributes['gensp'].lower()
    genotype = attributes['infra_id']
    my_yaml['identifier'] = identifier  # set key
    my_yaml['genotype'] = genotype  # set infraspecific id
    my_yaml['scientific_name'] = sci_name  # Genus species
    my_yaml['scientific_name_abbrev'] = sci_name_abr  # gensp
    readme_handle = open(my_readme, 'w')
    yaml.dump(my_yaml, readme_handle)
#    for m in my_yaml:
    readme_handle.close()


@click.command()
@click.option('--target', type=str,
              help='''Formatted Directory from Previous Normalizer Steaps''')
@click.option('--template', type=str,
              help='''YAML file to be used as the template for the README''')
@click.option('--log_file', metavar = '<FILE>', 
              default='./normalizer_index.log',
              help='''File to write log to. (default:./normalizer_index.log)''')
@click.option('--log_level', metavar = '<LOGLEVEL>', default='INFO',
              help='''Log level: DEBUG, INFO, WARNING, ERROR, CRITICAL (default:INFO)''')
def cli(target, template, log_file, log_level):
    '''Determines what typae of index to apply to input target'''
    logger = setup_logging(log_file, log_level)
    if not (target and template):
        logger.error('--target and --template arguments are required')
        sys.exit(1)
    target = os.path.abspath(target)  # get full path
    if not os.path.isdir(target):
        logger.error('target dmust be a directory')
        sys.exit(1)
    organism_dir = os.path.basename(os.path.dirname(target))
    target_dir = os.path.basename(target)
    organism_attributes = organism_dir.split('_')  # Genus_species
    target_attributes = target_dir.split('.')
    if len(organism_attributes) != 2:
        logger.error('Parent directory {} is not Genus_species'.format(organism_dir))
        sys.exit(1)
    if len(target_attributes) < 3:
        logger.error('Target directory {} is not delimited correctly'.format(target_dir))
        sys.exit(1)
    gensp = '{}{}'.format(organism_attributes[0][:3], organism_attributes[1][:2])  # abbreviation in yaml README file
    attributes = {'genus': organism_attributes[0],
                  'species': organism_attributes[1],
                  'gensp': gensp,
                  'key': target_attributes[-1],
                  'infra_id': target_attributes[0],
                  'target_dir': target}
#    print('{}\t{}\t{}'.format(template, attributes, logger))
    generate_readme(template, attributes, logger)  # write readme template
