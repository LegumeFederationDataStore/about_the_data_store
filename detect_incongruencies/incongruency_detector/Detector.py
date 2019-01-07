#!/usr/bin/env python

import os
import sys
import logging
import re
import subprocess
import hashlib
import requests
import json
import importlib
import importlib.util
from glob import glob
from . import specification_checks
from .helper import check_file, return_filehandle

class Detector:
    '''Class to detect datastore file incongruencies with

       https://github.com/LegumeFederation/datastore
    '''
    def __init__(self, target, **kwargs):
        '''Check for logger, check for gt'''
        subprocess.check_call('gt -help', shell=True)  # check for gt in env
        self.options = kwargs
        self.checks = {}  # object that determines which checks are skipped
        self.checks['genome_main'] = kwargs.get('disable_genome_main', True)
        self.checks['gene_models_main'] = kwargs.get(
                                            'disable_gene_models_main', True)
        self.checks['perform_gt'] = kwargs.get('disable_gt', True)
        self.checks['fasta_headers'] = kwargs.get(
                                            'disable_fasta_headers', True)
        self.canonical_types = ['genome_main', 
                                'gene_models_main', 
                                'ADDMORESTUFF']  #  types for detector
        self.canonical_parents = {'genome_main': None,
                                  'gene_models_main': 'genome_main'}
        self.log_level = kwargs.get('log_level', 'INFO')
        self.log_file = kwargs.get('log_file', './incongruencies.log')
        self.output_prefix = kwargs.get('output_prefix', './incongruencies')
        self.target_objects = {}  # store all target pairings self.get_targets
        self.fasta_ids = {}
        self.reporting = {}
        self.target = os.path.abspath(target)
        self.target_readme = ''
        self.setup_logging()  # setup logging sets self.logger
        logger = self.logger
        self.get_target_type()  # sets self.target_name and self.target_type
        if not self.target_type:  # target type returned False not recognized
            logger.error('Target type not recognized for {}'.format(
                                                                  self.target))
            sys.exit(1)
        logger.info('Target type looks like {}'.format(self.target_type))
        self.get_targets()
#        logger.debug(''.format(self.target_objects))
        logger.info('Performing Checks for the Following:\n')
        for t in self.target_objects:  # for each object set validate
            logger.debug('{}'.format(t))
            logger.debug('{}'.format(self.target_objects[t]))
            logger.info('Parent {}:'.format(t))
            for c in self.target_objects[t]['children']:
                logger.info('Child {}'.format(c))
#            for c in self.target_objects[t]['children']:
#                logger.info('{}'.format(c))
        logger.info('Initialized Detector\n')

    def setup_logging(self):
        '''Return logger based on log_file and log_level'''
        log_level = getattr(logging, self.log_level.upper(), logging.INFO)
        msg_format = '%(asctime)s|%(name)s|[%(levelname)s]: %(message)s'
        logging.basicConfig(format=msg_format, datefmt='%m-%d %H:%M',
                            level=log_level)
        log_handler = logging.FileHandler(self.log_file, mode='w')
        formatter = logging.Formatter(msg_format)
        log_handler.setFormatter(formatter)
        logger = logging.getLogger('detect_incongruencies')
        logger.addHandler(log_handler)
        self.logger = logger

    def get_target_type(self):
        '''Determines what type of target this is, is it a file or directory?

           if its a directory, is it an organism directory or a data directory?
        '''
        logger = self.logger
        target_name = os.path.basename(self.target)
        self.target_name = target_name
        self.target_type = False
        logger.debug(self.target_name)
        if target_name.endswith('.gz'):  # all datastore files end in .gz
            self.target_type = 'file'
        elif(len(target_name.split('_')) == 2 and \
             len(target_name.split('.')) < 4):
            self.target_type = 'organism_dir'  # will always be Genus_species
        elif len(target_name.split('.')) >= 4:  # standard naming minimum
            self.target_type = 'data_dir'

    def get_targets(self):
        '''Gets and discovers target files relation to other files

           if the target is a directory, the program will discover 
           
           all related files that can be checked
        '''
        logger = self.logger
        target = self.target
        target_type = self.target_type
        if target_type == 'file':  # starting with a file
            self.add_target_object()
            return
        elif target_type == 'data_dir' or target_type == 'organism_dir':
            self.get_all_files()  # works for both data and organism
            return
        #elif target_type == 'organism_dir':  # entire organism
        #    self.get_all_files()
        #    return

    def get_all_files(self):
        '''Walk down filetree and recursively return all files
        
           These will be checked using add_target_object
        '''
        logger = self.logger
        target = self.target
        for root, directories, filenames in os.walk(target):
            for filename in filenames:  # we only care about the files
                my_target = os.path.join(root, filename)  # make path
                logger.info('Checking file {}'.format(my_target))
                self.target = my_target  # set target
                self.target_name = os.path.basename(self.target)
                self.add_target_object()  # add target if canonical

    def add_target_object(self):
        '''Uses parent child logic to create a datastructure for objects'''
        logger = self.logger
        target = self.target
        logger.debug(target)
        if not check_file(target):
            logger.error('Could not find file: {}'.format(target))
            sys.exit(1)
        organism_dir = os.path.dirname(os.path.dirname(target))  # org dir
        target_attributes = self.target_name.split('.')
        if len(target_attributes) < 3:
            logger.error('File {} does not seem to have attributes'.format(
                                                                      target))
            return
        target_format = target_attributes[-2]  # get gff, fna, faa all gx
        canonical_type = target_attributes[-3]  # check content type
        if canonical_type not in self.canonical_types:  # regject
            logger.error('Type {} not recognized in {}.  Skipping'.format(
                                                          canonical_type,
                                                    self.canonical_types))
            return
        target_key = target_attributes[-4]  # get key
        target_ref_type = self.canonical_parents[canonical_type]
        logger.info('Getting target files reference if necessary...')
        if len(target_attributes) > 7 and target_ref_type:  # check parent
            logger.info('Target Derived from Some Reference Searching...')
            ref_glob = '{}/{}*/*{}.*.gz'.format(organism_dir, 
                                      '.'.join(target_attributes[1:3]),
                                      target_ref_type)
            my_reference = self.get_reference(ref_glob)
            if my_reference not in self.target_objects:
                self.target_objects[my_reference] = {'type': target_ref_type,
                                                     'readme': '',
                                                     'children': {}}
                self.target_objects[my_reference]['children'][target] = {
                                                    'type': canonical_type}
            else:
                if target not in self.target_objects[my_reference]['children']:
                    self.target_objects[my_reference]['children'][target] = {
                                                     'type': canonical_type}
        else:
            logger.info('Target has no Parent, it is a Reference')
            if not target in self.target_objects:
                self.target_objects[target] = {'type': canonical_type,
                                               'children': {}}

    def get_reference(self, glob_target):
        '''Finds the FASTA reference for some prefix'''
        logger = self.logger
        if len(glob(glob_target)) > 1:  # too many references....?
            logger.error('Multiple references found {}'.format(glob_target))
            sys.exit(1)
        reference = glob(glob_target)
        if not reference:  # if the objects parent could not be found
            logger.error('Could not find ref glob'.format(glob_target))
            sys.exit(1)
        reference = glob(glob_target)[0]
        if not os.path.isfile(reference):  # if cannot find reference file
            logger.error('Could not find main target {}'.format(reference))
            sys.exit(1)
        logger.info('Found reference {}'.format(reference))
        return reference

    def detect_incongruencies(self):
        '''Detects all incongruencies in self.target_objects
        
           report all incongruencies to self.output_prefix

           Probably refactor this at some point to use a self.dict

           This is to allow just iteration and lookup to occur and remove

           basically all this conditional foo
        '''
        logger = self.logger
        targets = self.target_objects  # get objects from class init
        for reference in targets:
            self.target = reference
            ref_method = getattr(specification_checks, 
                                 targets[reference]['type'])
            if not ref_method:
                logger.error('Check for {} does not exist'.format(
                                                targets[reference]['type']))
                continue
            logger.debug(ref_method)
            my_detector = ref_method(self, **self.options)
            my_detector.run()
            logger.debug('{}'.format(targets[reference]))
            if self.target_objects[reference]['children']:  # process children
                children = self.target_objects[reference]['children']
                for c in children:
                    logger.info('Performing Checks for {}'.format(c))
                    self.target = c
                    child_method = getattr(specification_checks,
                                           children[c]['type'])
                    if not child_method:
                        logger.error('Check for {} does not exist'.format(
                                                         children[c]['type']))
                        continue
                    logger.debug(child_method)
                    my_detector = child_method(self, **self.options)
                    my_detector.run()
                    logger.info('{}'.format(c))
