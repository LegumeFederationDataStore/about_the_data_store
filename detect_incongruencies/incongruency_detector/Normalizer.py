#!/usr/bin/env python

import os
import sys
import logging
import re
import subprocess
import gzip
from ruamel.yaml import YAML
from .helper import check_file, return_filehandle, create_directories


class Normalizer:
    '''Class to normalize datastore file incongruencies with

       https://github.com/LegumeFederation/datastore/issues/23
    '''
    def __init__(self, target, **kwargs):
        '''setup logger or load one.  set genome and annotation.

           set attributes and reporting
        '''
        subprocess.check_call('which gt', shell=True)  # check in env
        subprocess.check_call('which bgzip', shell=True)  # check in env
        subprocess.check_call('which tabix', shell=True)  # check in env
        subprocess.check_call('which samtools', shell=True)  # check in env
        self.target = target
        self.options = kwargs
        self.gnm = kwargs.get('gnm')
        self.ann = kwargs.get('ann')
        self.genus = kwargs.get('genus')
        self.species = kwargs.get('species')
        self.infra_id = kwargs.get('infra_id')
        self.unique_key = kwargs.get('unique_key')
        self.readme_tempalte = kwargs.get('readme_template')
        self.yaml=YAML(typ='safe')   # round-trip default
        self.readme_specification = ''  # will be filled if yaml validates
        if not kwargs.get('logger'):
            msg_format = '%(asctime)s|%(name)s|[%(levelname)s]: %(message)s'
            logging.basicConfig(format=msg_format, datefmt='%m-%d %H:%M',
                                level=logging.INFO)
            log_handler = logging.FileHandler(
                                   './detect_incongruencies.log')
            formatter = logging.Formatter(msg_format)
            log_handler.setFormatter(formatter)
            logger = logging.getLogger('detect_incongruencies')
            logger.addHandler(log_handler)
            self.logger = logger
        else:
            self.logger = kwargs.get('logger')
        if not self.gnm or not self.species or not self.genus or not self.unique_key or not self.infra_id:
            error_message = '''
Minimum requirements for normalizer are:

    gnm, species, genus, infraspecific or line id, and unique_key

    ann is specified if file is gff3
'''
            self.logger.error(error_message)
            sys.exit(1)
        if check_file(self.readme_template):
            logger.info('Checking README Speficiation: {}'.format(
                                                        self.readme_template))
            self.readme_specification = yaml.load(open(self.readme_template, 
                                                       'rt'))
            for k in readme_specification:
                logger.debug('{}: {}'.format(k, readme_specification[k]))
        self.gensp = '{}{}'.format(self.genus[:3].lower(), 
                                   self.species[:2].lower())
#        self.fasta_ids = {}
#        self.genome_attributes = {'filename': '', 'version': '',
#                                  'prefix': '', 'type': '', 'build': '',
#                                  'compression': ''}
#        self.incongruencies = {'Genome Name': [], 'Genome Headers': [],
#                               'Annotation Name': [], 'GFF File': []}
        self.logger.info('Initialized Normalizer')

    def normalize_fasta(self):
        '''Normalizes FASTA files for data store standard'''
        re_header = re.compile("^>(\S+)\s*(.*)")  # match header
        logger = self.logger
        infra_id = self.infra_id
        prefix = self.gensp
        gnm = self.gnm
        key = self.unique_key
        new_file_raw = '{}.{}.{}.{}.genome_main.fna'.format(
                                             prefix, infra_id, gnm, key)
        species_dir = './{}_{}'.format(self.genus.capitalize(), self.species)
        new_file_dir = '{}/{}.{}.{}'.format(species_dir, infra_id, gnm,
                                               key)  # setup for ds
        create_directories(new_file_dir)  # make genus species dir for output
        fasta_file_path = os.path.abspath(
                                  '{}/{}'.format(new_file_dir, new_file_raw))
        new_fasta = open(fasta_file_path, 'w')
        name_prefix = '{}.{}.{}'.format(prefix, infra_id, gnm)
        fh = return_filehandle(self.target)
        with fh as gopen:
            for line in gopen:
                line = line.rstrip()
                if re_header.match(line):  # header line
                    parsed_header = re_header.search(line).groups()
                    logger.debug(line)
                    logger.debug(parsed_header)
                    hid = ''
                    desc = ''
                    new_header = ''
                    if isinstance(parsed_header, str):  # check tuple
                        hid = parsed_header
                        new_header = '>{}.{}'.format(name_prefix, hid)
                        logger.debug(hid)
                    else:
                        if len(parsed_header) == 2:  # header and description
                            hid = parsed_header[0]  # header
                            desc = parsed_header[1]  # description
                            new_header = '>{}.{} {}'.format(name_prefix, hid, 
                                                            desc)
#                    normalized_header = '>'
                            logger.debug(hid)
                            logger.debug(desc)
                            logger.debug(new_header)
                        else:
                            logger.error('header {} looks odd...'.format(line))
                            sys.exit(1)
                    new_fasta.write(new_header + '\n')  # write new header
                else:  # sequence lines
                    new_fasta.write(line + '\n')
        new_fasta.close()
        if not check_file(fasta_file_path):
            logger.error('Output file {} not found for normalize fasta'.format(
                                                              fasta_file_path))
            sys.exit(1)  # new file not found return False
        return fasta_file_path

    def type_rank(self, hierarchy, feature_type):
        if feature_type in hierarchy:
            return hierarchy[feature_type]['rank']
        else:
            return 1000  # no feature id, no children, no order
    
    def update_hierarchy(self, hierarchy, feature_type, parent_types):
        '''Breaks input gff3 line into attributes.
    
           Determines feature type hierarchy using type and attributes fields
        '''
        if not parent_types:
            if feature_type not in hierarchy:
                hierarchy[feature_type] = {'children': [], 'parents': [], 
                                           'rank': 0}
            hierarchy[feature_type]['rank'] = 1
        else:
            if feature_type not in hierarchy:
                hierarchy[feature_type] = {'children': [], 
                                           'parents': [], 'rank': 0}
            hierarchy[feature_type]['parents'] = parent_types
            for p in parent_types:
                if hierarchy.get(p):
                    if feature_type not in hierarchy[p]['children']:
                        hierarchy[p]['children'].append(feature_type)
                else:
                    hierarchy[p] = {'children': [feature_type], 'parents': [],
                                    'rank': 0}
    
    def normalize_gff3(self):
        logget = self.logger
        infra_id = self.infra_id
        prefix = self.gensp
        gnm = self.gnm
        ann = self.ann
        key = self.unique_key
        new_file_raw = '{}.{}.{}.{}.{}.gene_models_main.gff3'.format(
                                             prefix, infra_id, gnm, ann, key)
        species_dir = './{}_{}'.format(self.genus.capitalize(), self.species)
        new_file_dir = '{}/{}.{}.{}.{}'.format(species_dir, infra_id, gnm,
                                               ann, key)  # setup for ds
        create_directories(new_file_dir)  # make genus species dir for output
        gff_file_path = os.path.abspath(
                                  '{}/{}'.format(new_file_dir, new_file_raw))
        new_gff = open(gff_file_path, 'w')
        get_id = re.compile('ID=([^;]+)')  # gff3 get id string
        get_name = re.compile('Name=([^;]+)')  # get name string
        get_parents = re.compile('Parent=([^;]+)')  # get parents
        gff3_lines = []  # will be sorted after loading and hierarchy
        sub_tree = {}  # feature hierarchy sub tree
        type_hierarchy = {}  # type hierarchy for features, will be ranked
        prefix_name = False
        fh = return_filehandle(self.target)
        with fh as gopen:
            for line in gopen:
                line = line.rstrip()
                if not line:
                    continue
                if line.startswith('#'):  # header
                    if line.startswith('##sequence-region'):  # replace header
                        fields = re.split('\s+', line)
                        ref_name = '{}.{}.{}.{}'.format(prefix, infra_id,
                                                        gnm, fields[1])
                        line = re.sub(fields[1], ref_name, line)
                    new_gff.write('{}\n'.format(line))
                    continue
                fields = line.split('\t')
                gff3_lines.append(fields)
                feature_id = get_id.search(fields[-1])
                parent_ids = get_parents.search(fields[-1])
                if parent_ids:
                    parent_ids = parent_ids.group(1).split(',')
                else:
                    parent_ids = []
                if not feature_id:
                    continue
                feature_id = feature_id.group(1)
                sub_tree[feature_id] = {'type': fields[2],  # parent child
                                        'parent_ids': parent_ids}
        for f in sub_tree:
            feature_type = sub_tree[f]['type']
            parent_types = []
            for p in sub_tree[f]['parent_ids']:
                parent_type = sub_tree[p]['type']
                if not parent_type:  # must have type
                    logger.error('could not find type for parent {}'.format(
                                                                parent_type))
                    sys.exit(1)  # error no type
                parent_types.append(parent_type)
            self.update_hierarchy(type_hierarchy, feature_type, parent_types)
        ranking = 1  # switch to stop ranking in while below
        while ranking:
            check = 0  # switch to indicate no features unranked
            for t in sorted(type_hierarchy, key=lambda k:type_hierarchy[k]['rank'], reverse=True):
                if type_hierarchy[t]['rank']:  # rank !=0
                    for c in type_hierarchy[t]['children']:
                        if not type_hierarchy[c]['rank']:
                            type_hierarchy[c]['rank'] = type_hierarchy[t]['rank'] + 1
                else:  # rank == 0
                    check = 1
            if not check:  # escape loop when all features are ranked
                ranking = 0
        feature_prefix = '{}.{}.{}.{}'.format(prefix, infra_id, gnm, ann)
        for l in sorted(gff3_lines, key=lambda x:(x[0], int(x[3]), self.type_rank(type_hierarchy, x[2]))):  # rank by chromosome, start, type_rank and stop
            l[0] = '{}.{}.{}.{}'.format(prefix, infra_id, gnm, l[0])  # rename
            l = '\t'.join(l)
            feature_id = get_id.search(l)
            feature_name = get_name.search(l)
            feature_parents = get_parents.search(l)
            if feature_id:  # if id set new id
                new_id = '{}.{}'.format(feature_prefix, feature_id.group(1))
                l = get_id.sub('ID={}'.format(new_id), l)
            if feature_name and prefix_name:  # if name and flag set new name
                new_name = '{}.{}'.format(feature_prefix, feature_name.group(1))
                l = get_name.sub('Name={}'.format(new_name), l)
            if feature_parents:  # parents set new parent ids
                parent_ids = feature_parents.group(1).split(',')
                new_ids = []
                for p in parent_ids:  # for all parents
                    new_id = '{}.{}'.format(feature_prefix, p)
                    new_ids.append(new_id)
                new_ids = ','.join(new_ids)
                l = get_parents.sub('Parent={}'.format(new_ids), l)
            new_gff.write('{}\n'.format(l))
        new_gff.close()
        if not check_file(gff_file_path):
            logger.error('Output file {} not found for normalize gff'.format(
                                                               gff_file_path))
            sys.exit(1)  # new file not found return False
        return gff_file_path  # return path to new gff file

    def index_fasta(self, fasta):
        '''Indexes the supplied fasta file using bgzip and samtools'''
        logger = self.logger
        logger.info('compressing {}'.format(fasta))
        cmd = 'bgzip -f --index {}'.format(fasta)  # bgzip command
        subprocess.check_call(cmd, shell=True)
        compressed_fasta = '{}.gz'.format(fasta)  # will now have gz after bgzip
        logger.info('Indexing {}'.format(compressed_fasta))
        cmd = 'samtools faidx {}'.format(compressed_fasta)  # samtools faidx
        subprocess.check_call(cmd, shell=True)
        return compressed_fasta  # return new compressed and indexed gff

    def index_gff3(self, gff3):
        '''Indexes the supplied gff3 file using bgzip and tabix'''
        logger = self.logger
        logger.info('compressing {}'.format(gff3))
        cmd = 'bgzip -f {}'.format(gff3)  # bgzip command
        subprocess.check_call(cmd, shell=True)
        compressed_gff = '{}.gz'.format(gff3)  # will now have gz after bgzip
        logger.info('Indexing {}'.format(compressed_gff))
        cmd = 'tabix -p gff {}'.format(compressed_gff)  # tabix index command
        subprocess.check_call(cmd, shell=True)
        return compressed_gff  # return new compressed and indexed gff

    def normalize(self):
        '''Determines what type of normalization to apply to input target'''
        file_types = ['fna', 'fasta', 'fa', 'gff', 'gff3', 'ADD MORE']
        fasta = ['fna', 'fasta', 'fa']
        gff3 = ['gff', 'gff3']
        logger = self.logger
        target = self.target
        file_attributes = target.split('.')
        new_file = False
        if len(file_attributes) < 2:
            error_message = '''Target {} does not have a type or attributes.  File must end in either gz, bgz, fasta, fna, fa, gff, or gff3.'''.format(target)
            logger.error(error_message)
            sys.exit(1)
        file_type = file_attributes[-1]
        if file_type == 'gz' or file_type == 'bgz':
            file_type = file_attributes[-2]
        if file_type not in file_types:
            logger.error('File {} is not a type in {}'.format(target, 
                                                              file_types))
            sys.exit(1)
        if file_type in fasta:
            logger.info('Target is a FASTA file formatting and indexing...')
            new_file = self.normalize_fasta()  # returns path to new file
            logger.info('Normalizing done, indexing {}'.format(new_file))
            new_file = self.index_fasta(new_file)
            logger.info('Indexing done, final file: {}'.format(new_file))
        if file_type in gff3:
            logger.info('Target is a gff3 file formatting and indexing...')
            new_file = self.normalize_gff3()  # returns path to new file
            logger.info('Sorting and Normalizing done, indexing {}'.format(
                                                                    new_file))
            new_file = self.index_gff3(new_file)
            logger.info('Indexing done, final file: {}'.format(new_file))
        if not new_file:
            logger.error('Normalizer FAILED.  See Log.')
            sys.exit(1)
        return new_file
