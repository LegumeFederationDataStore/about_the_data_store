#!/usr/bin/env python

import os
import sys
import logging
import re
import subprocess
from Normalizer import Normalizer
from file_helpers import check_file, return_filehandle


class Detector:
    '''Class to detect datastore file incongruencies with

       https://github.com/LegumeFederation/datastore/issues/23
    '''
    def __init__(self, **kwargs):
        '''setup logger or load one.  set genome and annotation.

           set attributes and reporting
        '''
        if not kwargs.get('logger'):
            msg_format = '%(asctime)s|%(name)s|[%(levelname)s]: %(message)s'
            logging.basicConfig(format=msg_format, datefmt='%m-%d %H:%M',
                                level=logging.DEBUG)
            log_handler = logging.FileHandler(
                                   './detect_incongruencies.log')
            formatter = logging.Formatter(msg_format)
            log_handler.setFormatter(formatter)
            logger = logging.getLogger('detect_incongruencies')
            logger.addHandler(log_handler)
            self.logger = logger
        else:
            self.logger = kwargs.get('logger')
        self.genome = kwargs.get('genome')
        self.annotation = kwargs.get('annotation')
        self.gt_path = kwargs.get('gt_path')
        self.normalize = kwargs.get('normalize')
        if self.annotation and not self.gt_path:
            logger.error('gt_path, /path/to/genome_tools must be provided')
            sys.exit(1)
        if self.gt_path:
            self.gt_path = os.path.abspath(self.gt_path)
        if self.normalize:
            self.normalizer = Normalizer(**kwargs)
        else:
            self.normalizer = None
        self.fasta_ids = {}
        self.genome_attributes = {'filename': '', 'version': '',
                                  'prefix': '', 'type': '', 'build': '',
                                  'compression': ''}
        self.incongruencies = {'Genome Name': [], 'Genome Headers': [],
                               'Annotation Name': [], 'GFF File': []}
        self.logger.info('Initialized Detector')

    def check_genome_main(self, attr):
        '''accepts a list of genome attributes split by "."

           https://github.com/LegumeFederation/datastore/issues/23

           checks these file attributes to ensure they are correct
        '''
        logger = self.logger
        if len(attr) != 7:  # should be 7 fields
            logger.error('File did not have 7 fields! {}'.format(attr))
            sys.exit(1)
        if len(attr[0]) != 5:  # should be 5 letter prefix
            logger.error('File must have 5 letter prefix, not {}'.format(
                                                                    attr[0]))
            sys.exit(1)
        if not attr[2].startswith('gnm'):  # should be gnm type
            logger.error('File should have gnm in field 3, not {}'.format(
                                                                     attr[2]))
            sys.exit(1)
        gnm_v = attr[2].replace('gnm', '')
        try:
            int(gnm_v)
        except ValueError:  # best way to check for int in python2
            logger.error('gnm version must be integer not {}'.format(gnm_v))
            sys.exit(1)
        if (len(gnm_v.split('.')) > 1):  # check for float
            logger.error('gnm version must be integer not {}'.format(gnm_v))
            sys.exit(1)
        if not attr[5] == 'fna':  # should be fna type
            logger.error('File should be fna not {}'.format(attr[5]))
            sys.exit(1)
        if not attr[6] == 'gz':  # should be gzip compressed
            logger.error('Last field should be gz, not {}'.format(attr[6]))
            sys.exit(1)

    def check_gene_models_main(self, attr):
        '''accepts a list of annotation attributes split by "."

           https://github.com/LegumeFederation/datastore/issues/23

           checks these file attributes to ensure they are correct
        '''
        logger = self.logger
        if len(attr) != 8:  # should be 8 fields
            logger.error('File did not have 7 fields! {}'.format(attr))
            sys.exit(1)
        if len(attr[0]) != 5:  # should be 5 letter prefix
            logger.error('File must have 5 letter prefix, not {}'.format(
                                                                    attr[0]))
            sys.exit(1)
        if not attr[2].startswith('gnm'):  # should be gnm type
            logger.error('File should have gnm in field 3, not {}'.format(
                                                                     attr[2]))
            sys.exit(1)
        gnm_v = attr[2].replace('gnm', '')
        try:
            int(gnm_v)
        except ValueError:  # best way to check for int in python2
            logger.error('gnm version must be integer not {}'.format(gnm_v))
            sys.exit(1)
        if (len(gnm_v.split('.')) > 1):  # check for float
            logger.error('gnm version must be integer not {}'.format(gnm_v))
            sys.exit(1)
        if not attr[3].startswith('ann'):  # should be gnm type
            logger.error('File should have ann in field 4, not {}'.format(
                                                                     attr[2]))
            sys.exit(1)
        ann_v = attr[3].replace('ann', '')
        try:
            int(ann_v)
        except ValueError:  # best way to check for int in python2
            logger.error('ann version must be integer not {}'.format(ann_v))
            sys.exit(1)
        if (len(ann_v.split('.')) > 1):  # check for float
            logger.error('ann version must be integer not {}'.format(gnm_v))
            sys.exit(1)
        if not attr[6] == 'gff3':  # should be gff3 type
            logger.error('File should be gff3 not {}'.format(attr[6]))
            sys.exit(1)
        if not attr[7] == 'gz':  # should be gzip compressed
            logger.error('Last field should be gz, not {}'.format(attr[6]))
            sys.exit(1)

    def check_fasta(self, fasta, attr):
        '''Confirms that headers in fasta genome_main conform with standard

           https://github.com/LegumeFederation/datastore/issues/23
        '''
        logger = self.logger
        f_ids = self.fasta_ids
        true_header = '.'.join(attr[:3])
        fh = return_filehandle(fasta)  # get file handle, text/gz
        re_header = re.compile("^>(\S+)\s*(.*)")  # grab header and description
        passed = True
        with fh as gopen:
            for line in gopen:
                line = line.rstrip()
                if not line:
                    continue
                if re_header.match(line):  # check for fasta header
                    hid = re_header.search(line)
                    if hid:
                        logger.debug(hid.groups(0))
                        if isinstance(hid, basestring):  # check for tuple
                            hid = hid.groups(0)
                        else:
                            hid = hid.groups(0)[0]  # get id portion of header
                    else:
                        logger.error('Header {} looks odd...'.format(line))
                        sys.exit(1)
                    logger.debug(hid)
                    f_ids[hid] = 1
                    standard_header = true_header + '.' + hid
                    if not hid.startswith(true_header):
                        logger.warning(('Inconsistency {} '.format(hid) +
                                        'Should be {}'.format(standard_header))
                                      )
                        passed = False
        return passed

    def check_gff3_seqid(self, seqid):
        '''Confirms that column 1 "seqid" exists in genome_main if provided'''
        f_ids = self.fasta_ids  # fasta_ids generated from check_Reference
        if seqid not in f_ids:
            return False
        return True

    def check_seqid_attributes(self, gff):
        '''Confirms that gff3 seqid exists in genome_main if provided

           checks ID and Name from gff3 attributes field

           https://github.com/LegumeFederation/datastore/issues/23
        '''
        logger = self.logger
        fh = return_filehandle(gff)
        file_name = os.path.basename(gff)
        true_id = file_name.split('.')[:4]  # ID should start with this string
        true_name = file_name.split('.')[0]  # maybe this should include infra
        get_id_name = re.compile("^ID=(.+?);.*Name=(.+?);")
        lines = 0
        with fh as gopen:
            for line in gopen:
                line = line.rstrip()
                lines += 1
                if not line or line.startswith('#'):
                    continue
                columns = line.split('\t')  # get gff3 fields
                seqid = columns[0]  # seqid according to the spec
                seqid = seqid.rstrip()
                if self.fasta_ids:  # if genome_main make sure seqids exist
                    if not self.check_gff3_seqid(seqid):  # fasta header check
                        logger.debug(seqid)
                        logger.error('{} not found in genome_main'.format(
                                                                        seqid))
                feature_type = columns[3]  # get type
                attributes = columns[8]  # attributes ';' delimited
                if feature_type != 'gene':  # only check genes (for now)
                    continue
                if not get_id_name.match(attributes):  # check for ID and Name
                    logger.error('No ID and Name attributes. line {}'.format(
                                                                        lines))
                else:
                    groups = get_id_name.search(attributes).groups()
                    if len(groups) != 2:  # should just be ID and Name
                        logger.error('too many groups detected: {}'.format(
                                                                      groups))
                    (feature_id, feature_name) = groups
                    if not feature_id.startswith(true_id):  # check id
                        logger.error('feature id, should start with ' +
                                     '{} line {}'.format(true_id, lines))
                    if not feature_name.startswith(true_name):
                        logger.error('feature name, should start with ' +
                                     '{} line {}'.format(true_name, lines))

    def check_gff3(self, gff):
        '''Confirms that gff3 files pass gt validation

           https://github.com/LegumeFederation/datastore/issues/23
        '''
        logger = self.logger
        gt_path = self.gt_path
        logger.info('checking gff3 seqids and attributes...')
        self.check_seqid_attributes(gff)
        gff = os.path.basename(gff)
        gt_report = './{}_gt_gff3validator_report.txt'.format(gff)
        gt_cmd = '({}/gt gff3validator {} 2>&1) > {}'.format(gt_path, gff,
                                                             gt_report)
        logger.debug(gt_cmd)
        exit_val = subprocess.call(gt_cmd, shell=True)  # get gt exit_val
        logger.debug(exit_val)
        return exit_val

    def parse_filenames(self, f):
        '''parse attributes of filenames according to

           https://github.com/LegumeFederation/datastore/issues/23

           Will detect file type and error if it cannot
        '''
        #  maybe do a dir check thing here
        logger = self.logger
        logger.info('Checking filename standards')
        f_path = os.path.dirname(f)
        parent_dir = f_path.split('/')[-1]
        file_name = os.path.basename(f)
        file_attr = file_name.split('.')  # file name attributes split by "."
        if file_attr[-3] == 'genome_main':  # position of type
            logger.info('{} looks like a genome_main, processing...'.format(
                                                                   file_name))
            self.check_genome_main(file_attr)  # file naming is correct
            logger.info('Filename checks out.  Checking reference headers...')
            passed = self.check_fasta(f, file_attr)  # headers follow standard
            return passed
        if file_attr[-3] == 'gene_models_main':  # position of type
            logger.info('{} looks like gene_models_main, processing...'.format(
                                                                   file_name))
            self.check_gene_models_main(file_attr)  # check file naming
            logger.info('Filename checks out.  Checking GFF3 file...')
            exit_val = self.check_gff3(f)  # gff follows standard
            return exit_val

    def detect_incongruencies(self, **kwargs):
        '''Initiate and control workflow'''
        fh = ''
        logger = self.logger
        genome = self.genome
        annotation = self.annotation
        normalizer = self.normalizer
        if genome:
            logger.info('Genome will be checked...')
            genome = os.path.abspath(genome)
            if not check_file(genome):
                logger.error('Could not find {}'.format(genome))
                sys.exit(1)
            passed = self.parse_filenames(genome)
            if not passed and normalizer:
                logger.info('Normalizing {}'.format(genome))
                normalizer.normalize_genome_main(genome)
        if annotation:
            annotation = os.path.abspath(annotation)
            if not check_file(annotation):
                logger.error('Could not find {}'.format(annotation))
                sys.exit(1)
            exit_val = self.parse_filenames(annotation)  # gt exit value
            if exit_val:
                logger.warning('{} Failed gff3validator'.format(annotation))
            if normalizer and exit_val:  # gt said it wasn't clean, tidy
                logger.info('tiding gff3 file {}'.format(annotation))
                normalizer.tidy_gff3(annotation)
#         logger.info('Collecting report...') #  implement reporting at end
        logger.info('DONE')
