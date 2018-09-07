#!/usr/bin/env python

import os
import sys
import logging
import re
import subprocess
import hashlib
import requests
import json
from glob import glob
from os import path
sys.path.append( path.dirname( path.dirname( path.abspath(__file__) ) ) )
from .helper import check_file, return_filehandle

class genome_main:

    def __init__(self, detector, **kwargs):
        self.detector = detector
        self.target = detector.target
        self.logger = detector.logger
        self.fasta_ids = {}
        detector.fasta_ids = self.fasta_ids

    def run(self):
        '''Runs checks'''
        logger = self.logger
        target = self.target
        logger.info('Performing Naming Checks for {}\n'.format(target))
        if not self.check_genome_main():
            logger.error('Naming Check FAILED')
            return False
        logger.info('Naming Looks Correct\n')
        logger.info('Checking Genome FASTA\n')
        if not self.check_genome_fasta():
            logger.error('FASTA Check FAILED')
            return False
        logger.info('ALL GENOME CHECKS PASSED\n')
        return True

    def check_genome_main(self):
        '''accepts a list of genome attributes split by "."

           https://github.com/LegumeFederation/datastore/issues/23

           checks these file attributes to ensure they are correct
        '''
        target = self.target
        logger = self.logger
        attr = os.path.basename(target).split('.')  # split on delimiter
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
        logger.info('Genome Naming Looks Correct')
        return True

    def check_genome_fasta(self):
        '''Confirms that headers in fasta genome_main conform with standard

           PUT SOME RULE REFERENCE HERE
        '''
        logger = self.logger
        fasta = self.target  # get fasta file
        attr = os.path.basename(fasta).split('.')  # get attributes for naming
        self.fasta_ids = {}  # initialize fasta ids for self
        f_ids = self.fasta_ids  # set to overwrite for each reference
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
                        if isinstance(hid, str):  # check for tuple
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


class gene_models_main:

    def __init__(self, detector, **kwargs):
        self.detector = detector
        self.target = detector.target
        self.logger = detector.logger
        self.fasta_ids = detector.fasta_ids

    def run(self):
        '''Run checks'''
        logger = self.logger
        target = self.target
        logger.info('Checking Gene Models File Naming...')
        if not self.check_gene_models_main():
            logger.error('Gene Model File Naming FAILED')
            return False
        logger.info('Naming Looks Correct\n')
        logger.info('Validating GFF3 with gt...')
        if not self.gt_gff3_validate():
            logger.error('gt GFF3 Validation FAILED')
            return False
        logger.info('GFF3 is Valid\n')
        logger.info('Checking Congruency Between Genome and Gene Models...')
        if not self.check_seqid_attributes():
            logger.error('Genome and Gene Models are not Congruent FAILED')
            return False
        logger.info('ALL GENE MODELS CHECKS PASSED\n')
        return True

    def check_gene_models_main(self):
        '''accepts a list of annotation attributes split by "."

           https://github.com/LegumeFederation/datastore/issues/23

           checks these file attributes to ensure they are correct
        '''
        target = self.target
        logger = self.logger
        attr = os.path.basename(target).split('.')  # split on delimiter 
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
        return True

    def gt_gff3_validate(self):
        '''Confirms that gff3 files pass gt validation

           https://github.com/LegumeFederation/datastore/issues/23
        '''
        gff = self.target
        logger = self.logger
        gff_name = os.path.basename(gff)
        gt_report = './{}_gt_gff3validator_report.txt'.format(gff_name)
        gt_cmd = '(gt gff3validator {} 2>&1) > {}'.format(gff,
                                                         gt_report)
        logger.debug(gt_cmd)
        exit_val = subprocess.call(gt_cmd, shell=True)  # get gt exit_val
        logger.debug(exit_val)
        if exit_val:
            return False
        return True
    
    def check_seqid_attributes(self):
        '''Confirms that gff3 seqid exists in genome_main if provided

           checks ID and Name from gff3 attributes field

           https://github.com/LegumeFederation/datastore/issues/23
        '''
        gff = self.target
        logger = self.logger
        fasta_ids = self.fasta_ids  # list of FASTA IDS from Reference
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
                logger.debug(line)
                logger.debug(seqid)
                if self.fasta_ids:  # if genome_main make sure seqids exist
                    if not seqid in fasta_ids:  # fasta header check
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
        return True
