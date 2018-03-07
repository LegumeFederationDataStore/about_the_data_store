#!/usr/bin/env python

import os
import sys
import logging
import re
import subprocess
import hashlib
from glob import glob
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
        self.directory = kwargs.get('directory')
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
        f_ids = {}
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
        gff_name = os.path.basename(gff)
        gt_report = './{}_gt_gff3validator_report.txt'.format(gff_name)
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

    def parse_checksum(self, md5_file, check_me):
        '''Get md5 checksum for file and compare to expected'''
        logger = self.logger
        fh = return_filehandle(md5_file)
        hash_md5 = hashlib.md5()
        check_sum_target = ''
        switch = 0
        with fh as copen:
            for line in copen:
                line = line.rstrip()
                if not line or line.startswith('#'):
                    continue
                fields = line.split(' ')
                check_sum = fields[0]
                filename = fields[1]
                if not check_sum and filename:
                    logger.error('Could not find sum and name for {}'.format(
                                                                        line))
                if filename == os.path.basename(check_me):
                    logger.info('Checksum found for {}'.format(filename))
                    check_sum_target = check_sum
                    switch = 1
        if not switch:
            logger.error('Could not find checksum for {}'.format(check_me))
            sys.exit(1)
        with open(check_me, "rb") as copen:
            for chunk in iter(lambda: copen.read(4096), b""):  # 4096 buffer
                hash_md5.update(chunk)
        target_sum = hash_md5.hexdigest()  # get sum
        logger.debug(target_sum)
        logger.debug(check_sum_target)
        if target_sum != check_sum_target:  # compare sums
            logger.error(('Checksum for file {} {} '.format(check_me, 
                                                           target_sum) + 
                          'did not match {}'.format(check_sum_target)))
            sys.exit(1)
        logger.info('Checksums checked out, moving on...')

    def check_dir_type(self, directory):
        '''Check the directory type and perform workflow based on type'''
        main_file = ''
        file_type = ''
        logger = self.logger
        dir_name = os.path.basename(directory)  # get dirname only
        dir_type = dir_name.split('.')[-2]  # get gnm, ann, etc
        glob_str = '{}/*'.format(directory)
        if dir_type.startswith('gnm'):
            logger.info('This is a genome directory.  Looking for genome_main')
            glob_str += 'genome_main.fna.gz'  # glob main fasta
            fasta = glob(glob_str)
            if len(fasta) != 1:
                logger.warning('Multiple/0 genome_main found {}'.format(fasta))
                return False
            main_file = fasta[0]
            logger.info('Found {}'.format(main_file))
            file_type = 'genome'  # set type for return
        elif dir_type.startswith('ann'):
            logger.info(('This is an annotation directory. ' +  
                         'Looking for gene_models_main'))
            glob_str += 'gene_models_main.gff3.gz'  # glob main gff
            gff = glob(glob_str)
            if len(gff) != 1:
                logger.warning('Multiple/0 gene_models_main for {}'.format(gff))
                return False
            main_file = gff[0]
            logger.info('Found {}'.format(main_file))
            file_type = 'annotation'  # set type for return
        else:
            logger.warning(('Format {} not recognized, '.format(dir_type) +
                          'should be ann or gnm.  Skipping...'))
            return False
        logger.info('Searching for checksum...')
        check_glob = '{}/CHECKSUM.*.md5'.format(directory)
        check_sum = glob(check_glob)
        if len(check_sum) != 1:
            logger.warning('Multiple/0 checksums for {}'.format(main_file))
            return False
        check_sum_file = check_sum[0]
        self.parse_checksum(check_sum_file, main_file)
        return (main_file, file_type)

    def run_genome(self, genome):
        '''Run genome workflow'''
        logger = self.logger
        normalizer = self.normalizer
        logger.info('Genome {} will be checked...'.format(genome))
        genome = os.path.abspath(genome)
        if not check_file(genome):
            logger.error('Could not find {}'.format(genome))
            sys.exit(1)
        passed = self.parse_filenames(genome)
        if not passed and normalizer:
            logger.info('Normalizing {}'.format(genome))
            normalizer.normalize_genome_main(genome)

    def run_annotation(self, annotation):
        '''Run annotation workflow'''
        logger = self.logger
        normalizer = self.normalizer
        annotation = os.path.abspath(annotation)
        logger.info('Annotation {} will be checked...'.format(annotation))
        if not check_file(annotation):
            logger.error('Could not find {}'.format(annotation))
            sys.exit(1)
        exit_val = self.parse_filenames(annotation)  # gt exit value
        if exit_val:
            logger.warning('{} Failed gff3validator'.format(annotation))
        if normalizer and exit_val:  # gt said it wasn't clean, tidy
            logger.info('tiding gff3 file {}'.format(annotation))
            normalizer.tidy_gff3(annotation)

    def get_files(self, directory):
        '''Get all related files, start with gnm return list of dicts
        
           for processing
        '''
        logger = self.logger
        gnm_glob = '{}/*/*.gnm[0-9]*'.format(directory)
        gnm_list = glob(gnm_glob)
        related_files = []
        logger.info('Globbing files to find main elements...')
        for g in gnm_list:
            if g.endswith('genome_main.fna.gz'):  # add more checks, transcriptome etc
#            elif g.endswith('transcriptome_main.fna.gz')  # others
                logger.info('Found genome {}'.format(g))
                filename = os.path.basename(g)
                prefix = '.'.join(filename.split('.')[:3])
                ann_glob = '{}/*/{}.ann[0-9]*gene_models_main.gff3.gz'.format(
                                                                    directory,
                                                                    prefix)
                anns = glob(ann_glob)
                if not anns:
                    logger.debug('No genemodels main for {}'.format(filename))
                    anns = None
                elif len(anns) > 1:
        #            logger.warning('Multiple gene models fround for {}'.format(
        #                                                             filename))
        #            continue
                    t_anns = []
                    for a in anns:
                        t_anns.append(os.path.dirname(a))
                    anns = t_anns
                else:
                    logger.info('Found gene_models {}'.format(anns[0]))
                    anns = [os.path.dirname(anns[0])]
                files_obj = {'genome': os.path.dirname(g), 
                             'annotation': anns}
                related_files.append(files_obj)
            else:
                logger.debug('File {} not cannonical will not check'.format(g))
                continue
        return related_files

    def detect_incongruencies(self, **kwargs):
        '''Initiate and control workflow'''
        fh = ''
        logger = self.logger
        genome = self.genome
        annotation = self.annotation
        directory = self.directory
        if directory:
            directory = os.path.abspath(directory)
            directories = self.get_files(directory)  # get object for loop
            logger.debug(directories)  # see what is going into loop
            for d in directories:
                genome = d.get('genome')  # get and check
                annotation = d.get('annotation') # get and check
                logger.info('Checking Genome:{} and Annotation:{}'.format(
                                                                 genome,
                                                                 annotation))
                if genome:
                    main_file, file_type = self.check_dir_type(genome)
                    if file_type == 'genome':
                        self.run_genome(main_file)
                    else:
                        logger.warning('Assembly does not look like a genome')
                        continue
                if annotation:
                    for a in annotation:
                        main_file, file_type = self.check_dir_type(a)
                        if file_type == 'annotation':
                            self.run_annotation(main_file)
                        else:
                            logger.warning('Annotation looks odd...')
                            continue
#                else:
#                    logger.warning('Did not recognize type {}'.format(
#                                                                 file_type))
#                    continue
                logger.info('Done Checking, Proceeding to next target...')
            return 
        if genome:
            self.run_genome(genome)
        if annotation:
            self.run_annotation(annotation)
            
#         logger.info('Collecting report...') #  implement reporting at end
        logger.info('DONE')
