#!/usr/bin/env python

import os
import sys
import logging
import re
import subprocess
import gzip
from file_helpers import check_file, return_filehandle


class Normalizer:
    '''Class to normalize datastore file incongruencies with

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
        if self.annotation and not self.gt_path:
            logger.error('gt_path, /path/to/genome_tools must be provided')
            sys.exit(1)
        if self.gt_path:
            self.gt_path = os.path.abspath(self.gt_path)
#        self.fasta_ids = {}
#        self.genome_attributes = {'filename': '', 'version': '',
#                                  'prefix': '', 'type': '', 'build': '',
#                                  'compression': ''}
#        self.incongruencies = {'Genome Name': [], 'Genome Headers': [],
#                               'Annotation Name': [], 'GFF File': []}
        self.logger.info('Initialized Normalizer')

    def normalize_genome_main(self, genome):
        '''accepts the genome prefix and a key from the legfed_registry

           https://github.com/LegumeFederation/datastore/issues/23
        '''
        logger = self.logger
        file_name = os.path.basename(genome)  # get filename 
        prefix = '.'.join(file_name.split('.')[:3])  # get correct prefix
        logger.debug(prefix)
        fh = return_filehandle(genome)
        re_header = re.compile("^>(\S+)\s*(.*)")  # match header
        fasta_out = gzip.open('./{}.normalized'.format(file_name), 'wb')
        with fh as gopen:
            for line in gopen:
                line = line.rstrip()
                if re_header.match(line):
                    parsed_header = re_header.search(line).groups()
                    logger.debug(line)
                    logger.debug(parsed_header)
                    hid = ''
                    desc = ''
                    new_header = ''
                    if isinstance(parsed_header, basestring):  # check tuple
                        hid = parsed_header
                        new_header = '>{}.{}'.format(prefix, hid)
                        logger.debug(hid)
                    else:
                        if len(parsed_header) == 2:  # header and description
                            hid = parsed_header[0]  # header
                            desc = parsed_header[1]  # description
                            new_header = '>{}.{} {}'.format(prefix, hid, desc)
#                    normalized_header = '>'
                            logger.debug(hid)
                            logger.debug(desc)
                            logger.debug(new_header)
                        else:
                            logger.error('header {} looks odd...'.format(line))
                            sys.exit(1)
                    fasta_out.write(new_header + '\n')  # write new header
                else:  # sequence lines
                    fasta_out.write(line + '\n')
        fasta_out.close()

    def tidy_gff3(self, gff):
        '''Performs gt gff -sort -tidy as parameterized by:

           https://github.com/genometools/genometools/wiki/speck-User-manual
        '''
        logger = self.logger
        gt_path = self.gt_path  # path to gt tool
        gff = os.path.basename(gff)
        gt_report = './{}_gt_gff3_tidy_report.txt'.format(gff)  # gt outfile
        tidy_out = './{}.tidy'.format(gff)  # tidied gff
        gt_cmd = ('{}/gt gff3 -sort -tidy -retainids -force '.format(gt_path) +
                  '-o {} -gzip {} 2> {}'.format(tidy_out, gff, gt_report))
        logger.debug(gt_cmd)
        exit_val = subprocess.call(gt_cmd, shell=True)
        if exit_val:
            logger.error('Exit value of gt gff3 tidy !=0: {}'.format(exit_val))
            sys.exit(1)
        tidy_out = tidy_out + '.gz'
        exit_val = self.check_gff3(tidy_out)
        if exit_val:
            logger.error('Tidy failed validation!')
            sys.exit(1)
        logger.debug(exit_val)
