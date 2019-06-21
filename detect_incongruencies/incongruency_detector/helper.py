#!/usr/bin/env python

import os, sys
import gzip
import errno
import logging
import subprocess


def check_subprocess_dependencies():
    '''attempts to run samtools, tabix and bgzip'''
    subprocess.check_call('which bgzip', shell=True)  # check in env
    subprocess.check_call('which tabix', shell=True)  # check in env
    subprocess.check_call('which samtools', shell=True)  # check in env
    return True


def setup_logging(log_file, log_level):
    '''Return logger based on log_file and log_level'''
    log_level = getattr(logging, log_level.upper(), logging.INFO)
    msg_format = '%(asctime)s|%(name)s|[%(levelname)s]: %(message)s'
    logging.basicConfig(format=msg_format, datefmt='%m-%d %H:%M',
                        level=log_level)
    log_handler = logging.FileHandler(log_file, mode='w')
    formatter = logging.Formatter(msg_format)
    log_handler.setFormatter(formatter)
    logger = logging.getLogger('detect_incongruencies')
    logger.addHandler(log_handler)
    return logger


def return_filehandle(open_me):
    '''return file handle for gz compressed or text file'''
    magic_dict = { # headers for compression
                  b'\x1f\x8b\x08': 'gz',
 #                 '\x42\x5a\x68': 'bz2',
 #                 '\x50\x4b\x03\x04': 'zip'
                 }
    max_bytes = max(len(t) for t in magic_dict)
    with open(open_me, 'rb') as f:
        s = f.read(max_bytes)
    for m in magic_dict:
        if s.startswith(m): #check file header for match with m
            t = magic_dict[m]
            if t == 'gz':
                return gzip.open(open_me, 'rt')
  #          elif t == 'bz2':
  #              return bz2.open(open_me)
  #          elif t == 'zip':
  #              return zipfile.open(open_me)
    return open(open_me)


def check_file(f):
    '''check for file using os.path.isfile'''
    try:
        os.path.isfile(f)
    except OSError:
        raise
    return os.path.isfile(f)


def create_directories(dirpath):
    '''make directory path'''
    try:
        os.makedirs(dirpath)
    except OSError as e:
        if e.errno != errno.EEXIST:  # ignore if error is exists else raise
            raise


if __name__ == '__main__':
    print('import me to use check_file and return_filehandle.' + 
          'This should be in a class probably')
    sys.exit(1)
