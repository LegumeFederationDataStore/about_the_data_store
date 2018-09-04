#!/usr/bin/env python

import os, sys
import logging
import gzip


def return_filehandle(open_me):
    '''return file handle for gz compressed or text file'''
    magic_dict = { # headers for compression
                  '\x1f\x8b\x08': 'gz',
 #                 '\x42\x5a\x68': 'bz2',
 #                 '\x50\x4b\x03\x04': 'zip'
                 }
    max_bytes = max(len(t) for t in magic_dict)
    with open(open_me) as f:
        s = f.read(max_bytes)
    for m in magic_dict:
        if s.startswith(m): #check file header for match with m
            t = magic_dict[m]
            if t == 'gz':
                return gzip.open(open_me)
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



if __name__ == '__main__':
    print('import me to use check_file and return_filehandle.' + 
          'This should be in a class probably')
    sys.exit(1)
