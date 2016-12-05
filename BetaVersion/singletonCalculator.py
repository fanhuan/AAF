#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  singletonCalculator.py
#
#  Copyright 2016 Huan Fan <hfan22@wisc.edu>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#

import sys,os,gzip,time
from optparse import OptionParser
import multiprocessing as mp
import numpy as np

def smartopen(filename,*args,**kwargs):
    '''opens with open unless file ends in .gz, then use gzip.open

        in theory should transparently allow reading of files regardless of compression
        '''
    if filename.endswith('.gz'):
        return gzip.open(filename,'rt',*args,**kwargs)
    else:
        return open(filename,'rt',*args,**kwargs)

def present(x,n=1):
    if int(x) >= n:
        return '1'
    else:
        return '0'

def rc(seq):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U':'A'}
	return "".join(complement.get(base, base) for base in reversed(seq))

def countTotal(lines):
    line_list = []
    for line in lines:
        line_list.append([int(present(i)) for i in line.split()[1:]])
    line_total = np.sum(line_list,axis = 0)
    return line_total


Usage = "%prog [options] shared_kmer_table kmer_diversity_file k"
version = '%prog 20161204.1'
parser = OptionParser(Usage, version = version)

parser.add_option("-t", dest = "nThreads", type = int, default = 1,
                  help = "number of threads to use, default = 1")

(options, args) = parser.parse_args()
kmer_table = smartopen(sys.argv[1])
total = smartopen(sys.argv[2])
k = int(sys.argv[3])
nThreads = options.nThreads
output = open(os.path.basename(sys.argv[1]).split('.')[0]+'_singleton.wc','w')

# Read in the kmer diversity file
ntotal = []
for line in total:
    ntotal.append(float(line.split()[1]))

###Read header
sl = []                 #species list
for line in kmer_table:
    if line.startswith('#-'):
        if line.startswith('#-k'):
            k_table = int(line.split()[1])
            if k_table != k:
                print('the k givin is different from the one in the kmer table, abort')
                sys.exit()
    elif line.startswith('#sample'):
        ll = line.split()
        sl.append(ll[1])
    else: #if there is no header
        line_size = sys.getsizeof(line)
        line1 = line.split()
        if line1[0].startswith(tuple('ATCG')):
        	sn = len(line1) - 1
        elif line1[0].startswith(tuple('0123456789')):
            sn = len(line1)
        break

if len(sl) != sn:
    print('not enough columns in the the kmer_table')
    sys.exit()
if len(ntotal) != sn:
    print('the number of samples in the diversity file does not match the kmer table, abort')
    sys.exit()
###Compute the number of lines to process per thread

chunkLength = int(1024 ** 3 / nThreads / line_size)
print('chunkLength =', chunkLength)
intable = []
###Compute pattern dictionary
nJobs = 0
pool = mp.Pool(nThreads)
results = []
print(time.strftime('%c'), 'start running jobs')

while True:
    if nJobs == nThreads:
        pool.close()
        pool.join()
        for job in results:
            intable_part = []
            intable_part = job.get()
            intable.append(intable_part)
        pool = mp.Pool(nThreads)
        nJobs = 0
        results = []
        print('{} running {} jobs'.format(time.strftime('%c'), nThreads))

    lines = []
    for nLines in range(chunkLength):
        if not line: #if empty
            break
        lines.append(line)
        line = kmer_table.readline()
    if not lines: #if empty
        break

    job = pool.apply_async(countTotal, args=[lines])
    results.append(job)
    nJobs += 1

if nJobs:
    print('{} running last {} jobs'.format(time.strftime('%c'), len(results)))
    pool.close()
    pool.join()
    i = 0
    for job in results:
        i = i + 1
        intable_part = []
        intable_part = job.get()
        intable.append(intable_part)

kmer_table.close()
total_intable = np.sum(intable, axis = 0)
singleton = np.sum([ntotal, [-a for a in total_intable]],axis = 0)
for i in range(sn):
    output.write('{}: {}\n'.format(sl[i],int(singleton[i]/k)))
output.close()
