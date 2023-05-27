#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  aaf_phylogeny.py
#
#  Copyright 2021 Huan Fan <huan.fan@wisc.edu>
#
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

import sys, gzip, bz2, os, time
from concurrent.futures import ProcessPoolExecutor as PPE
from optparse import OptionParser
from AAF import smartopen, is_exe, countShared_single, aaf_dist, run_command, aaf_kmercount

usage = "usage: %prog [options]"
version = '%prog 20220817.1'
parser = OptionParser(usage = usage, version = version)
parser.add_option("-k", dest = "kLen", type = int, default = 25,
                  help = "k for reconstruction, default = 25")
parser.add_option("-n", dest = "filter", type = int, default = 1,
                  help = "k-mer filtering threshold, default = 1")
parser.add_option("-d", dest = "dataDir", default = 'data',
                  help = "directory containing the data, default = data/")
parser.add_option("-G", dest = "memSize", type = int, default = 4,
                  help = "total memory limit (in GB), default = 4")
parser.add_option("-t", dest = "nThreads", type = int, default = 1,
                  help = "number of threads to use, default = 1")
parser.add_option("-l", dest = "long", action = 'store_true',
                  help = "use fitch_kmerX_long instead of fitch_kmerX")

(options, args) = parser.parse_args()

nThreads = options.nThreads
n = options.filter
memPerThread = int(options.memSize / float(nThreads))
kl = options.kLen
dataDir = options.dataDir
memSize = options.memSize

if not memPerThread:
    print('Not enough memory, decrease nThreads or increase memSize')
    sys.exit()


###check the data directory:
if not os.path.isdir(dataDir):
    print('Cannot find data directory {}'.format(dataDir))
    sys.exit(2)


###check for the executable files:
#kmer_merge
if os.system('which kmer_merge > /dev/null'):
    filt = './kmer_merge'
    if not is_exe(filt):
        print('kmer_merge not found. Make sure it is in your PATH or the')
        print('current directory, and that it is executable')
        sys.exit(1)
else:
    filt = 'kmer_merge'


###Run aaf_kmercount to get pkdat for each species
samples = aaf_kmercount(dataDir,kl,n,nThreads,memSize/nThreads)
###Merge output wc files
divFile = dataDir+'.wc'
handle = open(divFile, 'w')
handle.close()

for sample in samples:
    countfile = sample + '.wc'
    os.system('cat {} >> {}'.format(countfile, divFile))
    os.remove(countfile)

###Run kmer_merge
outFile = dataDir+'.dat.gz'
handle = smartopen(outFile, 'w')
handle.write(('#-k {}\n'.format(kl)).encode('latin-1'))
handle.write(('#-n {}\n'.format(n)).encode('latin-1'))
for i, sample in enumerate(samples):
    handle.write(('#sample{}: {}\n'.format(i + 1, sample)).encode('latin-1'))
handle.close()

command = "{} -k s -c -d '0' -a 'T,M,F'".format(filt)
cut = []
for i, sample in enumerate(samples):
    command += " '{}.pkdat.gz'".format(sample)
    cut.append(str((i + 1) * 2))

command += ' | cut -f {} | gzip >> {}'.format(','.join(cut), outFile)
print('\n', time.strftime('%c'))
print(command)
os.system(command)
print(time.strftime('%c'))

#Calculate the distance and generate the tree!
aaf_dist(outFile,divFile,nThreads,samples,kl)
command = 'mv aaf.tre {}.tre'.format(dataDir)
# command = 'mv aaf.tre %s_k%s.tre'
#'cat %s %s > %s' % (fq1,fq2,Args.workdir + output_file +'.fq')
command = 'mv aaf.dist {}.dist'.format(dataDir)
