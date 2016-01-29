#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  aaf_phylokmer.py
#  
#  Copyright 2015,2016 Huan Fan
#  <hfan22@wisc.edu> & Yann Surget-Groba <yann@xtbg.org.cn>
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
import multiprocessing as mp
from optparse import OptionParser

def smartopen(filename,*args,**kwargs):
    if filename.endswith('gz'):
        return gzip.open(filename,*args,**kwargs)
    elif filename.endswith('bz2'):
        return bz2.BZ2File(filename,*args,**kwargs)
    else:
        return open(filename,*args,**kwargs)

def is_exe(fpath):
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

def runJob(command, sim):
    print command
    if not sim:
        os.system(command)
    return 1



usage = "usage: %prog [options]"
version = '%prog 20160104.1'
parser = OptionParser(usage = usage, version = version)
parser.add_option("-k", dest = "kLen", type = int, default = 25, 
                  help = "k-mer length, default = 25")
parser.add_option("-t", dest = "nThreads", type = int, default = 1, 
                  help = "number of threads to use, default = 1")
parser.add_option("-n", dest = "filter", type = int, default = 1,
                  help = "k-mer filtering threshold, default = 1")
parser.add_option("-o", dest = "outFile", default = 'phylokmer.dat.gz',
                  help = "output file, default = phylokmer.dat.gz")
parser.add_option("-d", dest = "dataDir", default = 'data',
                  help = "directory containing the data, default = data/")
parser.add_option("-G", dest = "memSize", type = int, default = 4,
                  help = "total memory limit (in GB), default = 4")
parser.add_option("-W", dest = "withKmer", action = 'store_true',
                  help = "include k-mers in the shared k-mer table, otherwise not")
parser.add_option("-s", dest = "sim", action = 'store_true',
                  help = "only print commands, do not run them")

(options, args) = parser.parse_args()

nThreads = options.nThreads
n = options.filter
memPerThread = int(options.memSize / float(nThreads))
if not memPerThread:
    print 'Not enough memory, decrease nThreads or increase memSize'
    sys.exit()
    

###check the data directory:
if not os.path.isdir(options.dataDir):
    print 'Cannot find data directory {}'.format(options.dataDir)
    sys.exit(2)


###check for the executable files:
#kmer_countx
if options.kLen > 25:
    if os.system('which kmer_countx > /dev/null'):
        kmerCount = './kmer_countx'
        if not is_exe(kmerCount):
            print 'kmer_countx not found. Make sure it is in your PATH or the'
            print 'current directory, and that it is executable'
            sys.exit(1)
    else:
        kmerCount = 'kmer_countx'

#kmer_count
else:
    if os.system('which kmer_count > /dev/null'):
        kmerCount = './kmer_count'
        if not is_exe(kmerCount):
            print 'kmer_count not found. Make sure it is in your PATH or the'
            print 'current directory, and that it is executable'
            sys.exit(1)
    else:
        kmerCount = 'kmer_count'

#kmer_merge
if os.system('which kmer_merge > /dev/null'):
    filt = './kmer_merge'
    if not is_exe(filt):
        print 'kmer_merge not found. Make sure it is in your PATH or the'
        print 'current directory, and that it is executable'
        sys.exit(1)
else:
    filt = 'kmer_merge'

if not options.sim:
    if os.path.exists(os.path.join(options.dataDir, options.outFile)):
        s = raw_input('{} is already in your data directory, overwrite it? Y/N '
                      .format(options.outFile))
        if s == 'Y' or s == 'y':
          print'{} is going to be overwritten'.format(options.outFile)
        else:
          print'No overwritting, exit'
          sys.exit(2)

###Get sample list:
samples = []
for fileName in os.listdir(options.dataDir):
    if os.path.isdir(os.path.join(options.dataDir, fileName)):
        samples.append(fileName)
    else:
        sample = fileName.split(".")[0]
        if sample in samples:
            sample = fileName.split(".")[0]+fileName.split(".")[1]
            if sample in samples:
                print 'Error, redundant sample or file names. Aborting!'
                sys.exit(3)
        os.system("mkdir {}/{}".format(options.dataDir,sample))
        os.system("mv {}/{} {}/{}/".format(options.dataDir,fileName,options.dataDir,sample))
        samples.append(sample)
samples.sort()

print 'SPECIES LIST:'
for sample in samples:
    print sample

###Prepare kmer_count jobs
jobList = []
for sample in samples:
	outFile = '{}.pkdat.gz'.format(sample)
	command = '{} -l {} -n {} -G {} -o {} -f '.format(kmerCount, options.kLen,
               n, memPerThread, outFile)
	command1 = ''
	for inputFile in os.listdir(os.path.join(options.dataDir, sample)):
        	inputFile = os.path.join(options.dataDir, sample, inputFile)
        	handle = smartopen(inputFile)
        	firstChar = handle.read(1)
        	if firstChar == '@':
            		seqFormat = 'FQ'
        	elif firstChar == '>':
            		seqFormat = 'FA'
        	else:
            		print 'Error, file {} is not FA or FQ format. Aborting!'.\
                   		format(inputFile)
            		sys.exit(3)
		command1 += " -i '{}'".format(inputFile)
    	command += '{}{}> {}.wc'.format(seqFormat,command1,sample)
    	jobList.append(command)
jobList = jobList[::-1]

###Run jobs
pool = mp.Pool(nThreads)
jobs = []
nJobs = 0
batch = 0
count = 0
nBatches = len(jobList) / nThreads
if len(jobList) % nThreads:
    nBatches += 1

while 1:
    if nJobs == nThreads:
        batch += 1
        print '\n', time.strftime('%c')
        print "running batch {}/{}".format(batch, nBatches)
        for job in jobs:
            pool.apply_async(runJob, args=[job, options.sim])
        pool.close()
        pool.join()
        pool = mp.Pool(nThreads)
        nJobs = 0
        jobs = []
    if jobList:
        command = jobList.pop()
        jobs.append(command)
        #job = pool.apply_async(runJob, args=[command, options.sim])
        nJobs += 1
    else:
        break 
    count += 1

if nJobs:
    print '\n', time.strftime('%c')
    print "running last batch"
    for job in jobs:
        pool.apply_async(runJob, args=[job, options.sim])
    pool.close()
    pool.join()

###Merge output wc files
if not options.sim:
    if options.outFile.endswith('.gz'):
        divFile = options.outFile.rstrip('.gz')+'.wc'
    else:
        divFile = options.outFile + '.wc'
    handle = open(divFile, 'w')
    handle.close()
    for sample in samples:
        kmerFile = sample + '.wc'
        os.system('cat {} >> {}'.format(kmerFile, divFile))
        os.remove(kmerFile)

###Run kmer_merge
if options.outFile.endswith('.gz'):
    outFile = options.outFile
else:
    outFile = options.outFile+'.gz'
if not options.sim:
    handle = smartopen(outFile, 'w')
    print >> handle, '#-k {}'.format(options.kLen)
    print >> handle, '#-n {}'.format(n)
    for i, sample in enumerate(samples):
        print >> handle, '#sample{}: {}'.format(i + 1, sample)
    handle.close()

command = "{} -k s -c -d '0' -a 'T,M,F'".format(filt)
cut = []
if options.withKmer:
	cut.append('1')
for i, sample in enumerate(samples):
    command += " '{}.pkdat.gz'".format(sample)
    cut.append(str((i + 1) * 2))

command += ' | cut -f {} | gzip >> {}'.format(','.join(cut), outFile)

print '\n', time.strftime('%c')
print command
if not options.sim:
    os.system(command)
print time.strftime('%c')
