#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  nonparametric_bootstrap_s2.py
#
#  Copyright 2015 Huan Fan <hfan22@wisc.edu>
#
#  This is a special nonparametric bootstrap script for AAF that starts with
#  phylokmer.dat.gz, meaning only doing the second stage of bootstrap. You can
#  set S1==0 in nonparametric_bootstrap.py but it will still generate the phylokmer.dat.gz
#  for one time.
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

import sys, gzip, bz2, os, time, random, math
import multiprocessing as mp
from optparse import OptionParser
#from Bio import SeqIO
from AAF import smartopen, is_exe, countTotal_shared, present
import numpy as np


usage = "usage: %prog [options]"
version = '%prog 20161204.1'
parser = OptionParser(usage = usage,version = version)
parser.add_option("-t", dest = "nThreads", type = int, default = 1,
                  help = "number of threads to use, default = 1")
#parser.add_option("-n", dest = "filter", type = int, default = 1,
#                  help = "kmer filtering threshold, default = 1")
parser.add_option("--fs", dest = "countfs", default = "phylokmer_singleton.wc",
                  help = "count of singletons, default = phylokmer_singleton.wc")
parser.add_option("-i", dest = "iptf", default = 'phylokmer.dat.gz',
                  help = "shared kmer table, default = phylokmer.dat.gz")
parser.add_option("-G", dest = "memSize", type = int, default = 4,
                  help = "total memory limit (in GB), default = 4")
parser.add_option("-l", dest = "long", action = 'store_true',
                  help = "use fitch_kmerX_long instead of fitch_kmerX")


(options, args) = parser.parse_args()

nThreads = options.nThreads
memory = options.memSize
###check for the executable files:

if os.system('which fitch_kmerX > /dev/null'):
    if options.long:
        fitch = './fitch_kmerX_long'
    else:
        fitch = './fitch_kmerX'
    if not is_exe(fitch):
        print(fitch+' not found. Make sure it is in your PATH or the')
        print('current directory, and that it is executable')
        sys.exit()
else:
    if options.long:
        fitch = 'fitch_kmerX_long'
    else:
        fitch = 'fitch_kmerX'

#check input files
try:
	kmerTable = smartopen(options.iptf)
except IOError:
	print('Cannot open file', options.iptf)
	sys.exit()

try:
	singleton = open(options.countfs)
except IOError:
	print('Cannot open file', options.countf)
	sys.exit()

###Read header and get sample list
samples = []              #species list
line = kmerTable.readline()
ll = line.split()
kl = int(ll[1])       #kmer length

while True:
    line = kmerTable.readline()
    if line.startswith('#-'):
        continue
    elif line.startswith('#sample'):
        ll = line.split()
        samples.append(ll[1])
    else:
        break
samples.sort()
sn = len(samples)    #species number
#get ntotal
Nsingle = [0.0] * sn

for i in range(sn):
    Nsingle[i] = np.random.poisson(int(singleton.readline().split()[1]))


nl = 0 #number of kmers in the shared kmer table
lines = kmerTable.readlines()
kmerTable.close()
simTable = open('simTable','wt')
table_len = len(lines)
for x in range(int(table_len/kl)):
    line_num = random.randint(0,table_len-1)
    simTable.write(lines[line_num])
simTable.close()
simTable = open('simTable','rt')
line = simTable.readline()
line_size = sys.getsizeof(line)
chunkLength = int(1024 ** 3 / nThreads / line_size)
print('chunkLength =', chunkLength)
intable = []
nshare = [[0] * sn for i in range(sn)]
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
            result = job.get()
            intable.append(result[0])
            shared = result[1]
            for i in range(sn):
                for j in range(i + 1, sn):
                    nshare[i][j] += shared[i][j]

        pool = mp.Pool(nThreads)
        nJobs = 0
        results = []
        print('{} running {} jobs'.format(time.strftime('%c'), nThreads))

    lines = []
    for i in range(chunkLength):
        if not line: #if empty
            break
        lines.append(line)
        line = simTable.readline()
    if not lines: #if empty
        break

    job = pool.apply_async(countTotal_shared, args=[lines,sn])
    results.append(job)
    nJobs += 1

if nJobs:
    print('{} running last {} jobs'.format(time.strftime('%c'), len(results)))
    pool.close()
    pool.join()
    for job in results:
        result = job.get()
        intable.append(result[0])
        shared = result[1]
        for i in range(sn):
            for j in range(i + 1, sn):
                nshare[i][j] += shared[i][j]

total_intable = np.sum(intable, axis = 0)
ntotal = np.sum ([total_intable,Nsingle],axis = 0)
dist = [[0] * sn for i in range(sn)]

for i in range(sn):
    for j in range(i + 1, sn):
        mintotal = min(ntotal[i], ntotal[j])
        if nshare[i][j] == 0:
            dist[j][i] = dist[i][j] = 1
        else:
            distance = (-1 / float(kl) * math.log(nshare[i][j] / float(mintotal)))
            dist[j][i] = dist[i][j] = distance
            nshare[j][i] = nshare[i][j]

###Write infile
try:
    infile = open('infile','wt')
except IOError:
    print('Cannot open infile for writing')
    sys.exit()
infile.write('{} {}'.format(sn, sn))
namedic = {}
for i in range(sn):
    lsl = len(samples[i])
    if lsl >= 10:
        ssl = samples[i][:10]
        appendix = 1
        while ssl in namedic:
            if appendix < 10:
                ssl = samples[i][:9]+str(appendix)
            elif appendix > 9:
                ssl = samples[i][:8]+str(appendix)
            appendix += 1
    else:
        ssl = samples[i] + ' ' * (10 - lsl)
    namedic[ssl] = samples[i]
    infile.write('\n{}'.format(ssl))
    for j in range(sn):
        infile.write('\t{}'.format(dist[i][j]))

infile.close()

###Run fitch_kmer
print('{} building tree'.format(time.strftime("%c")))
if os.path.exists("./outfile"):
    os.system("rm -f outfile outtree")
command = 'printf "K\n{}\nY" | {} > /dev/null'.format(int(kl),fitch)
os.system(command)
fh = open('outtree','rt')
fh1 = open(os.path.basename(options.iptf).split('.')[0]+'_bootstrap.tre','wt')

for line in fh:
    for key in namedic:
        key_new = key.rstrip()+":"
        if key_new in line:
            newline = line.replace(key_new,namedic[key].rstrip()+":",1)
            line = newline
    fh1.write(line) #This can be either line or new line because when it exits
#the for loop, line==newline
fh.close()
fh1.close()
