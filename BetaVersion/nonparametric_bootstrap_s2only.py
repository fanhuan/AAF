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
from Bio import SeqIO

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

def countLine(lines, sn): #count both shared the total kmers, for total kmer table
    shared = [[0] * sn for i in xrange(sn)]
    for line in lines:
        line = line.split()
	if len(line) == sn+1:
		line = line[1:]
        line = [int(i) for i in line]
        for i in xrange(sn):
            for j in xrange(i + 1, sn):
                if line[i]*line[j] != 0:
                    shared[i][j] += 1
    return shared

def aaf_distance(outFile,t,m,samples,kl,s2,ntotal):
    try:
        iptf = smartopen(outFile)
    except IOError:
        print 'Cannot open file', outFile
        sys.exit()
    try:
        infile = open('infile','w')
    except IOError:
        print 'Cannot open infile for writing'
        sys.exit()
    ###Initialize shared kmers matrix
    nshare = [[0] * sn for i in xrange(sn)]
    ###Compute the number of lines to process per thread
    line = iptf.readline()
    line_size = sys.getsizeof(line)
    chunkLength = int(memory * 1024 ** 3 / nThreads / line_size)
    
    ###Compute shared kmer matrix
    nJobs = 0
    pool = mp.Pool(nThreads)
    results = []
    print '\n', time.strftime('%c'), 'start calulation distances'
    while True:
        if nJobs == nThreads:
            #print '{} running {} jobs'.format(time.strftime('%c'), nThreads)
            print '.',
            pool.close()
            pool.join()
            for job in results:
                shared = job.get()
                for i in xrange(sn):
                    for j in xrange(i + 1, sn):
                        nshare[i][j] += shared[i][j]
            
            pool = mp.Pool(nThreads)
            nJobs = 0
            results = []
        
        lines = []
        for nLines in xrange(chunkLength):
            if not line: #if empty
                break
            lines.append(line)
            line = iptf.readline()
        if not lines: #if empty
            break
        job = pool.apply_async(countLine, args=[lines, sn])
        results.append(job)
        nJobs += 1
    
    if nJobs:
        print '\nrunning last jobs'
        pool.close()
        pool.join()
        for job in results:
            shared = job.get()
            for i in xrange(sn):
                for j in xrange(i + 1, sn):
                    nshare[i][j] += shared[i][j]
    
    iptf.close()
    
    ###Compute distance matrix
    dist = [[0] * sn for i in xrange(sn)]
    for i in xrange(sn):
        for j in xrange(i + 1, sn):
            mintotal = float(min(ntotal[i], ntotal[j]))
            if nshare[i][j] == 0:
                dist[j][i] = dist[i][j] = 1
            else:
                distance = (-1.0 / float(kl)) * math.log(nshare[i][j] / mintotal)
                dist[j][i] = dist[i][j] = distance
                nshare[j][i] = nshare [i][j]

    ###Write infile
    infile.write('{} {}'.format(sn, sn))
    namedic = {}
    for i in xrange(sn):
        lsl = len(samples[i])
        if lsl >= 10:
            ssl = samples[i][:10]
            appendix = 1
            while ssl in namedic:
                ssl = samples[i][:9]+str(appendix)
                appendix += 1
        else:
            ssl = samples[i] + ' ' * (10 - lsl)
        namedic[ssl] = samples[i]
        infile.write('\n{}'.format(ssl))
        for j in xrange(sn):
            infile.write('\t{}'.format(dist[i][j]))
    infile.close()
    
    ###Run fitch_kmer
    print time.strftime("%c"), 'building tree'
    if os.path.exists("./outfile"):
        os.system("rm -f outfile outtree")
    command = 'printf "K\n{}\nY" | {} > /dev/null'.format(kl,fitch)
    print command
    os.system(command)
    print time.strftime("%c"), 'end, chunk size =', chunkLength
    #command = 'mv infile infile_{}'.format(s2)
    #os.system(command)
    return

usage = "usage: %prog [options]"
version = '%prog 20151105.1'
parser = OptionParser(usage = usage,version = version)
parser.add_option("-t", dest = "nThreads", type = int, default = 1,
                  help = "number of threads to use, default = 1")
parser.add_option("-n", dest = "filter", type = int, default = 1,
                  help = "kmer filtering threshold, default = 1")
parser.add_option("-f", dest = "countf", default = "kmer_diversity.wc",
                  help = "k-mer diversity file, default = kmer_diversity.wc")
parser.add_option("-i", dest = "iptf", default = 'phylokmer.dat.gz',
                  help = "shared kmer table, default = phylokmer.dat.gz")
parser.add_option("-G", dest = "memSize", type = int, default = 4,
                  help = "total memory limit (in GB), default = 4")
parser.add_option("--S2",dest = "stage2", type = int, default = 0,
                  help = "number of resampling of each kmer table, default = 0")
parser.add_option("-s", dest = "sim", action = 'store_true',
                  help = "only print commands, do not run them")

(options, args) = parser.parse_args()

nThreads = options.nThreads
memory = options.memSize

###check for the executable files:

#fitch_kmerX
if os.system('which fitch_kmerX > /dev/null'):
    fitch = './fitch_kmerX'
    if not is_exe(fitch):
        print 'fitch_kmerX not found. Make sure it is in your PATH or the'
        print 'current directory, and that it is executable'
        sys.exit()
else:
    fitch = 'fitch_kmerX'

#consense
if os.system('which consense > /dev/null'):
    consense = './consense'
    if not is_exe(fitch):
        print 'consense not found. Make sure it is in your PATH or the'
        print 'current directory, and that it is executable'
        sys.exit()
else:
    consense = 'consense'

#check input files
try:
	kmerTable = smartopen(options.iptf)
except IOError:
	print 'Cannot open file', options.iptf
	sys.exit()

if not os.path.isfile(options.countf):
	print 'Cannot find file', options.countf
	sys.exit()

try:
	total = open(options.countf)
except IOError:
	print 'Cannot open file', options.countf
	sys.exit()



###Set filter option
n = options.filter

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
Ntotal = [0.0] * sn

for i in xrange(sn):
    Ntotal[i] = float(total.readline().split()[1])

print 'SPECIES LIST:'
for sample in samples:
    print sample
if os.path.exists("./consensus_trees_table_nonparametric"):
    os.system("rm -f consensus_trees_table_nonparametric consensus_tree_table_nonparametric.tre consense_outfile_table_nonparametric")

lines= kmerTable.readlines()
table_len = len(lines)
for s2 in xrange(options.stage2):
    print "{} out of {} times of bootstrap over the real total kmer table.".format(s2+1,options.stage2)
    simTable = open('simTable_{}'.format(s2),'w')
#    for x in xrange(table_len):
    for x in xrange(table_len/kl):
        line_num = random.randint(0,table_len-1)
        simTable.write(lines[line_num])
    simTable.close()
    aaf_distance('simTable_{}'.format(s2),nThreads,memory,samples,kl,s2,Ntotal)
    os.system('cat outtree >> consensus_trees_table_nonparametric')
kmerTable.close()


#consense
if os.path.exists("./outfile"):
    os.system("rm -f outfile outtree")

command = 'printf "consensus_trees_table_nonparametric\nY" | {} > /dev/null'.format(consense)
print command
os.system(command)
os.system("mv outfile consense_outfile_table_nonparametric")
os.system("mv outtree consensus_tree_table_nonparametric.tre")

os.system("rm simTable_*")

