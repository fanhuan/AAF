#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  nonparametric_bootstrap.py
#
#  Copyright 2014 Huan Fan <hfan22@wisc.edu>
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
    total = [0.0] * sn
    for line in lines:
        line = line.split()
        line = [int(i) for i in line]
        for i in xrange(sn):
            if line[i] > 0:
                total[i] += 1
            for j in xrange(i + 1, sn):
                if line[i]*line[j] != 0:
                    shared[i][j] += 1
    return total,shared

def aaf_distance(outFile,t,m,samples,kl):
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
    sn = len(samples)    #species number
    nshare = [[0] * sn for i in xrange(sn)]
    ntotal = [0.0] * sn
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
                total,shared = job.get()
                for i in xrange(sn):
                    ntotal[i] += total[i]
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
            total,shared = job.get()
            for i in xrange(sn):
                ntotal[i] += total[i]
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
                distance = (-1.0 / kl) * math.log(nshare[i][j] / mintotal)
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
    command = 'printf "K\n{}\nY" | {} > /dev/null'.format(int(kl),fitch)
    print command
    os.system(command)
    print time.strftime("%c"), 'end, chunk size =', chunkLength
    return

usage = "usage: %prog [options]"
version = '%prog 20141013.1'
parser = OptionParser(usage = usage,version = version)
parser.add_option("-k", dest = "kLen", type = int, default = 25,
                  help = "kmer length, default = 25")
parser.add_option("-t", dest = "nThreads", type = int, default = 1,
                  help = "number of threads to use, default = 1")
parser.add_option("-n", dest = "filter", type = int, default = 1,
                  help = "kmer filtering threshold, default = 1")
parser.add_option("-f", dest = "seqFormat", default = "FA",
                  help = "format of input files, FA|FQ, default = FA")
parser.add_option("-o", dest = "outFile", default = 'phylokmer.dat.gz',
                  help = "k-mer table name, default = phylokmer.dat.gz")
parser.add_option("-d", dest = "dataDir", default = 'data',
                  help = "directory containing the data, default = data/")
parser.add_option("-G", dest = "memSize", type = int, default = 4,
                  help = "total memory limit (in GB), default = 4")
parser.add_option("--S1",dest = "stage1", type = int, default = 0,
                  help = "number of resampling of the reads, default = 0")
parser.add_option("--S2",dest = "stage2", type = int, default = 1,
                  help = "number of resampling of each kmer table, default = 0")
parser.add_option("-s", dest = "sim", action = 'store_true',
                  help = "only print commands, do not run them")

(options, args) = parser.parse_args()

nThreads = options.nThreads
memory = options.memSize
###check the data directory:
if not os.path.isdir(options.dataDir):
    print 'Cannot find data directory {}'.format(options.dataDir)
    sys.exit(2)

if os.path.exists(options.dataDir + '/' + options.outFile):
    s = raw_input('{} is already in your data directory, overwrite it? Y/N '
                  .format(options.outFile))
    if s == 'Y' or s == 'y':
        print('{} is going to be overwritten'.format(options.outFile))
    elif s == 'N' or s == 'n':
        print('Quitting the program, please rerun with new output file name')
        sys.exit(2)
    else:
        print('Wrong keyboard input, exit')
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
    merge = './kmer_merge'
    if not is_exe(merge):
        print 'kmer_merge not found. Make sure it is in your PATH or the'
        print 'current directory, and that it is executable'
        sys.exit(1)
else:
    merge = 'kmer_merge'

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

###Set filter option
n = options.filter

###Get sample list:
samples = []
for fileName in os.listdir(options.dataDir):
    if os.path.isdir(os.path.join(options.dataDir, fileName)):
        samples.append(fileName)
samples.sort()

print 'SPECIES LIST:'
for sample in samples:
    print sample
if os.path.exists("./consensus_trees_total_nonparametric"):
    os.system("rm -f consensus_trees_read_nonparametric consensus_trees_table_nonparametric consensus_trees_total_nonparametric")
#os.system('touch consensus_trees_read_nonparametric')
#os.system('touch consensus_trees_total_nonparametric')

#############################
# STAGE 1: bootstrap over reads
#############################
if options.stage1 > 0:
    for s1 in xrange(options.stage1):
        print "{} out of {} times of bootstrap over reads.".format(s1+1,options.stage1)
        ###Prepare kmer_count jobs
        jobList = []
        if os.path.exists("./boot_Folder"):
            os.system("rm -r ./boot_Folder")
        os.system('mkdir boot_Folder')
        for sample in samples:
            os.system('mkdir boot_Folder/{}'.format(sample))
            outFile = os.path.join(options.dataDir, '{}.pkdat.gz'.format(sample))
            command = '{} -l {} -n {} -G {} -o {}'.format(kmerCount, options.kLen,
                                                n, memory/nThreads, outFile)
            for inputFile in os.listdir(os.path.join(options.dataDir, sample)):
                simFile = os.path.join(os.getcwd(),'boot_Folder',sample, inputFile)
                sim_file = smartopen(simFile, 'w')
                inputFile = os.path.join(options.dataDir, sample, inputFile)
                handle = smartopen(inputFile)
                firstChar = handle.read(1)
                seqdict = {}
                seq_id = 0
                if firstChar == '@':
                    seqFormat = 'fastq-illumina'
                elif firstChar == '>':
                    seqFormat = 'fasta'
                else:
                    seqFormat = False
                if seqFormat:
                    command += " -i '{}'".format(simFile)
                    for seq_record in SeqIO.parse(handle, seqFormat):
                        seqdict[seq_id] = seq_record.seq
                        seq_id += 1
                else:
                    print 'Error, file {} is not in {} format. Aborting!'.\
                       format(inputFile, options.seqFormat)
                    sys.exit(3)
                for i in xrange(len(seqdict)):
                    sim_id = random.randint(0,len(seqdict)-1)
                    sim_file.write('>{}_{}\n{}\n'.format(sample,sim_id,seqdict[sim_id]))
                sim_file.close()
            jobList.append(command)
        jobList = jobList[::-1]
        ###Run jobs
        pool = mp.Pool(nThreads)
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
                pool.close()
                pool.join()
                pool = mp.Pool(nThreads)
                nJobs = 0
            if jobList:
                command = jobList.pop()
                job = pool.apply_async(runJob, args=[command, options.sim])
                nJobs += 1
            else:
                break
            count += 1

        if nJobs:
            print '\n', time.strftime('%c')
            print "running last batch"
            pool.close()
            pool.join()
        ###Run aaf_merge
        outFile = os.path.join(options.dataDir, options.outFile)
        command = "{} -k s -c -d '0' -a 'O,M,F'".format(merge) #total kmer table
        cut = []
        for i, sample in enumerate(samples):
            command += " '{}.pkdat.gz'".format(os.path.join(options.dataDir, sample))
            cut.append(str((i + 1) * 2))
        if options.outFile.endswith('.gz'):
            command += ' | cut -f {} | gzip > {}'.format(','.join(cut), outFile)
        else:
            command += ' | cut -f {} > {}'.format(','.join(cut), outFile)
        print '\n', time.strftime('%c')
        print command
        if not options.sim:
            os.system(command)
        print time.strftime('%c')
        aaf_distance(outFile,nThreads,memory,samples,options.kLen)
        os.system('cat outtree >> consensus_trees_read_nonparametric')
        #############################
        # STAGE 2: bootstrap over table
        #############################
        kmerTable = smartopen(outFile)
        lines= kmerTable.readlines()
        table_len = len(lines)
        for s2 in xrange(options.stage2):
            print "{} out of {} times of bootstrap over table for the read bootstrap {}.".format(s2+1,options.stage1,s1+1)
            simTable = open('simTable_{}_{}'.format(s1,s2),'w')
            for x in xrange(table_len):
                line_num = random.randint(0,table_len-1)
                simTable.write(lines[line_num])
            simTable.close()
            aaf_distance('simTable_{}_{}'.format(s1,s2),nThreads*2,nThreads*2,samples,options.kLen)
            os.system('cat outtree >> consensus_trees_total_nonparametric')
        kmerTable.close()

#Skipping stage 1 (when s1==0)
elif options.stage1 == 0:
    ###Prepare kmer_count jobs
    jobList = []
    for sample in samples:
        outFile = os.path.join(options.dataDir, '{}.pkdat.gz'.format(sample))
        command = '{} -l {} -n {} -G {} -o {}'.format(kmerCount, options.kLen,
                                                    n, memory/nThreads, outFile)
        for inputFile in os.listdir(os.path.join(options.dataDir, sample)):
            inputFile = os.path.join(options.dataDir, sample, inputFile)
            handle = smartopen(inputFile)
            firstChar = handle.read(1)
            seqdict = {}
            seq_id = 0
            if firstChar == '@':
              seqFormat = 'FQ'
            elif firstChar == '>':
              seqFormat = 'FA'
            else:
              seqFormat = False
            if seqFormat:
                command += " -i '{}'".format(inputFile)
            else:
                print 'Error, file {} is not in {} format. Aborting!'.\
                      format(inputFile, options.seqFormat)
                sys.exit(3)
        jobList.append(command)
    jobList = jobList[::-1]
    ###Run jobs
    pool = mp.Pool(nThreads)
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
            pool.close()
            pool.join()
            pool = mp.Pool(nThreads)
            nJobs = 0
        if jobList:
            command = jobList.pop()
            job = pool.apply_async(runJob, args=[command, options.sim])
            nJobs += 1
        else:
            break
        count += 1

    if nJobs:
        print '\n', time.strftime('%c')
        print "running last batch"
        pool.close()
        pool.join()
    ###Run kmer_merge
    outFile = os.path.join(options.dataDir, options.outFile)
    command = "{} -k s -c -d '0' -a 'O,M,F'".format(merge) #total kmer table
    cut = []
    for i, sample in enumerate(samples):
        command += " '{}.pkdat.gz'".format(os.path.join(options.dataDir, sample))
        cut.append(str((i + 1) * 2))
    if options.outFile.endswith('.gz'):
        command += ' | cut -f {} | gzip > {}'.format(','.join(cut), outFile)
    else:
        command += ' | cut -f {} > {}'.format(','.join(cut), outFile)
    print '\n', time.strftime('%c')
    print command
    if not options.sim:
        os.system(command)
    print time.strftime('%c')
    #############################
    # STAGE 2: bootstrap over table
    #############################
    kmerTable = smartopen(outFile)
    lines= kmerTable.readlines()
    table_len = len(lines)
    for s2 in xrange(options.stage2):
        print "{} out of {} times of bootstrap over the real total kmer table.".format(s2+1,options.stage2)
        simTable = open('simTable_{}'.format(s2),'w')
        for x in xrange(table_len/options.kLen):
            line_num = random.randint(0,table_len-1)
            simTable.write(lines[line_num])
        simTable.close()
        aaf_distance('simTable_{}'.format(s2),nThreads,memory,samples,options.kLen)
        os.system('cat outtree >> consensus_trees_table_nonparametric')
    kmerTable.close()


#consense
if os.path.exists("./outfile"):
    os.system("rm -f outfile outtree")
if options.stage1 == 0:
    command = 'printf "consensus_trees_table_nonparametric\nY" | {} > /dev/null'.format(consense)
    print command
    os.system(command)
    os.system("mv outfile consense_outfile_table_nonparametric")
    os.system("mv outtree consensus_tree_table_nonparametric.tre")
elif options.stage2 == 0:
    command = 'printf "consensus_trees_read_nonparametric\nY" | {} > /dev/null'.format(consense)
    print command
    os.system(command)
    os.system("mv outfile consense_outfile_read_nonparametric")
    os.system("mv outtree consensus_tree_read_nonparametric.tre")
else:
    command = 'printf "consensus_trees_read_nonparametric\nY" | {} > /dev/null'.format(consense)
    print command
    os.system(command)
    os.system("mv outfile consense_outfile_read_nonparametric")
    os.system("mv outtree consensus_tree_read_nonparametric.tre")
    command = 'printf "consensus_trees_total_nonparametric\nY" | {} > /dev/null'.format(consense)
    print command
    os.system(command)
    os.system("mv outfile consense_outfile_total_nonparametric")
    os.system("mv outtree consensus_tree_total_nonparametric.tre")
os.system("rm simTable_*")
