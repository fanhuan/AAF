#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  aaf_kmercount.py
#
#  Copyright 2017 Huan Fan <hfan22@wisc.edu>
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

import sys, os, time, math, psutil
from concurrent.futures import ProcessPoolExecutor as PPE
import numpy as np


version = '%prog 20171211.1'


'''
run_command
function included:
smartopen
is_exe
present
countTotal
countShared_single
countTotal_shared
aaf_kmercount
aaf_dist

'''
def run_command(command):
    print(command)
    print(time.strftime('%c'))
    try:
        os.system(command)
    except:
        print("Unexpected error:", sys.exc_info()[0])
        raise

def smartopen(filename, mode = 'rt'):
    import gzip, bz2
    if filename.endswith('gz'):
        return gzip.open(filename, mode)
    elif filename.endswith('bz2'):
        return bz2.BZ2File(filename, mode)
    else:
        return open(filename,mode)

def is_exe(fpath):
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

def present(x,n=1):
    if int(x) >= n:
        return '1'
    else:
        return '0'

def countTotal(lines):
    line_list = []
    for line in lines:
        line_list.append([int(present(i)) for i in line.split()[1:]])
    line_total = np.sum(line_list,axis = 0)
    return line_total

def countShared_single(line):
    '''takes one line from the shared kmer table(skt) and count nshare.
    Previously this function is called countShared and it takes a chunk of
    skt.
    '''
    line = line.split()
    if line[0][0].isdigit():
        sn = len(line)
        flag = 'd'
    else:
        sn = len(line) - 1
        flag = 'k'
    shared = [[0] * sn for i in range(sn)]
    if flag == 'k':
        line = line[1:]
    line = [int(i) for i in line]
    for i in range(sn):
        for j in range(i + 1, sn):
            if line[i]*line[j] != 0:
                shared[i][j] += 1
    return shared

def countTotal_shared(lines,sn):
    line_list = []
    shared = [[0] * sn for i in range(sn)]
    for line in lines:
        line = line.split()
        if len(line) == sn+1:
            line = line[1:]
        line = [int(i) for i in line]
        line_list.append([int(present(i)) for i in line])
        for i in range(sn):
            for j in range(i + 1, sn):
                if line[i]*line[j] != 0:
                    shared[i][j] += 1
    line_total = np.sum(line_list,axis = 0)
    return (line_total,shared)

def aaf_kmercount(dataDir,k,n,nThreads,memPerThread):
    #check excutables
    if k > 25:
        if os.system('which kmer_countx > /dev/null'):
            kmerCount = './kmer_countx'
            if not is_exe(kmerCount):
                print('kmer_countx not found. Make sure it is in your PATH or the')
                print('current directory, and that it is executable')
                sys.exit(1)
        else:
            kmerCount = 'kmer_countx'

    else:
        if os.system('which kmer_count > /dev/null'):
            kmerCount = './kmer_count'
            if not is_exe(kmerCount):
                print('kmer_count not found. Make sure it is in your PATH or the')
                print('current directory, and that it is executable')
                sys.exit(1)
        else:
            kmerCount = 'kmer_count'

    ###Get sample list:
    samples = []
    for fileName in os.listdir(dataDir):
        if os.path.isdir(os.path.join(dataDir, fileName)):
            samples.append(fileName)
        else:
            if not fileName.startswith('.'):
                sample = fileName.split(".")[0]
                if sample in samples:
                    sample = fileName.split(".")[0]+fileName.split(".")[1]
                    if sample in samples:
                        print('Error, redundant sample or file names. Aborting!')
                        sys.exit(3)
                os.system("mkdir {}/{}".format(dataDir,sample))
                os.system("mv {}/{} {}/{}/".format(dataDir,fileName,dataDir,sample))
                samples.append(sample)
    samples.sort()
    print(time.strftime('%c'))
    print('SPECIES LIST:')
    for sample in samples:
        print(sample)

    ###Prepare kmer_count jobs
    jobList = []
    for sample in samples:
        outFile = '{}.pkdat.gz'.format(sample)
        command = '{} -l {} -n {} -G {} -o {} -f '.format(kmerCount, k, n,
                    memPerThread, outFile)
        command1 = ''
        for inputFile in os.listdir(os.path.join(dataDir, sample)):
            inputFile = os.path.join(dataDir, sample, inputFile)
            handle = smartopen(inputFile)
            firstChar = handle.read(1)
            if firstChar == '@':
                seqFormat = 'FQ'
            elif firstChar == '>':
                seqFormat = 'FA'
            else:
                print('Error, file {} is not FA or FQ format. Aborting!'.format(inputFile))
                sys.exit(3)
            command1 += " -i '{}'".format(inputFile)
        command += '{}{}> {}.wc'.format(seqFormat,command1,sample)
        jobList.append(command)
    ###Run jobs
    with PPE(max_workers = nThreads) as executor:
        executor.map(run_command,jobList)
    return samples


def aaf_dist(datfile,countfile,nThreads,samples,kl,long=False):
    #check executables
    if os.system('which fitch_kmerX > /dev/null'):
        if long:
            fitch = './fitch_kmerX_long'
        else:
            fitch = './fitch_kmerX'
        if not is_exe(fitch):
            print(fitch+' not found. Make sure it is in your PATH or the')
            print('current directory, and that it is executable')
            sys.exit()
    else:
        if long:
            fitch = 'fitch_kmerX_long'
        else:
            fitch = 'fitch_kmerX'
    #process the .dat.gz file
    try:
        iptf = smartopen(datfile,'rt')
    except IOError:
        print('Cannot open file', datafile)
        sys.exit()

    if not os.path.isfile(countfile):
        print('Cannot find file', countfile)
        sys.exit()

    try:
        total = open(countfile,'rt')
    except IOError:
        print('Cannot open file', countfile)
        sys.exit()

    try:
        infile = open('infile','wt')
    except IOError:
        print('Cannot open infile for writing')
        sys.exit()

    ###Read header
    sl = []                 #species list
    line = iptf.readline()
    ll = line.split()
    if kl != float(ll[1]):
        print("The recorded k in the shared kmer table file is not the same with the k supplied to aaf_dist; exiting now.")       #kmer length
        sys.exit()
    while True:
        line = iptf.readline()
        if line.startswith('#-'):
            continue
        elif line.startswith('#sample'):
            ll = line.split()
            sl.append(ll[1])
        else:
            break
    if sl != samples:
        print("The recorded sample list in the shared kmer table file is not the same with the one supplied to aaf_dist; exiting now.")       #kmer length
        sys.exit()
    ###Initialize shared kmers matrix
    sn = len(samples)    #species number
    nshare = [[0] * sn for i in range(sn)]

    ### It turns out to be very slow if we give very big chunks. So we will be
    ### Using only 1G of RAM intotal. As a result, we could use all the cores
    ### available, which was not possible for the kmer_count step.
    cpu_num = psutil.cpu_count()

    ###Compute the number of lines to process per thread (chunk size)
    line = iptf.readline()
    line_size = sys.getsizeof(line)
    chunkLength = int(1024 ** 3 / (cpu_num/2) / line_size)
    print('chunkLength = {}'.format(chunkLength))
    while True:
        lines = []
        for nLines in range(chunkLength):
            if not line: #if empty
                break
            lines.append(line)
            line = iptf.readline()
        if not lines: #if empty
            break
        ###Compute shared kmer matrix
        with PPE(max_workers = int(cpu_num/2)) as executor:
            for result in executor.map(countShared_single,lines):
                for i in range(sn):
                    for j in range(i + 1, sn):
                        nshare[i][j] += result[i][j]

    iptf.close()

    ###Compute distance matrix
    ntotal = [0.0] * sn

    for i in range(sn):
        ntotal[i] = float(total.readline().split()[1])
    dist = [[0] * sn for i in range(sn)]

    for i in range(sn):
        for j in range(i + 1, sn):
            mintotal = min(ntotal[i], ntotal[j])
            if nshare[i][j] == 0:
                dist[j][i] = dist[i][j] = 1
            else:
                distance = (-1 / float(kl) * math.log(nshare[i][j] / mintotal))
                #print(mintotal,nshare[i][j])
                dist[j][i] = dist[i][j] = distance
                nshare[j][i] = nshare[i][j]

    total.close()

    ###Write infile
    infile.write('{} {}'.format(sn, sn))
    namedic = {}
    for i in range(sn):
        lsl = len(sl[i])
        if lsl < 10:
            ssl = sl[i] + ' ' * (10 - lsl)
        else:
            ssl = sl[i][:10]
            appendix = 1
            while ssl in namedic:
                ssl = ssl[:-len(str(appendix))] + str(appendix)
                appendix += 1
        namedic[ssl] = sl[i]
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
    fh1 = open(datfile.split('.')[0]+'.tre','wt')

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
    command = 'mv infile {}.dist'.format(datfile.split('.')[0])
    os.system(command)

    os.system('rm -f outfile outtree')

    print('{} end'.format(time.strftime("%c")))
