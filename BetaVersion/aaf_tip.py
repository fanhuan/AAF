#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  aaf_distance.py
#
#  Copyright 2013, 2014 Huan Fan <hfan22@wisc.edu> & Yann Surget-Groba
#  <yann@xtbg.org.cn>

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

from optparse import OptionParser
import os, sys, subprocess, math, re

Usage = "%prog [options] -i <input tree file> -k <kmer size> " + \
        "--tip <tip setting file>"
version = '%prog 20171109.1'
parser = OptionParser(usage = Usage, version = version)
parser.add_option("-i", dest = "iptf", help = "tree file to be trimmed", \
                  default = "aaf.tre")
parser.add_option("-k", dest = "klen", type = int, \
                  help = "kmer size used for constructing the input tree")
parser.add_option("--tip", dest = "tip_file", default = "tip_file_test.txt", \
                  help = "tip setting file, default = tip_file_test.txt")
parser.add_option("-n", dest = "filter", action = 'store_true', \
                  help = "k-mer filtering was on for tree construction")
parser.add_option("-f", dest = "countf", default = "kmer_diversity.wc", \
                  help = "k-mer diversity file, default = kmer_diversity.wc")
(options, args) = parser.parse_args()

kl = options.klen
n = options.filter
sl = [] #species list
info = {}

#process count file
with open(options.countf) as total:
    lines = total.readlines()

sn = len(lines)
ntotal = [0.0] * sn

i = 0
for line in lines:
    ntotal[i] = float(line.split()[1])



# process tip file
with open(options.tip_file, 'rU') as tip_info:
    tips = tip_info.readlines()

if sn == len(tips) - 1:
    tips = tips[1:]
    for line in tips:
        line = line.split()
        if len(line) == 4:
            sl.append(line[0])
            info[line[0]] = line[1:]
        else:
            print("tip-info not in the right format, exit")
            sys.exit()
else:
    print("kmer diversity file does not match with the tip info file, exit")
    sys.exit()

if not info:
    print("tip-info file empty, exit")
    sys.exit()

tip = {}

'''
Previous wrong equations:
    if n: # with filter
        tip[key] = -1 / float(kl) * math.log(1 - math.exp(-L * (1 - e) ** \
                   float(kl)) - L * math.exp(-L * (1 - e) ** float(kl)) * \
                   (1 - e) ** float(kl))
    else: # no filter
        tip[key] = 1 / float(kl) * math.log((1 - math.exp(-L * (1 - e) ** \
                   float(kl)) + L * (1 -(1 - e) ** float(kl))) / ((1 - \
                   math.exp(-L * (1 - e) ** float(kl))) ** 2))
'''
for key in info:
    c = float(info[key][0])
    r = float(info[key][1])
    e = float(info[key][2])
    L = c / r * (r - float(kl) + 1)
    Pta = L * (1 - (1 - e) * float(kl))
    Pr = 1 - math.exp(-L)
    Prf = 1 - (1 + L) * math.exp(-L)
    #Eq.6
    Pe = 1 - (math.exp(-L * (1-e)**float(kl)) - math.exp(-L))/(1 - math.exp(-L))
    #Eq.7
    Pef = 1 - (math.exp(-L * (1-e)**float(kl)) - (1+L) * math.exp(-L))/(1-(1+L) \
        * math.exp(-L)) - (1-e)**float(kl) * L * math.exp(-L*(1-e)**float(kl)) \
        /(1-(1+L)* math.exp(-L))
    if n: # with filter
        tip[key] = 0.5*(1/float(kl)) * math.log(1/(Prf * Pef))
    else: # no filter
        #Eq.10
        tip[key] = 0.5/float(kl) * math.log((Pr*Pe + Pta)/(Pr**2 * Pe**2))
# calculate the tip to trim
tip_total = 0
for i in range(sn):
    for j in range(i + 1, sn):
        if ntotal[i] > ntotal[j]:
            tip_total += tip[sl[j]]
        else:
            tip_total += tip[sl[i]]
tip_ave = tip_total / (sn * (sn - 1))

#trim the tips in input tree file
tree_handle = open(options.iptf)
tree = tree_handle.read()
tree_handle.close()
tree_list = tree.split(',')
tree_tips = []
for species in tree_list:
    match = re.search('[A-za-z]*\:\d+\.\d*', species)
    bl = match.group(0).split(':')[1] #bl for branchlength
    bl_tip = str(float(bl) - tip_ave)
    species = species.replace(bl, bl_tip)
    tree_tips.append(species)

Dtip = open('tip_' + options.iptf.split('/')[-1], 'w')
Dtip.write(','.join(tree_tips))
Dtip.close()
