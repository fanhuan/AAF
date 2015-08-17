#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  pkdat2hist.py
#
#  Copyright 2013, 2014 Huan Fan <hfan22@wisc.edu>

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

from optparse import OptionParser
import os, gzip, collections, operator

def smartopen(filename,*args,**kwargs):
    '''opens with open unless file ends in .gz, then use gzip.open

    in theory should transparently allow reading of files regardless of compression
    '''
    if filename.endswith('.gz'):
        return gzip.open(filename,*args,**kwargs)
    else:
        return open(filename,*args,**kwargs)
Usage = "pkdat2hist.py [ -i <input filename>]"
parser = OptionParser(Usage)
parser.add_option( "-i", dest="iptf", help="input the pkdat file")
(options, args) = parser.parse_args()

iptf = smartopen(options.iptf)
#n=options.n
output = open(options.iptf.split('.')[0]+'.hist','w')

counts_int = []
for line in iptf:
    counts_int.append(int(line.split()[1]))
#abundance=collections.Counter(counts).most_common(n) #this returns a list
abundance_dic=collections.Counter(counts_int) #this returns a dic
n=len(abundance_dic)
abundance=abundance_dic.most_common(n)
sorted_ab=sorted(abundance, key=operator.itemgetter(0)) #sorted according to frequency, ascending.
for item in sorted_ab:
    if item[0] < 10001:
        output.write("%d\t%d\n" % (int(item[0]),item[1]))

output.close()
