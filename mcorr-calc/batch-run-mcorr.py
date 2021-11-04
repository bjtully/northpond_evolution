#!/usr/bin/env python3

'''
Batch run mcorr
'''

import subprocess

timepoints = {}
for line in open("timepoints-of-interest.txt", "r"):
	if line[0] != "#":
		data = line.strip().split("\t")
		timepoints[data[0]] = data[1].split(",")

for genome in timepoints:
	for tp in timepoints[genome]:
		subprocess.call(['mcorr-bam', genome+'.gff', '82A_'+tp+'.sorted_filtered.bam', genome+'_'+tp])
        subprocess.call(['mcorr-fit', genome+'_'+tp+'.csv', genome+'_'+tp])
        