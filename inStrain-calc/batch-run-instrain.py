#!/usr/bin/env python3

'''
Batch run inStrain for MAGs in relevant time points
-Using parameters from the rest of the manuscript
 - ANI = 0.95%
 - coverage = 20
 - min_freq = 0.1

'''

import subprocess

timepoints = {}
for line in open("timepoints-of-interest.txt", "r"):
    if line[0] != "#":
        data = line.strip().split("\t")
        timepoints[data[0]] = data[1].split(",")

for genome in timepoints:
    for tp in timepoints[genome]:
        subprocess.call(['inStrain', 'profile', '82A_'+tp+'.sorted_filtered.bam', genome+'.fasta', '-o', genome+'-'+tp+'.IS', '-p', '20', '-g', genome+'.prodigal.genecalls.fna', '-l', '0.95', '--min_mapq', '1', '-c', '20', '-f', '0.1'])
