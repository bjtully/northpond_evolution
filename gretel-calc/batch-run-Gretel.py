#!/usr/bin/env python3

'''
Process North Pond Evolution MAGs to determine possible haplotypes
for each gene in the genome using Gretel
'''

import glob
import subprocess
import os
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor, as_completed


def contig_prep(record, mag, tps):
	#Process the contig FASTA for ID numbers to run gretel-snpper and create necessary
	#VCF files
	#Iterates through all timepoints from the timepoints of interest list
	
	#As part of the MASTER input during the Gretel run, Gretel expects a FASTA
	#file containing only the target contig, otherwise, filling in the sequence
	#fails
	if os.path.exists('./'+mag+'/'+record.id+".fasta") == False:
		outcontigfasta = open('./'+mag+'/'+record.id+".fasta", "w")
		outcontigfasta.write(record.format("fasta"))
	for i in tps:
		outvcfname = '82A_'+i+'-'+record.id+'.vcf'
		if os.path.exists('./'+mag+'/'+outvcfname+'.gz') == False:
			outvcffile = open('./'+mag+'/'+outvcfname, "w")
			subprocess.call(['gretel-snpper', '--bam', '82A_'+i+'.sorted_filtered.bam', '--contig', record.id], stdout=outvcffile)
			subprocess.call(['bgzip', './'+mag+'/'+outvcfname])
			subprocess.call(['tabix', './'+mag+'/'+outvcfname+'.gz'])
	return

def run_gretel(mag, tps, contig, geneid, s, e, rev):
	for i in tps:
		greteldataout = "CONTIG-"+contig+"-GENECALL-"+str(geneid)+"-TIMEPOINT-"+i
		#Gretel needs an existing output directory
		if os.path.exists('./'+mag+'/'+greteldataout) == False:
			subprocess.call(['mkdir', './'+mag+'/'+greteldataout])
			subprocess.call(['gretel', '82A_'+i+'.sorted_filtered.bam', 
				'./'+mag+'/82A_'+i+'-'+contig+'.vcf.gz', contig, 
				'-s', s, '-e', e, '--master', './'+mag+'/'+contig+'.fasta', 
				'-o', './'+mag+'/'+greteldataout])
		#If Gretel fails, remove the output directory -- maybe we should be saving which gene calls
		#and timepoints do not have alternate paths?
		#if os.path.exists('./'+mag+'/'+greteldataout+'/gretel.crumbs') == False:
		#	subprocess.call(['rm', '-r', './'+mag+'/'+greteldataout+"/"])
		#Interested in protein sequences with variances
		#else:
			#Identified by Prodigal, if true need to reverse complement first
			if os.path.exists('./'+mag+'/'+greteldataout+'/gretel.crumbs') == True:
				if rev == True:
					subprocess.call(['seqmagick', 'convert', '--reverse-complement', 
						'./'+mag+'/'+greteldataout+'/out.fasta', 
						'./'+mag+'/'+greteldataout+'/'+greteldataout+'.reverse.fasta'])
					subprocess.call(['seqmagick', 'convert', '--translate', 'dna2protein', 
						'./'+mag+'/'+greteldataout+'/'+greteldataout+'.reverse.fasta', 
						'./'+mag+'/'+greteldataout+'/'+greteldataout+'.proteins.faa'])
				if rev == False:
					subprocess.call(['cp', './'+mag+'/'+greteldataout+'/out.fasta', 
						'./'+mag+'/'+greteldataout+'/'+greteldataout+'.forward.fasta'])
					subprocess.call(['seqmagick', 'convert', '--translate', 'dna2protein', 
						'./'+mag+'/'+greteldataout+'/'+greteldataout+'.forward.fasta', 
						'./'+mag+'/'+greteldataout+'/'+greteldataout+'.proteins.faa'])
	return

if __name__ == "__main__":
	#this contains a simple text file with the genome ID and the NP timepoints of interest
	for line in open("timepoints-of-interest.txt", "r"):
		if line[0] != "#":
			data = line.strip().split("\t")
			magid = data[0]
			timepoints = data[1].split(",")
			#Track each set of results in their own directory
			subprocess.call(['mkdir', "./"+magid])
			#It seems like anvio gene calls are identical to Prodigal, but Prodigal resets at 1 on each contig
			#and anvio starts at zero and counts continuously, will use anvio numbering so
			#that results might be able to be paired in the future
			anvigeneid = 0
			convert = {}
			for line in open(magid+".prodigal.genecalls.fasta", "r"):
				if line[0] == ">":
					data = line.strip()[1:].split(" # ")
					convert[data[0]] = anvigeneid
					anvigeneid += 1

			with ThreadPoolExecutor(max_workers=15) as executor:
				futures1 = []
				for record in SeqIO.parse(open(magid+".fasta", "r"), "fasta"):
					futures1.append(executor.submit(contig_prep, record, magid, timepoints))
				futures2 = []
				for future1 in as_completed(futures1):
					result1 = future1.result()
					error1 = future1.exception()
					if error1 is None:
						for line in open(magid+".prodigal.genecalls.fasta", "r"):
							if line[0] == ">":
								reverse = False
								data = line.strip()[1:].split(" # ")
								#Skip genes determined to be partial by Prodigal
								if "1" in data[4].split(";")[1]:
									continue
								targetcontig = "_".join(data[0].split("_")[:2])
								geneid = convert[data[0]]
								start = data[1]
								end = data[2]
								if data[3] == "-1":
									reverse = True
								futures2.append(executor.submit(run_gretel, magid, timepoints, targetcontig, geneid, start, end, reverse))
					else:
						print(error1)
				for future2 in as_completed(futures2):
					result2 = future2.result()
					error2 = future2.exception()
					if error2 is None:
						print(result2)
					else:
						print(error2)