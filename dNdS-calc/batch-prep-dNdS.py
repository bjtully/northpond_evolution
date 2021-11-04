#!/usr/bin/env python3

'''
Process the Gretel outputs to calculate dN/dS for haplotypes
across time points
 - PAL2NAL

Also, for each gene, want to calculate the number of 100% unique
haplotypes and isoforms
 - cd-hit

'''

import glob
import subprocess
import os
from Bio import SeqIO
from Bio.Seq import Seq
from concurrent.futures import ThreadPoolExecutor, as_completed

def run_cdhit(geneid, seqrecord, magid):
	storedgenes = []
	storedproteins = []
	###### COUNTS
	totalgretelhaplo = 0 #total predicted gretel haplotypes
	translateableiso = 0 #only haplotypes that result in protein without stop codons
	haploclusters = 0 #number of 100% cd-hit haplotype clusters
	isoclusters = 0 #number of 100% cd-hit isoforms clusters
	######
	#Stores the reference gene
	seqrecord.id = "GENECALL-"+geneid+"-REFERENCE"
	seqrecord.description = ""
	storedgenes.append(seqrecord)
	#Stores the reference protein after translation
	protrecord = seqrecord.translate(table=11)
	protrecord.id = "GENECALL-"+geneid+"-REFERENCE"
	protrecord.description = ""
	storedproteins.append(protrecord)
	tphaplo_in = glob.glob("./"+magid+"/"+"*GENECALL-"+geneid+"-*/"+"CONTIG*.fasta")
	tpisoform_in = glob.glob("./"+magid+"/"+"*GENECALL-"+geneid+"-*/"+"CONTIG*.faa")
	genedirectory = ""
	if len(tpisoform_in) > 0:
		#Skip proteins that have stop codons. Also do not store corresponding
		#genecalls
		#Rename to <20 characters for codeml downstream
		genedirectory = tphaplo_in[0].split("/")[2][:-13]
		subprocess.call(['mkdir', "./dNdS-"+magid+"/"+genedirectory])
		skipgenes = []
		for i in tpisoform_in:
			for record in SeqIO.parse(open(i, "r"), "fasta"):
				if "*" not in record.seq[:-1]:
					record.id = "GENECALL-"+i.split("-")[3]+"-"+i.split("/")[2][-2:]+"_"+record.id.split("__")[0]
					record.description = ""
					storedproteins.append(record)
				else:
					skipgenes.append(record.id)
		for h in tphaplo_in:
			for record in SeqIO.parse(open(h, "r"), "fasta"):
				totalgretelhaplo += 1
				if record.id not in skipgenes:
					record.id = "GENECALL-"+h.split("-")[3]+"-"+h.split("/")[2][-2:]+"_"+record.id.split("__")[0]
					record.description = ""
					storedgenes.append(record)
		combinedproteins = "./dNdS-"+magid+"/"+genedirectory+"/"+genedirectory+".proteins.faa"
		combinedgenes = "./dNdS-"+magid+"/"+genedirectory+"/"+genedirectory+".genecalls.fasta"
		if os.path.exists(combinedgenes) == False:
			SeqIO.write(storedproteins, combinedproteins, "fasta")
			SeqIO.write(storedgenes, combinedgenes, "fasta")
		
		#### CD-HIT

		subprocess.call(['cd-hit', '-i', combinedproteins, '-o', "./dNdS-"+magid+"/"+genedirectory+"/ISOFORMS", '-c', '1', '-d', '0', '-g', '1'])
		subprocess.call(['cd-hit', '-i', combinedgenes, '-o', "./dNdS-"+magid+"/"+genedirectory+"/HAPLOTYPES", '-c', '1', '-d', '0', '-g', '1'])
		####CD-HIT HAPLOTYPE CLUSTER OUTPUT
		for line in open("./dNdS-"+magid+"/"+genedirectory+"/HAPLOTYPES.clstr", "r"):
			if line[0] == ">":
				haploclusters += 1
		####CD-HIT ISOFORM CLUSTER OUTPUT
		for line in open("./dNdS-"+magid+"/"+genedirectory+"/ISOFORMS.clstr", "r"):
			if line[0] == ">":
				isoclusters += 1
			if line[0] != ">":
				translateableiso += 1

		###CODEML-prep
		### Each run of CODEML requires its own control file
		### After running pal2nal alignment stage, create a control file for each gene
		proteinaln = "./dNdS-"+magid+"/"+genedirectory+"/"+genedirectory+".proteins.aln"
		subprocess.call(['muscle', '-in', combinedproteins, '-out', proteinaln])
		codonalign = "./dNdS-"+magid+"/"+genedirectory+"/"+genedirectory+".codonalign.fa"
		codonalignout = open(codonalign, "w")
		subprocess.call(['pal2nal.pl', proteinaln, combinedgenes, '-codontable', '11', '-output', 'fasta'], stdout=codonalignout)
		genecontrolfile = "./dNdS-"+magid+"/"+genedirectory+"/codeml.ctl"
		genecontrolwrite = open(genecontrolfile, "w")
		for line in open("codeml-2.ctl", "r"):
			if "seqfile" in line:
				genecontrolwrite.write("      seqfile = "+codonalign+"\n")
			if "outfile" in line:
				genecontrolwrite.write("      outfile = "+"./dNdS-"+magid+"/"+genedirectory+"/mlc"+"\n")
			if "seqfile" not in line and "outfile" not in line:
				genecontrolwrite.write(line)

	return {geneid:[totalgretelhaplo, translateableiso, haploclusters, isoclusters]}



if __name__ == "__main__":
	for line in open("timepoints-of-interest.txt", "r"):
		if line[0] != "#":
			data = line.strip().split("\t")
			magid = data[0]
			#timepoints = data[1].split(",")
			#Track each set of results in their own directory
			subprocess.call(['mkdir', "./dNdS-"+magid])
			greteloutput = {}
			
			with ThreadPoolExecutor(max_workers=20) as executor:
				futures1 = []
				for record in SeqIO.parse(open(magid+".anvio.genecalls.fasta", "r"), "fasta"):
					futures1.append(executor.submit(run_cdhit, record.id, record, magid))
				for future1 in as_completed(futures1):
					#print(future1.result())
					greteloutput.update(future1.result())
					error1 = future1.exception()
			haplotypeinfo = open("./"+magid+".gretelsummary.tsv", "w")
			haplotypeinfo.write("GeneID\tGretel Haplotypes\tTranslatable haplotypes\tHaplotype Clusters\tIsoform Clusters\n")
			for k in greteloutput:
				haplotypeinfo.write(k+"\t"+"\t".join(str(x) for x in greteloutput[k])+"\n")

					