#!/usr/bin/env python3

'''
Using the code from Crits-Christoph as a guide,
computing the Hudson Fst score for genes that have
at least 3 alleles in pairs of time points.
'''

import numpy as np 
import pandas as pd 
from tqdm import tqdm
import allel
import argparse
import os
from collections import defaultdict
import copy
from Bio import SeqIO
import glob

def create_gene_index(gene_fasta): #Code borrowed from https://github.com/alexcritschristoph/soil_popgen
    gene_index = []
    complete_genes = 0.0
    partial_genes = 0.0
    for record in SeqIO.parse(gene_fasta, 'fasta'):
        gene = str(record.id)
        if 'partial=00' in record.description:
            complete_genes += 1
            gene_scaf = "_".join(gene.split("_")[:-1])
            # NOTE: PRODIGAL USES A 1-BASED INDEX AND WE USE 0, SO CONVERT TO 0 HERE 
            gene_start = int(record.description.split("#")[1].strip())-1
            gene_end = int(record.description.split("#")[2].strip())-1
            gene_index.append({"name":gene, 'scaf': gene_scaf, "start": gene_start, "end": gene_end})

        else:
            partial_genes += 1

    print("Notice: " + str(round(partial_genes *100 / (complete_genes + partial_genes))) + "% of genes were incomplete and snps in these genes were marked I. Pi is not calculated for incomplete genes.")
    return gene_index


def process_SNV_tsv(infile): #process inStrain SNV files
    allele_counts = {}
    for line in tqdm(infile, desc= "Reading SNV tsv file"):
        if line[:8] != "scaffold":
            data = line.strip().split()
            if int(data[3]) == 2:
                allele_np = np.array(data[10:14])
                try: #try statement prevents an error unique to SNVs in NORP147 where allele_count = 2, but ref_freq is missing causing a shifted data list
                    allele_np = allele_np.astype('i1')
                except ValueError:
                    continue
                allele_counts[data[0] + ":" + data[1]] = (allele_np)
    return allele_counts


def main():
    outdata = {}
    prodigalfiles = glob.glob("*prodigal.genecalls.fna") #This assumes that the first element of the prodigal matches the genomeID
    for genecalls in prodigalfiles:
        genomeid = genecalls.split(".")[0]
        timepoints = glob.glob("./inStrain-Results/"+genomeid+"*IS/output/*SNVs.tsv")
        
        allele_count_dict = {}
        for IS_SNV_file in timepoints: #process each inStrain SNV file
            tp = IS_SNV_file.split("/")[-1].split(".")[0].split("-")[-1]
            snvtp = open(IS_SNV_file)
            allele_counts = process_SNV_tsv(snvtp)
            snvtp.close()
            allele_count_dict[tp] = allele_counts

        #Pairwise analysis
        for i,x in enumerate(allele_count_dict.keys()):
            for j,y in enumerate(allele_count_dict.keys()):
                if i > j:
                    
                    allele_counts1 = allele_count_dict[x]
                    allele_counts2 = allele_count_dict[y]

                    combined_snp_table = {'scaffold':[], 'position':[]}
                    for k in allele_counts1:
                        data = k.split(":")
                        combined_snp_table['scaffold'].append(data[0])
                        combined_snp_table['position'].append(int(data[1]))
                        if k not in allele_counts2.keys():
                            allele_np = np.array([0,0,0,0])
                            allele_np = allele_np.astype('i1')
                            allele_counts2[k] = allele_np
                    for k in allele_counts2:
                        if k not in allele_counts1.keys():
                            data = k.split(":")
                            combined_snp_table['scaffold'].append(data[0])
                            combined_snp_table['position'].append(int(data[1]))
                            allele_np = np.array([0,0,0,0])
                            allele_np = allele_np.astype('i1')
                            allele_counts1[k] = allele_np

                    SNPTable = pd.DataFrame(combined_snp_table)
                    FstTable = defaultdict(list)

                    filegenecalls = open(genecalls, "r")
                    for gene in tqdm(create_gene_index(filegenecalls), desc="calculating fst"):
                        snps = SNPTable[(SNPTable['scaffold'] == gene['scaf']) & (SNPTable['position'].astype(int) >= gene['start']) & (SNPTable['position'].astype(int) <= gene['end'])]
                        snp_list = []
                        for index, row in snps.iterrows():
                            snp_list.append(row['scaffold'] + ":" + str(row['position']))

                        if len(snp_list) >= 3:
                            allele_counts_1 = []
                            allele_counts_2 = []
                            for snp in snp_list:
                                allele_counts_1.append(allele_counts1[snp])
                                allele_counts_2.append(allele_counts2[snp])
                                
                            allel1 = allel.AlleleCountsArray(allele_counts_1)
                            allel2 = allel.AlleleCountsArray(allele_counts_2)
                            
                            fst_h = allel.moving_hudson_fst(allel1, allel2, size=len(snp_list))[0] #allel.moving_hudson_fst(a1,a2, size=3)
                            nd_1 = np.sum(allel.mean_pairwise_difference(allel1)) / (1 + gene['end'] - gene['start'])
                            nd_2 = np.sum(allel.mean_pairwise_difference(allel2)) / (1 + gene['end'] - gene['start'])
                            
                            FstTable['gene'].append(gene['name'])
                            FstTable['snp_num'].append(len(snp_list))
                            #print(fst_h)
                            if float(fst_h) < 0.0 or str(fst_h) == "nan":
                                FstTable['fst'].append(0)
                            if 1 >= float(fst_h) >= 0.0:
                                FstTable['fst'].append(fst_h)
                            if float(fst_h) > 1:
                                FstTable['fst'].append(1)
                            FstTable['pi_1'].append(nd_1)
                            FstTable['pi_2'].append(nd_2)
                            FstTable['cov_1'].append(np.mean(np.sum(allele_counts_1, axis=1)))
                            FstTable['cov_2'].append(np.mean(np.sum(allele_counts_2, axis=1)))

                    FstTable = pd.DataFrame(FstTable)
                    outid = genomeid+"_"+x+"-"+y
                    outdata[outid] = str(np.mean(FstTable['fst']))
                    FstTable.to_csv(outid+'.Fst.tsv', index=False, sep='\t')

    for i in outdata:
        print(i, outdata[i])

if __name__ == "__main__":
	main()