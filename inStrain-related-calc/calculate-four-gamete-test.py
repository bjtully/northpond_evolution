#!/usr/bin/env python3

'''
Process inStrain linkage and SNV data to to determine the results of
the four-gamete test for each genome as in Bouma-Gregson et al (2021)

From the SNV file, filter out all SNVs where number of alleles != 2
Process the linkage file, calculate the frequency of linked bi-allelic
SNVs that have 2 haplotypes (AB, ab), 3 haplotypes, or 4 haplotypes
(AB,ab, Ab, aB).

The frequency of the type 4H haplotypes is a metric of recombination
'''

import glob

# genomes = [
#     'NORP6',
#     'NORP57',
#     'NORP100',
#     'NORP167',
#     'NORP83',
#     'NORP139',
#     'NORP147',
#     'NORP163',
#     'NORP169',
#     'NORP246'
# ]

genomes = ['NORP83']

print("Genome-TP\tH1\tH2\tH3\tH4\tTotalBiallelicSites")
for genome in genomes:
    instrainDirectories = glob.glob("./inStrain-Results/"+genome+"-T*.IS") #identify inStrain output files
    for timepoint in instrainDirectories:
        #track SNV sites that have an allele_count != 2
        notbiallelicsites = {}
        gametefreq = {genome:{}}
        sampleID = timepoint.split("/")[2].split(".")[0]
        for line in open(timepoint+"/output/"+sampleID+".IS_SNVs.tsv", "r"): #read in inStrain SNV table
            data = line.strip().split("\t")
            if data[0] != "scaffold":
                if int(data[3]) != 2: #Determin which SNV sites are not biallelic
                    try:
                        notbiallelicsites[data[0]].append(str(data[1]))
                    except KeyError:
                        notbiallelicsites[data[0]] = [str(data[1])]
        for line in open(timepoint+"/output/"+sampleID+".IS_linkage.tsv", "r"): #read in inStrain linkage table
            data = line.strip().split("\t")
            if data[0] != "scaffold":
                if data[0] in notbiallelicsites.keys(): #exclude sites that are not biallelic
                    if str(data[1]) in notbiallelicsites[data[0]] or str(data[2]) in notbiallelicsites[data[0]]:
                        continue
                    else: #if biallelic calculate the haplotypes
                        allelecounts = [int(x) for x in data[12:16]]
                        numalleles = 0
                        for i in allelecounts:
                            if i > 0:
                                numalleles += 1
                        try:
                            gametefreq[genome][numalleles] += 1
                        except KeyError:
                            gametefreq[genome].update({numalleles:1})
                else: #for biallelic sites determine how many haplotypes are in data
                    allelecounts = [int(x) for x in data[12:16]]
                    numalleles = 0
                    for i in allelecounts:
                        if i > 0:
                            numalleles += 1
                    try:
                        gametefreq[genome][numalleles] += 1
                    except KeyError:
                        gametefreq[genome].update({numalleles:1})
        
        haplotypes = [1,2,3,4]
        total = 0
        for i in gametefreq[genome]: #calculate total and determine fraction of different haplotypes
            total += gametefreq[genome][i]
        outHfreq = []
        for x in haplotypes:
            try:
                outHfreq.append(str(round(float(gametefreq[genome][x])/float(total)*100, 2)))
            except KeyError:
                outHfreq.append('0')
        print(sampleID+"\t"+"\t".join(outHfreq)+"\t"+str(total))
