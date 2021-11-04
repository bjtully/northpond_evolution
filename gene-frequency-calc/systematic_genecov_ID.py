#just makes a heatmap of all gene frequencies over time. have already done this in excel but want a prettier one.
#usage: systematic_genecov_ID.py [genefreq file '*_gene_cov_detection-GENE-COVERAGES-genefreq-script.txt] [list of samples separated by commas, no spaces]

import numpy as np
import pandas as pd
import sys

infile = sys.argv[1]
samplelist = sys.argv[2]

outfile = open(str(infile).replace('-script.txt', '_max-min1.txt'), 'w')

sample_list = str(samplelist).split(',')
#sample_list.append('key')

#read in the data file
data = pd.read_table(infile, index_col = 'key') 
#print data

#only grab data points from the samples of interest
filtered_data = data[sample_list]
#print filtered_data

#iterate through rows and calculate max/min and print genes that satisfy a condition
for index, row in filtered_data.iterrows():
     #print(index,row)
     #print max(row)
     #print min(row)
     if max(row)-min(row) >= 1:
     	outfile.write(str(index))
     	outfile.write('\n')