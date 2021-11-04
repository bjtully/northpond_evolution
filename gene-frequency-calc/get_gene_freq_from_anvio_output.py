#just makes a heatmap of all gene frequencies over time. have already done this in excel but want a prettier one.
#usage: gets gene frequencies from the anvi'o data (very simple script but i want to look at it)

#import statistics
import numpy as np
import pandas as pd
import sys

infile = sys.argv[1]

outfile = open(str(infile).replace('.txt', '-genefreq-script.txt'), 'w')

#read in the data file
data = pd.read_table(infile, index_col='key') 

#get everything except the key column
#nokey_data = data.loc[:, data.columns != 'key']

#get the median of each column
mediandata = data.median()
#print mediandata

#get the median of a specific column
#data.loc[:,"CTD_T0_SORTED_FILTERED"].median()

genefreq_data = (data)/(mediandata)
#print genefreq_data

genefreq_data.to_csv(outfile, index=True, sep='\t')