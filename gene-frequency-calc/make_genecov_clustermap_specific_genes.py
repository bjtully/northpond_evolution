#just makes a clustered heatmap of all gene frequencies over time for a specific MAG. 
#usage: make_genecov_clustermap.py [data file] [samples of interest, separated by commas, no spaces]

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys

infile = sys.argv[1]
samplelist = sys.argv[2]
sample_list = str(samplelist).split(',')

#read in the data file
data = pd.read_table(infile, index_col = 0) 

#only grab data points from the samples of interest
filtered_data = data[sample_list]


#only grab genes of interest
binsplit = str(infile).split('_')
bin = binsplit[0]
genelist = []
genefile = open(str(bin) + '_gene_cov_detection-GENE-COVERAGES-genefreq_max-min1.txt', 'r')
for item in genefile:
	item = item.strip('\n')
	genelist.append(int(item))
genefile.close()
filtered_data2 = filtered_data.loc[genelist]

#drop infinite values
fixeddata = filtered_data2.replace([np.inf, -np.inf], np.nan).dropna(axis=1)

#convert this data (anvi'o output) into longform (for seaborn)
tidydata = pd.melt(fixeddata.reset_index(), id_vars='key', var_name='sample', value_name='gene_freq') 

#pivot data for heatmap
heatmap_data = pd.pivot_table(tidydata, values='gene_freq', index=['sample'], columns='key')

#make heatmap
ax = sns.clustermap(heatmap_data, figsize=(25,10), row_cluster=False, cmap="YlGnBu", vmin=0, vmax=5)

#save plot
plt.savefig(str(infile).replace('.txt', 'clustermap-specific_genes.pdf')) 


