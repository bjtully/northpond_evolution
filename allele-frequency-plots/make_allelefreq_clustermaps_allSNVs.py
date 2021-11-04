#just makes a clustered heatmap of all gene frequencies over time for a specific MAG. 
#usage: make_allelefreq_clustermap.py [data file] [samples of interest, separated by commas, no spaces]

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

#only grab data points from samples of interest
filtered_data = data.loc[data['sample_id'].isin(sample_list)]


#drop infinite values
fixeddata = filtered_data.replace([np.inf, -np.inf], np.nan).dropna(axis=1)

#pivot data for heatmap
heatmap_data = pd.pivot_table(fixeddata, values='allele_freq', index=['sample_id'], columns='unique_pos_identifier')

#change size
fig, ax = plt.subplots(figsize=(20,7))   

#make heatmap
ax = sns.clustermap(heatmap_data, figsize=(25,10), row_cluster=False, cmap="YlGnBu", vmin=0, vmax=1)

#save plot
plt.savefig(str(infile).replace('.txt', '_clustermap-allelefreq_ALL_SNVs.pdf')) 


