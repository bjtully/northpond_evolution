#makes a bubble plot of RPKM and SNVs/kbp for only bins of interest in samples of interest

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys


#sample = sys.argv[1]

#read in the data file
data = pd.read_table('SNV-SAAV-bloomers-persisters-RPKMmin10-NEW-RPKM.txt') 

#only grab data points from the sample of interest
#filtered_data = data[(data['Bin'] == sample)]

#print filtered_data

#use seaborn to make the plot
sns.set_style("ticks") #no grids, just ticks
sns.despine(trim=True) #removes the top and side lines in the graph so it's just the x and y axis lines drawn
sns.set_context("paper", font_scale=0.90) #made it suitable for publishing a paper and also changes the font size
sns.set_color_codes()
current_palette = sns.color_palette('colorblind')
sns.palplot(current_palette)
fig = plt.figure(figsize=(10,6)) #i think this adjusts the size
ax = fig.add_subplot(1,1,1) #no idea what this does but i think it's related to the line above

# Set the scale of the x-and y-axes
#ax.set(yscale="log")
#ax.set(ylim=(10, 400))

#make a bubble plot, adjusts the marker size.
ax = sns.scatterplot(x='Sample', y='SNVs/kbp', hue='Bin', data=data, size='RPKM', sizes= (10,1000), edgecolor='black')


# Put the legend out of the figure
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

#change tick frequency
plt.xticks(np.arange(0, 10, 1.0))

ax.tick_params(labelsize=10) #makes the font smaller on the x axis
locs, labels = plt.xticks() #this and the next line make the labels rotate however many degrees
plt.setp(labels, rotation=00)
plt.savefig('bubble_plot_bins_of_interest-NEW-RPKM.pdf') 


