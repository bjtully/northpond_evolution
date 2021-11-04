#!/usr/bin/env python3

'''
Make plots of dN/dS vs dS for all MAGs
All comparisons to the reference sequences - designated as [1]
in the seqIDs
Exclude all comparisons that have dN = 0.0000 & dS = 0.0000

'''

import glob
import plotly.graph_objects as go
from statistics import mean

# ingenomes = [
#     'NORP83',
#     'NORP139',
#     'NORP147',
#     'NORP163',
#     'NORP169',
#     'NORP246',
#     'NORP6',
#     'NORP57',
#     'NORP100',
#     'NORP167',
# ]

ingenomes = ['NORP83']
dndsvalues = {}
for genome in ingenomes:
    dndsvalues[genome] = {}
    genelist = glob.glob("./dNdS-"+genome+"/*GENECALL*") #Read in directories for each gene as created by Gretel
    combineddnds = []
    combinedds = []
    combinednames = []
    for genedirectory in genelist:
        genename = genedirectory.split("/")[-1].split("-")[1] + "_" + genedirectory.split("/")[-1].split("-")[-1]
        dndsvalues[genome][genename] = {'dnds':[], 'ds':[], 'pair':[]}
        datarows = False
        for line in open(genedirectory+"/"+"rst", "r"): #Read the PAML/CODEML produced rst file and 
            if datarows == True:
                data = line.strip().split()
                if len(data) != 0:
                    if int(data[0]) == 1 or int(data[1]) == 1: #Excluding dN/dS results that are poor quality
                        if float(data[6]) != 99.0000: #Indicated infinite float calc
                            if float(data[4]) != 0.0000 and float(data[5]) != 0.0000: #Excludes dN == 0 and dS == 0
                                dndsvalues[genome][genename]['dnds'].append(float(data[6]))
                                combineddnds.append(float(data[6]))
                                dndsvalues[genome][genename]['ds'].append(float(data[5]))
                                combinedds.append(float(data[5]))
                                dndsvalues[genome][genename]['pair'].append(genename+"_"+data[0]+"vs"+genename+"_"+data[1])
                                combinednames.append(genename+"_"+data[0]+"vs"+genename+"_"+data[1])
            if line[:3] == "seq":
                datarows = True
            else:
                continue
    meandnds = []
    meands = []
    meanname = []
    for gene in dndsvalues[genome]: #Calculate a mean dN/dS value for each gene
        if len(dndsvalues[genome][gene]['ds']) > 0:
            if mean(dndsvalues[genome][gene]['ds']) > 0.01: #Exclude outputs with dS <= 0.01
                meanname.append(gene)
                meands.append(mean(dndsvalues[genome][gene]['ds']))
                meandnds.append(mean(dndsvalues[genome][gene]['dnds']))

    outfile = open(genome+".dNdS.geneaverage.tsv", "w")
    outfile.write("Gene ID\tMean dN/dS\tMean dS\n")
    if len(dndsvalues[genome].keys()) > 0:
        fig = go.Figure(data=go.Scatter(
                x = meands,
                #x = combinedds,
                y = meandnds,
                #y = combineddnds,
                mode = 'markers',
                text = meanname,
                #text = combinednames,
                marker_color='rgba(152, 0, 0, .8)'
                ))
        fig.update_layout(title=genome)
        fig.update_xaxes(title_text='dS')
        fig.update_yaxes(title_text='dN/dS')
        fig.update_traces(mode='markers', marker_line_width=2)
        fig.update_xaxes(range=[0, 1])
        fig.update_yaxes(range=[0, 5])
        fig.write_image(genome+".geneaverage.png")
        fig.write_html(genome+".geneaverage.html")
        for i in range(len(meanname)):
            if meands[i] < 1 and meandnds[i] < 5:
                outfile.write(meanname[i]+"\t"+str(round(meandnds[i],4))+"\t"+str(round(meands[i],4))+"\n")

