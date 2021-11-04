#identifies single-copy COGs from the Alneberg list that are present in a MAG based on anvio-provided COG annotations
#then goes to anvio-generated snv files and only pulls out variants that are from those genes.
#usage: get_variants_SCGs.py [COG_function file] [snv file]

import pandas as pd
import sys
COG_function_file = sys.argv[1]
SNV_file = sys.argv[2]

#making a list of the samples we care about for later
samplelist = ['s82A_T0_SORTED_FILTERED','s82A_T1_SORTED_FILTERED','s82A_T2_SORTED_FILTERED','s82A_T3_SORTED_FILTERED','s82A_T4_SORTED_FILTERED','s82A_T5_SORTED_FILTERED','s82A_T6_SORTED_FILTERED','s82A_T7_SORTED_FILTERED','s82A_T8_SORTED_FILTERED']

#first, we make a list of the single copy COGs.
SCGlist = []
alnebergfile = open('alneberg_single_copy_COGs.txt', 'r')
for line in alnebergfile:
	cols = line.split('\t')
	SCGlist.append(cols[0])
#print SCGlist

#then, open the COG file and ID the genes in this specific MAG that are SCGs and put them into a list
NORP_SCGlist = []
COGfile = open(COG_function_file, 'r')
for cogline in COGfile:
	columns = cogline.split('\t')
	genecaller_id = columns[0]
	COGnum = columns[2]
	if COGnum in SCGlist:
		NORP_SCGlist.append(genecaller_id)
#print NORP_SCGlist

#then, we go into the SNV file and only write the genes that are relevant, from the samples that we care about.
#note that i have to use pandas because the snv outfiles change column order WHICH IS STUPID
SNV_infile = open(SNV_file, 'r')
data = pd.read_table(SNV_infile)
filtered_data = data[data['corresponding_gene_call'].isin(NORP_SCGlist)] #only pulls out the data for the gene calls i'm interested in
#print filtered_data
filtered_data = filtered_data[filtered_data['sample_id'].isin(samplelist)] #only pulls out the samples we're interested in (note that this is all of the time points, not just the relevant ones for the MAG)
#print filtered_data

#ok, now i'm going to fix the stupidness of anvi'o and re-order the columns the way i want so i don't have to use pandas, ha.
col_order = ['entry_id','unique_pos_identifier','pos','pos_in_contig','sample_id','corresponding_gene_call','in_partial_gene_call','in_complete_gene_call','base_pos_in_codon','codon_order_in_gene','codon_number','gene_length','reference','consensus','competing_nts','departure_from_reference','departure_from_consensus','n2n1ratio','entropy','coverage','cov_outlier_in_split','cov_outlier_in_contig','A','C','G','N','T']
filtered_data = filtered_data[col_order]

filtered_data.to_csv(str(SNV_file).replace('.txt', '_SCGonly.txt'), index=False)