#this script takes the filtered SNV output and formats it for use in DESMAN
#usage: convert_anvio_to_DESMAN.py [snv.pos_entropy filtered file that i ran with get_variants_SCG.py]


###WHERE I'M AT ON THIS: JUST CHANGED THE INPUT FILE TO A STRING USING (.READ) AND IT SEEMS TO BE WORKING BETTER; BUT NOT QUITE THERE YET-- OUTPUT FILE LOOKS FUNKY-- KEEP THINKING ABOUT THIS.

import sys
snvfile = sys.argv[1]

infile = open(snvfile).read() #reads it as a string rather than as a file, might make this work better
outfile = open(str(snvfile).replace('.txt', '_DESMANfriendly.txt'), 'w')
outfile.write('Contig,Position,s82A_T0_SORTED_FILTERED-A,s82A_T0_SORTED_FILTERED-C,s82A_T0_SORTED_FILTERED-G,s82A_T0_SORTED_FILTERED-T,s82A_T1_SORTED_FILTERED-A,s82A_T1_SORTED_FILTERED-C,s82A_T1_SORTED_FILTERED-G,s82A_T1_SORTED_FILTERED-T,s82A_T2_SORTED_FILTERED-A,s82A_T2_SORTED_FILTERED-C,s82A_T2_SORTED_FILTERED-G,s82A_T2_SORTED_FILTERED-T,s82A_T3_SORTED_FILTERED-A,s82A_T3_SORTED_FILTERED-C,s82A_T3_SORTED_FILTERED-G,s82A_T3_SORTED_FILTERED-T,s82A_T4_SORTED_FILTERED-A,s82A_T4_SORTED_FILTERED-C,s82A_T4_SORTED_FILTERED-G,s82A_T4_SORTED_FILTERED-T,s82A_T5_SORTED_FILTERED-A,s82A_T5_SORTED_FILTERED-C,s82A_T5_SORTED_FILTERED-G,s82A_T5_SORTED_FILTERED-T,s82A_T6_SORTED_FILTERED-A,s82A_T6_SORTED_FILTERED-C,s82A_T6_SORTED_FILTERED-G,s82A_T6_SORTED_FILTERED-T,s82A_T7_SORTED_FILTERED-A,s82A_T7_SORTED_FILTERED-C,s82A_T7_SORTED_FILTERED-G,s82A_T7_SORTED_FILTERED-T,s82A_T8_SORTED_FILTERED-A,s82A_T8_SORTED_FILTERED-C,s82A_T8_SORTED_FILTERED-G,s82A_T8_SORTED_FILTERED-T' + '\n')

#first, iterate through file and make a gene list.
genelist = [] #added this so it would skip the first line
lines = infile.split('\n')
for line in lines[1:-1]: #ignore first and last lines
	#print line
	try:
		cols = line.split(',')
		genecall = cols[5]
		if genecall not in genelist:
			genelist.append(genecall)
	except IndexError:
		pass
#print genelist #ok! i made a gene list.

print genelist

#now iterate through the gene list, and iterate through the file and make a list of positions for each genecall
for item in genelist:
	#print 'now working on gene number ' + str(item)
	positionlist = [] #create a position list for this genecall
	for line2 in lines[1:-1]:
		#print line2
		cols2 = line2.split(',')
		genecall2 = cols2[5]
		#print genecall2
		position = cols2[3]
		if genecall2 == item:
			if position not in positionlist:
				positionlist.append(position)
	#print positionlist
		#except IndexError:
			#pass
			
	#ok, now we're still in the for gene in genelist loop, and we have an established position list. now iterate through the position list.
	for pos in positionlist:
		#print 'now working on position ' + str(pos)
		samplelist = [] #create a sample list for this position
		for line3 in lines[1:-1]:
	#	print line
			cols = line3.split(',')
			genecall = cols[5]
			position = cols[3]
			sampleid = cols[4]
			if genecall == item:
				if position == pos:
					samplelist.append(sampleid)
		#print 'sample list is '
		#print samplelist
	
#ok, now i've made all the lists, so i need to iterate through and see if it's the right gene call, right position, see if the sample is in there, otherwise, write 0,0,0,0
		outfile.write(str(item) + ','+ str(pos) + ',')
		#if this sample is present in the sample list, then go ahead and write output. do this for each sample.
		if 's82A_T0_SORTED_FILTERED' in samplelist:
			for line4 in lines[1:-1]:
				try:
					cols = line4.split(',')
					genecall = cols[5]
					position = cols[3]
					sampleid = cols[4]
					A = cols[22]
					C = cols[23]
					G = cols[24]
					N = cols[25]
					T = cols[26]
					if genecall == item:
						if position == pos:
							if sampleid == 's82A_T0_SORTED_FILTERED':
								outfile.write(str(A) + ',' + str(C) + ',' + str(G) + ',' + str(T) + ',')
				except IndexError:
					pass
		#otherwise, put in zeroes
		else:
			outfile.write('0,0,0,0,')
	
		#if this sample is present in the sample list, then go ahead and write output. do this for each sample.
		if 's82A_T1_SORTED_FILTERED' in samplelist:
			for line4 in lines[1:-1]:
				try:
					cols = line4.split(',')
					genecall = cols[5]
					position = cols[3]
					sampleid = cols[4]
					A = cols[22]
					C = cols[23]
					G = cols[24]
					N = cols[25]
					T = cols[26]
					if genecall == item:
						if position == pos:
							if sampleid == 's82A_T1_SORTED_FILTERED':
								outfile.write(str(A) + ',' + str(C) + ',' + str(G) + ',' + str(T) + ',')
				except IndexError:
					pass
		#otherwise, put in zeroes
		else:
			outfile.write('0,0,0,0,')
		if 's82A_T2_SORTED_FILTERED' in samplelist:
			for line5 in lines[1:-1]:
				try:
					cols = line5.split(',')
					genecall = cols[5]
					position = cols[3]
					sampleid = cols[4]
					A = cols[22]
					C = cols[23]
					G = cols[24]
					N = cols[25]
					T = cols[26]
					if genecall == item:
						if position == pos:
							if sampleid == 's82A_T2_SORTED_FILTERED':
								outfile.write(str(A) + ',' + str(C) + ',' + str(G) + ',' + str(T) + ',')
				except IndexError:
					pass
		#otherwise, put in zeroes
		else:
			outfile.write('0,0,0,0,')
		if 's82A_T3_SORTED_FILTERED' in samplelist:
			for line6 in lines[1:-1]:
				try:
					cols = line6.split(',')
					genecall = cols[5]
					position = cols[3]
					sampleid = cols[4]
					A = cols[22]
					C = cols[23]
					G = cols[24]
					N = cols[25]
					T = cols[26]
					if genecall == item:
						if position == pos:
							if sampleid == 's82A_T3_SORTED_FILTERED':
								outfile.write(str(A) + ',' + str(C) + ',' + str(G) + ',' + str(T) + ',')
				except IndexError:
					pass
		#otherwise, put in zeroes
		else:
			outfile.write('0,0,0,0,')			
		if 's82A_T4_SORTED_FILTERED' in samplelist:
			for line7 in lines[1:-1]:
				try:
					cols = line7.split(',')
					genecall = cols[5]
					position = cols[3]
					sampleid = cols[4]
					A = cols[22]
					C = cols[23]
					G = cols[24]
					N = cols[25]
					T = cols[26]
					if genecall == item:
						if position == pos:
							if sampleid == 's82A_T4_SORTED_FILTERED':
								outfile.write(str(A) + ',' + str(C) + ',' + str(G) + ',' + str(T) + ',')
				except IndexError:
					pass
		#otherwise, put in zeroes
		else:
			outfile.write('0,0,0,0,')		
		if 's82A_T5_SORTED_FILTERED' in samplelist:
			for line8 in lines[1:-1]:
				try:
					cols = line8.split(',')
					genecall = cols[5]
					position = cols[3]
					sampleid = cols[4]
					A = cols[22]
					C = cols[23]
					G = cols[24]
					N = cols[25]
					T = cols[26]
					if genecall == item:
						if position == pos:
							if sampleid == 's82A_T5_SORTED_FILTERED':
								outfile.write(str(A) + ',' + str(C) + ',' + str(G) + ',' + str(T) + ',')
				except IndexError:
					pass
		#otherwise, put in zeroes
		else:
			outfile.write('0,0,0,0,')
		if 's82A_T6_SORTED_FILTERED' in samplelist:
			for line9 in lines[1:-1]:
				try:
					cols = line9.split(',')
					genecall = cols[5]
					position = cols[3]
					sampleid = cols[4]
					A = cols[22]
					C = cols[23]
					G = cols[24]
					N = cols[25]
					T = cols[26]
					if genecall == item:
						if position == pos:
							if sampleid == 's82A_T6_SORTED_FILTERED':
								outfile.write(str(A) + ',' + str(C) + ',' + str(G) + ',' + str(T) + ',')
				except IndexError:
					pass
		#otherwise, put in zeroes
		else:
			outfile.write('0,0,0,0,')				
		if 's82A_T7_SORTED_FILTERED' in samplelist:
			for line10 in lines[1:-1]:
				try:
					cols = line10.split(',')
					genecall = cols[5]
					position = cols[3]
					sampleid = cols[4]
					A = cols[22]
					C = cols[23]
					G = cols[24]
					N = cols[25]
					T = cols[26]
					if genecall == item:
						if position == pos:
							if sampleid == 's82A_T7_SORTED_FILTERED':
								outfile.write(str(A) + ',' + str(C) + ',' + str(G) + ',' + str(T) + ',')
				except IndexError:
					pass
		#otherwise, put in zeroes
		else:
			outfile.write('0,0,0,0,')
		if 's82A_T8_SORTED_FILTERED' in samplelist:
			for line11 in lines[1:-1]:
				try:
					cols = line11.split(',')
					genecall = cols[5]
					position = cols[3]
					sampleid = cols[4]
					A = cols[22]
					C = cols[23]
					G = cols[24]
					N = cols[25]
					T = cols[26]
					if genecall == item:
						if position == pos:
							if sampleid == 's82A_T8_SORTED_FILTERED':
								outfile.write(str(A) + ',' + str(C) + ',' + str(G) + ',' + str(T))
				except IndexError:
					pass
		#otherwise, put in zeroes
		else:
			outfile.write('0,0,0,0')
		outfile.write('\n')
#
#infile.close()
outfile.close()
	