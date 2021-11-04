# compiles the posterior distribution from the "fit" output in DESMAN to see what the values look like over several replicates.
# usage: compile_posterior_dist_DESMAN.py [NOPR number, i.e. NORP83]
# use it inside the folder for that specific MAG

import sys, os, numpy, glob
path = r'./'
NORPnum = sys.argv[1]

outfile = open(str(NORPnum) + '_posterior_mean_dev.txt', 'w')
outfile.write('num_haplotypes' + '\t' + 'replicate_num' + '\t' + 'posterior_mean_deviance' + '\n')

#this step just makes a list of all the clusters so you can iterate through the directories
L = []
for directory in glob.glob(os.path.join(path, str(NORPnum) + '_g*')):
	#print filename
	L.append(directory)
print L

#this step iterates through each directory name and does stuff with the files within
for i in L:
	directory_name = str(i).split('_')
	haplotypes = directory_name[1].replace('g', '')
	try:
		replicate = directory_name[2]
	except IndexError:
		replicate = str(1)
	outfile.write(str(haplotypes) + '\t' + str(replicate) + '\t')
	#print i
	for root, dir, files in os.walk(i, 'fit.txt'): #walks through the files, only the ones matching the fit file
		#print files	
		for file in files:
			if file == 'fit.txt':
				infile = open(str(i) + '/' + str(file), 'r')
				for line in infile:
					columns = line.split(',')
					postdist = columns[4]
					#print postdist
					outfile.write(str(postdist))
				infile.close()
outfile.close()