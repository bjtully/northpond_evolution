# running DESMAN in triplicate
# usage: run_DESMAN_triplicate.py [NORP number, ie NORP83]
# run it within the folder for that specific MAG or all hell will break loose

import sys
import glob
import os
#path = r'/workspace/data/randerson/NorthPond'

NORPnum = sys.argv[1]

replist = [1,2,3]
glist = [2,3,4,5,6,7,8]

for rep in replist: #iterate through replicates
	for g in glist: #iterate through haplotype numbers
		cline = '~/STRONG/STRONG/DESMAN/build/scripts-3.6/desman ' + str(NORPnum) + '_outsel_var.csv -g ' + str(g) + ' -e ' + str(NORPnum) + '_outtran_df.csv -o ' + str(NORPnum) + '_g' + str(g) + '_' + str(rep)
		print cline
		run_program = os.popen(cline)
		status = run_program.close()