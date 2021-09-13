####################################################################################
##### purpose is to normalize wiggle files to account for different number of 
##### reads used in a given data set. This works as of 11/18/2016
####################################################################################

import sys
####################################################################################
def get_total(wig_file):
	"""For a given wiggle file, add up the total number of counts"""
	
	line_num = 0
	total = 0

	for line in wig_file:
		splt = line.split()
		num = splt[1]
		number = float(num.strip("\n"))
		total += number

	return total

####################################################################################
def new_wig_file(num_counts, gnm):
	"""Using a 'cut' file as input, adjusts the value at each nucleotide position
	based upon its percentage of the total counts. Then multiplies by 10^6 which 
	is an arbitrary, but useful number. Wiggle files can then be compared to others
	even if they have differe numbers of reads used to produce them."""

	out = open((sys.argv[1][:-4] + ".norm.wig") , "w")
	out.write("track\nvariableStep chrom=" + gnm +" span=1\n")

	line_ctr = 0

	for line in open(sys.argv[1], "r"):
		splt = line.split()
		num = splt[1]
		number = float(num.strip("\n"))
		pos = splt[0]
		new_number = (number/num_counts)*10000000
		out.write(str(pos) + "\t" + str(new_number) + "\n")

	out.close()

####################################################################################
# Function calls
####################################################################################
f = open(sys.argv[1], "r")
gnm_name = sys.argv[2]

total_counts = get_total(f)
new_wig_file(total_counts, gnm_name)

sys.exit()