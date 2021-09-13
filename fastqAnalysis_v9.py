# Requires: Bowtie2 (cmd line), Samtools (cmd line), Picard (scripts folder), and featureCounts (scripts folder)
# Execute with Python3 on Linux - using Python2 will yield an error

import os
import subprocess
############################################################################################################
def get_samples_from_fastq(fq_folder_name):
	'''	  
	  - All fastq file pairs should be formatted like "samplename_R1.fq" and "samplename_R2.fq"
	  - Collects all file names in fastq folder. 
	  - Identifies name of the sample from each fastq file pairs
	  - Returns list of sample names. These will be used for mapping.
	'''
	
	# collect sample names from the fastq file pairs in a folder
	list_of_names = []
	try:
		fnames = os.listdir(fq_folder_name)
		for i in fnames:
			i1 = i.split(".")[0].replace("_R1","")		# takes everything before '.' and removes _R1 or _R2
			i2 = i1.replace("_R2","")
			if i2 not in list_of_names:
				list_of_names.append(i2)

		print("Processing the following samples:")
		for z in list_of_names:
			print(z)
		print("")
		return list_of_names
	except FileNotFoundError:
		print("Error: No folder found named 'fastq/'. Please check the folder names and try again.")
		print()
		exit()

##########################################################################################################
def fastqcheck(samplename, fqlocation):
	"""
	  - Checks that the fastq filenames are formatted as expected (X_R1.fq, and X_R2.fq) 
	  - Confirms that both forward and reverse files are there. 
	  - If files are not formatted correctly, script exits.
	  - If successful returns None.
	"""

	f1 = samplename + "_R1.fq"
	f2 = samplename + "_R2.fq"	
	filenames1 = os.listdir(fqlocation)

	if f1 in filenames1:
		if f2 in filenames1:
			return f1, f2
		else:
			print("Fastq files are not formatted properly for sample " + samplename + ".")
			print("Exiting.")
			exit()

##########################################################################################################
def file_and_folder_check(rc_fldr, mapf, saf_file, gfldr, fasta_file, scripts_fldr, fq_fldr):
	'''
	  - Makes sure all of the expected folders and files are there before mapping
	  	- rc_fldr = readcounts output folder
	  	- gflder = genomes file folder
	  	- .saf file is a tab delimited gene coordinates file. Users can create their own using
	  	  existing .saf files by copy/pasting the file and replacing the data colums in Excel.
	'''
	# Checks on existence of folders/files/etc
	if not os.path.exists(mapf):		# make mapped folder
		os.makedirs(mapf)
	if not os.path.exists(rc_fldr):		# make readcounts folder
		os.mkdir(rc_fldr)
	if not os.path.exists(scripts_fldr):
		print("Script folder is missing. Replace and restart.")
		exit()
	if not os.path.exists(gfldr + saf_file):
		print("The .saf file is missing or improperly named. Add it to the genome files folder and re-try.")
		exit()
	if not os.path.exists(gfldr + fasta_file):
		print("The .fasta file is missing is improperly named. Add it to the genome files folder and re-try.")
		exit()

	# Collects systematic genome name from .fasta file. Returns bytes format. Decode to string.
	sn = subprocess.check_output(args=['head', '-n 1', gfldr + fasta_file]).strip()
	sysName = sn.decode().strip(">")
	
	# Gets sample names from fastq file pairs in the fastq folder
	fq_files = get_samples_from_fastq(fq_fldr)
	fq_file_names = []
	for name in fq_files:

		# Confirms expected formatting of fastq files from a given pair
		fq1, fq2 = fastqcheck(name, fq_fldr)
		
		if fq1 != "error":
			if fq2 != "error":
				fq_file_names.append([fq1,fq2,name.strip("_")])		# Stores names as [name_1.fq, name_2.fq, name]
			else:
				print("Warning! Cannot find fastq file #2: ", fq2)
				print("Check that the sample names are accurate, and that")
				print("files are have '..._R1.fq' or '..._R1.fastq'")
				exit()
		else:
			print("Cannot find fastq file #1: " + fq1)
			print("Check that the sample names are accurate, and that")
			print("files are have '..._R1.fq' or '..._R1.fastq'")
			exit()
	return fq_file_names, sysName

##########################################################################################################3
def process_file(file_pair, fq_folder, gnm_fldr, genome, fasta, mappedfolder, picard_use, saffile, rcf, sysGnm):
	"""
	  - Maps paired-end fastq files to the indicated genome (bottom) by issuing commands to the terminal.
	  - For each pair of paired-end fastq files in the fastq folder:
		- Maps both files to the user-specified genome (below) using Bowtie2
	  	- Converts the .sam file to .bam with Samtools
	  	- Makes a sorted .bam file with Samtools
	  	- If user requested it, Picard will remove PCR and optical duplicates > new sorted .bam file
	  	- Makes a wiggle file for visualization with Samtools
	  	- Makes a second wiggle file that is normalized to total counts
	  	- Uses featureCounts > read counts (quantitative table of feature:mapped read number for RNA and ChIP-Seq)
	"""

	fqf1 = file_pair[0]		# forward fastq file
	fqf2 = file_pair[1]		# reverse fastq file
	name = file_pair[2]		# generic file name

	print("Now processing: " + name)
	
	# Bowtie2 mapping parameters can be customized here
	c1 = "bowtie2 --no-mixed -x " + gnm_fldr + genome + " -1 " + fq_folder + fqf1 + " -2 " + fq_folder + fqf2 +" -S " + mappedfolder +  name + ".sam"		# no mixed option
	# The command issued to the command line (made above) is printed 
	print("CMD: ", c1)
	# The command is passed to the command line
	os.system(c1)
	# Samtools converts the .sam file to .bam
	c2 = 'samtools view -bS ' + mappedfolder + str(name) + ".sam > " + mappedfolder + str(name) + ".bam"
	print("CMD: ", c2)
	os.system(c2)
	# Samtools sorts the .bam file
	c3 = "samtools sort " + mappedfolder + str(name) + ".bam -o " + mappedfolder + str(name) + ".sorted.bam"
	print("CMD: ", c3)
	os.system(c3)

	# Picard use is selected through user input when map.py is executed
	if picard_use == "y":
		picard_use = "Y"
	if picard_use == "Y":
			os.system("java -jar scripts/picard.jar CollectAlignmentSummaryMetrics R=" + gnm_fldr + fasta + " I=" + name + ".sorted.bam O=" + mappedfolder + name + ".metrics")
			os.system("java -jar scripts/picard.jar MarkDuplicatesWithMateCigar REMOVE_DUPLICATES=true M=" + mappedfolder + name + ".metrics I=" + mappedfolder + name + ".sorted.bam O=" + mappedfolder + name + "_DR.sam") 
			os.system("samtools view -bS " + mappedfolder + name + "_DR.sam > " + mappedfolder + name + "_DR.bam")
			# Feature counts is called here
			C6a = "scripts/featureCounts -a " + gnm_fldr + saffile + " -F SAF -o " + rcf + "FC_" + name + ".txt " + mappedfolder + name + "_DR.bam"
			print("CMD: ", C6a)
			os.system(C6a)
			# The sorted .bam file is converted to an mpileup file (wig with some extra data columns)
			c3b = "samtools mpileup " + mappedfolder + name + "_DR.bam -o " + mappedfolder + name + ".mp"
			os.system(c3b)
	elif picard_use != "Y":
		# Feature counts is called here
		C6a = "scripts/featureCounts -a " + gnm_fldr + saffile + " -F SAF -o " + rcf + "FC_" + name + ".txt " + mappedfolder + name + ".sorted.bam"
		print("CMD: ", C6a)
		os.system(C6a)
		# The sorted .bam file is converted to an mpileup file (wig with some extra data columns)
		c3c = "samtools mpileup " + mappedfolder + name + ".sorted.bam -o " + mappedfolder + name + ".mp"
		print("CMD: ", c3c)
		os.system(c3c)

	# Linux "cut" command is used to select columns from the .mpileup file > wig
	c4 = "cut -f2,4 " + mappedfolder + name + ".mp > " + mappedfolder + name + ".cut"
	print(c4)
	print("CMD: ", c4)
	os.system(c4)

	# The normalize wig script scales the wig file.
	# Forumula is (value at each nucleotide / total counts)*10E6
	# This essentially normalizes to total reads
	c5 = "python scripts/normalize_wig_v3.py " +  mappedfolder + name + ".cut " + sysGnm
	print("CMD: ", c5)
	os.system(c5)	# output is .norm.wig
	print("Map.py: Processing complete: " + name)

##########################################################################################################
### User Customization
##########################################################################################################
# Genome name for mapping. This name needs to be recognizable by Bowtie2 meaning that the genome must
# have been pre-made using the bowtie2-build command and placed in the genomefiles folder. 
# For visualization of .wig files in MochiView:
#    - import the import the .fasta file into MochiView (Import: Sequence Set) 
#    - import the .mochi gene coordinates files (Import/Location Set (genes)/Format:MochiView)
g = "JH642" 						
##########################################################################################################
# Assumed file/folder names:
##########################################################################################################
ff = "fastq/" 						# Folder where trimmed/cleaned .fq files are
mf = "mapped/"						# Output mapping folder
rf = "read_counts/"					# Output folder for readcounts
gf = "genomefiles/"					# Input files for genomes
sf = "scripts/"						# Scripts/app folder
fasta_file = g + ".fasta"			# Name of the fasta formatted genome sequence file
saf_file = g + ".saf"				# Name of features for featureCounts
##########################################################################################################
# Execute:
##########################################################################################################
print()
print(" ********                     **                 **                         **                  **        ")
print("/**/////                     /**    ****        ****                       /**  **   **        //  ")
print("/**        ******    ****** ****** **//**      **//**   *******   ******   /** //** **   ****** **  ******")
print("/*******  //////**  **//// ///**/ /** /**     **  //** //**///** //////**  /**  //***   **//// /** **//// ")
print("/**////    ******* //*****   /**  //*****    ********** /**  /**  *******  /**   /**   //***** /**//***** ")
print("/**       **////**  /////**  /**   ////**   /**//////** /**  /** **////**  /**   **     /////**/** /////**")
print("/**      //******** ******   //**     /***  /**     /** ***  /**//******** ***  **      ****** /** ****** ")
print("//        //////// //////     //      ///   //      // ///   //  //////// ///  //      //////  // //////  ")
print("")
print(" *******  **                  ** **   ")
print("/**////**//  ******          /**// ")
print("/**   /** **/**///**  *****  /** ** *******   *****  ")
print("/******* /**/**  /** **///** /**/**//**///** **///**  ")
print("/**////  /**/****** /******* /**/** /**  /**/*******  ")
print("/**      /**/**///  /**////  /**/** /**  /**/**////  ")
print("/**      /**/**     //****** ***/** ***  /**//****** ")
print("//       // //       ////// /// // ///   //  //////  ")
print("")

try:
	answer = input("Should the fastq files be mapped to genome '" + g  + "'? (Y/N): ")
	if answer == "Y" or answer == "y":
		remove_duplicates = input("Do you want to remove PCR duplicates with Picard? (Y/N) ")
		
		# Collects sample names from fastq files and ensures all other files/folders are found
		# Collects systematic genome name from fasta file.
		fq_file_list, systematicName = file_and_folder_check(rf, mf, saf_file, gf, fasta_file, sf, ff)
		
		# Executes mapping program for each sample
		for fq_filepair in fq_file_list:
			process_file(fq_filepair, ff, gf, g, fasta_file, mf, remove_duplicates, saf_file, rf, systematicName)
		exit()
	elif answer == "n" or answer == "N" or answer == "No" or answer == "NO":
		print("User selected cancel. Exiting.")
	else:
		
		print("Wrong genome was indicated. Edit the script's user customization section to indicate the correct genome and re-start.")
		print("Exiting...")
		exit()
except NameError:
	print("Error: it appears you are using Python2. Use Python3.")
	print("Exiting.")
	print()
##########################################################################################################