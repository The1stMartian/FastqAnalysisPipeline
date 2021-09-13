![](https://github.com/The1stMartian/FastqAnalysisPipeline/blob/main/doc/logo.png)

The fastqAnalysis script is a simple Seq data analysis pipeline written in Python3 for Linux. It will map fastq files to the user-specified genome using Bowtie2 and output both sorted .bam files and normalized wiggle files that can be visualized in MochiView. It will also create a read-counts file based upon the features listed in the accompanying .saf file. Read count data is particularly useful for ChIP-Seq and RNA-Seq analyses.

Pre-requisites:
- Python3 [Installation instructions](https://docs.anaconda.com/anaconda/install/linux/).
- Bowtie2. [Installation instructions](https://www.metagenomics.wiki/tools/bowtie2/install).
- SamTools [Installation instructions](https://bioinformaticsreview.com/20210404/installing-samtools-on-ubuntu/#:~:text=%20Installing%20SAMtools%20on%20Ubuntu%20%201%20Preparing,We%20are%20in%20the%20same%20directory...%20More%20).


Pipeline:
1) Users should normalize all fastq file names. The script will be expecting paired-end reads with format:<br><br>
           sampleName1_R1.fq<br>
           sampleName1_R2.fq<br>
           sampleName2_R1.fq<br>
           sampleName2_R2.fq<br><br>
           <i>Place all files in a folder named "fastq" in the same directory as the fastqAnalysis.py script.</i>

2) If needed, create a genome for your model organism:
      - In the command line (linux/linux subsystem) bowtie2-build <fasta_file_path> <chosen_genome_name>
      - Place the resulting output files into the folder named "genomefiles"
      - Create a .saf file (feature coordinates/straind) for your genome by opening an existing .saf file and saving it as <chosen_genome_name>.saf
      - To input features:<br>
             - Download a features file from NCBI genome (example using bacteria)<br>
               https://www.ncbi.nlm.nih.gov/genome/browse#!/overview/<br>
             - Enter your species name and click enter<br>
             - Click "prokaryotes", NOT the species name link<br>
             - Find your species and under the Accession tab, click the link formatted like "GCA_ ..."<br>
             - Click the link for the FTP Directory for Genbank on the right<br>
             - Click the link that ends with "_feature_table.txt.gz"<br>
             - After downloading, open the file in excel<br>
             - Copy/Paste the coordinates, systematic gene name, and strand information columns into the .saf file.<br>
             - Make sure the new .saf file is in the genomefiles folder<br>
      - Create .mochi formatted gene coordinates files for later use:<br>
             - Copy an existing .mochi file as <chosen_genome_name>.mochi<br>
             - Copy/Paste feature information from your features file into the .mochi file and save.<br>

3) Execute the mapping program:
      <i>python fastqAnalysis_v8.py</i>

      ![Cmd line](https://github.com/The1stMartian/fastqAnalysis/blob/main/doc/cmd.png)


4) Set up MochiView:
- Download from the Johnson Lab at UCSF [MochiView Website](http://www.johnsonlab.ucsf.edu/mochi/)
- Ensure you have an updated version of Java installed (Windows or Linux)
- Instructions for launching the app can be found [here](http://www.johnsonlab.ucsf.edu/mochiview-downloads)
             
Set up your genome:
- Import your .fasta file:

![](https://github.com/The1stMartian/fastqAnalysis/blob/main/doc/mochi1.png)

- Import gene coordinate files

![](https://github.com/The1stMartian/fastqAnalysis/blob/main/doc/mochi2.png)

- Import your normalized .wig files

![](https://github.com/The1stMartian/fastqAnalysis/blob/main/doc/mochi3.png)

- Plot the data (Click "New Plot")

![](https://github.com/The1stMartian/fastqAnalysis/blob/main/doc/mochi8.png)

- Select feature annotations to be shown

![](https://github.com/The1stMartian/fastqAnalysis/blob/main/doc/mochi5.png)

- Select the data set(s) you want to visualize and pick line or column displays (and color)

![](https://github.com/The1stMartian/fastqAnalysis/blob/main/doc/mochi6.png)

- Plot the data

![](https://github.com/The1stMartian/fastqAnalysis/blob/main/doc/mochi9.png)

- Results:

![](https://github.com/The1stMartian/fastqAnalysis/blob/main/doc/mochi7.png)


### Troubleshooting:
- One common error is that the .fasta file genome name doesn't match up the .wig file headers
- To see if this is causing a problem, look at the fasta file header using either "head -n 1 <fastaFilePath>" or simply by opening the file in notepad. The file is large so notepad can be slow.
- Write down the name in the first line, i.e. for JH642 the first line is >NZ_CP007800, so the systematic name is "NZ_CP007800".
- Now open the .wig files and look at the header. The name should match. For JH642 it looks like:

     track
     variableStep chrom=NZ_CP007800 span=1

     <i>If your genome's systematic name doesn't appear like "chrom=<systematicName>", then replace whatever is there with the proper name. No spaces!</i>

- Now try re-importing the .wig file.

### Notes on wiggle file normalization: 
The normalization script essentially allows all wiggle files to be directly compared by normalizing to the total reads. Without normalization, slight differences in the number of reads per library would throw off the scale. For example, if the same library was sequenced on two occasions, the resulting fastq files would likely have a different number of reads, resulting in one library simply having higher values than the other in visualization software. To account for differences in depth, the normalization script divides the value (read depth) at every nucleotide position by the sum of all values, then multiplies the resulting number (which is typically very small) by 10^6 to produce a values over 1 which is simply a convenient number to work with. The multiplier can be adjusted as need be by modifying the normalization script.  