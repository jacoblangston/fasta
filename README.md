Reference Simulator
=====

Current de novo assembly software encounters difficulty while attempting to assemble heterogeneous data. This tool can create simulated genomic data with increasing levels of heterogeneity. Developers that are currently constructing algorithms to handle such data will be able to test their de novo assemblers throughout production. Not only can the degree of heterogeneity, but also the number, length, and composition of repeat regions of the simulated data can be controlled by the user. In this way, assembly developers can focus on certain biological characteristics more exclusively, as well as the ability of their assemblers to handle data of increasing complexity.
 
The Reference Simulator simulates genomic data and outputs the information into a fasta file format. Only diploid genomes can be simulated. While the generated data will not perfectly emulate real-world data, the user will be able to control such biological attributes as GC content and interspersed repeats in an effort to mirror genuine genomic data. The output of the de novo assembler can be compared to the simulated genome reference and the robustness of an algorithm may be assessed.

Usage: RefSim.py [options]

Options:
  
  -h                 	--help
						Show this help message and exit
 
  -n BASES           	--number=BASES
                        The number of base pairs
  
  -o OUTPUT          	--output=OUTPUT
                        The name of the output file
  
  -g GC              	--gcContentPercentage=GC
                        The percentage of the simulated reference genome that consists of guanines and cytosines

  -b REPEATBASE      	--repeatBase=REPEATBASE
                        The base to repeat in a repeat region containing a single base
  
  -R REPEATREGIONCOUNT  --repeatRegionCount=REPEATREGIONCOUNT
                        The number of times to repeat a single base region
  
  -L REPEATBASELENGTH 	--repeatBaseLength=REPEATBASELENGTH
                        The length of a single base region
  
  -s SEGMENT         	--segment=SEGMENT
                        A specific region to repeat
  
  -t SEGMENTCOUNT    	--segmentCount=SEGMENTCOUNT
                        The number of times to repeat a segment
  
  -H HETEROZYGOSITY  	--heterozygosity=HETEROZYGOSITY
                        Percent of heterozygosity
  
  -i INSERTIONPERCENT	--insertionPercent=INSERTIONPERCENT
                        The percent of insertions
  
  -d DELETIONPERCENT 	--deletionPercent=DELETIONPERCENT
                        The percent of deletions

  -C CONVERTFORMAT      --convertFormat=CONVERTFORMA
  			Format to convert fasta file to; clustal, fastq-illumina, fastq-sanger, fastq-solexa, phylip, phd, tab, or stockholm
						
  -r GENERATEREADS	--generateReads=GENERATEREADS
  			Specify the method to use for read generation; 454, exact, illumina-single, or illumina-paired
						
  -l READLENGTH		--readLength=READLENGTH
  			The length of the read
						
  -c READCOVERAGE	--readCoverage=READCOVERAGE
			Specify the amount of coverage
