Reference Simulator
=====

Current de novo assembly software encounters difficulty while attempting to assemble heterogeneous data. This tool can create simulated genomic data with increasing levels of heterogeneity. Developers that are currently constructing algorithms to handle such data will be able to test their de novo assemblers throughout production. Not only can the degree of heterogeneity, but also the number, length, and composition of repeat regions of the simulated data can be controlled by the user. In this way, assembly developers can focus on certain biological characteristics more exclusively, as well as the ability of their assemblers to handle data of increasing complexity.
 
The Reference Simulator simulates genomic data and outputs the information into a fasta file format. Only diploid genomes can be simulated. While the generated data will not perfectly emulate real-world data, the user will be able to control such biological attributes as GC content and interspersed repeats in an effort to mirror genuine genomic data. The output of the de novo assembler can be compared to the simulated genome reference and the robustness of an algorithm may be assessed.

Usage: RefSim.py [options]

Quick Start
========
* Download and install the lastest version of Python 2.x from http://www.python.org/download/
* Download and install the version of BioPython from http://biopython.org/wiki/Download that maps to the installed version of Python.  For instance if you
installed Python 2.7.x then you will need to download and install biopython-<version>.win32-py2.7.exe.
* Run git clone https://github.com/jacoblangston/fasta.git.

Options
========
  
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
                        
  -C CONVERTFORMAT      --convertFormat=CONVERTFORMAT
  			Format to convert fasta file to; clustal, fastq-illumina, fastq-sanger, fastq-solexa, phylip, phd, tab, or stockholm
		
  -r GENERATEREADS	--generateReads=GENERATEREADS
  			Specify the method to use for read generation; 454, exact, illumina-single, or illumina-paired
						
  -l READLENGTH		--readLength=READLENGTH
  			The length of the read
  									
  -c READCOVERAGE	--readCoverage=READCOVERAGE
			Specify the amount of coverage
			
  -f REFERENCEFILE	--referenceFile=REFERENCEFILE
  			The path of an existing reference file in FastA format.  Using this option will not generate a FastA reference file.
					
Usage
=====

Running the following will output a reference file in fasta format. The reference will be 100,000 bases in length where 47% will consist of G's and C's. GC content = (G's + C's)/Total Number of Bases.

python RefSim.py -n 100000 -o "output_file_name.fasta" -g 0.47 

Running the following will output a reference file that is 100,000 bases in length with the default GC content (41%) and four sets of user-specified repeat regions. There will be 45 repeat regions consisting of 6 A's in a row. There will be 7 repeat regions consisting of 4 G's in a row. There will be 15 repeat regions consisting of 3 T's in a row. Finally, there will be 2 repeat regions consisting of 75 N's in a row. The user can include as many single-base repeat regions as is possible given the specified length of the reference. Important note: un-specified single-base repeat regions will occur randomly throughout the reference due to chance.

python RefSim.py -n 100,000 -b "A,G,T,N" -R "45,7,15,2" -L "6,4,3,75" -o "output_file_name.fasta"

Running the following will output a reference file that is 100,000 bases in length. The reference will include the sequence "GAATTC" at three different, randomly chosen locations.

python RefSim.py -n 100000 -s "GAATTC" -t "3" -o "output_file_name.fasta"

Running the following will output a reference file that is 200,000 bases in length. The first 100,000 bases will be generated in the same fashion as if no heterozygosity had been selected. The second 100,000 bases will be an exact replica of the first 100,000 bases, but with the indicated alterations. The -H "0.09" option will result in 9% of the second 100,000 bases being randomly exchanged for an alternate A, T, G, or C. The -i "0.02" option will result in 2% of the second 100,000 bases consisting of randomly inserted bases. A frameshift for all remaining bases will result. The -d "0.06" option will result in 6% of the second 100,000 bases being randomly deleted, also resulting in frameshifts.

python RefSim.py -n 100000 -H "0.09" -i "0.02" -d "0.06" -o "output_file_name.fasta"

Running the following will output a reference file in fasta format that is 100,000 bases in length and that has all of the default parameters. A second file ending in ".fa" ("output_file_name.fa") will also be generated and will consist of simulated reads generated from the previously mentioned reference file. The reads will be simulated using the Roche/454 error model at a 30X coverage. The length of the reads will be 400 bases.

python RefSim.py -n 100000 -r "454" -l 400 -c 30 -o "output_file_name.fasta"

All of the discussed options can be run on the same line. 
