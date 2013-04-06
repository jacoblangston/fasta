# Author: Kelsey and Jacob Langston
# Url: https://github.com/jacoblangston/fasta
#
# Generates fasta and qual files using the format specified at:
# http://en.wikipedia.org/wiki/FASTA_format
# http://bioperl.org/wiki/Qual_sequence_format
#
# How to convert fasta to fastq using BioPython:
# http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec:SeqIO-fastq-conversion

import sys
from optparse import OptionParser
from os.path import basename
from FastaGenerator import *
from ReadsGenerator import *

def addOptions(parser):
    parser.add_option("-n", "--number", dest="bases", 
                      help = "The number of base pairs")
    parser.add_option("-o", "--output", dest="output",
                      help = "The name of the output file")
    parser.add_option("-g", "--gcContentPercentage", dest="gc",
                      help = "The percentage of the simulated reference genome that consists of guanines and cytosines")
    parser.add_option("-b", "--repeatBase", dest="repeatBase",
                      help = "The base to repeat in a repeat region containing a single base")
    parser.add_option("-R", "--repeatRegionCount", dest="repeatRegionCount",
                      help = "The number of times to repeat a single base region")
    parser.add_option("-L", "--repeatBaseLength", dest="repeatBaseLength",
                      help = "The length of a single base region")
    parser.add_option("-s", "--segment", dest="segment",
                      help = "A specific region to repeat")
    parser.add_option("-t", "--segmentCount", dest="segmentCount",
                      help = "The number of times to repeat a segment")
    parser.add_option("-H", "--heterozygosity", dest="heterozygosity",
                      help = "Percent of heterozygosity")
    parser.add_option("-i", "--insertionPercent", dest="insertionPercent",
                      help = "The percent of insertions")
    parser.add_option("-d", "--deletionPercent", dest="deletionPercent",
                      help = "The percent of deletions")
    parser.add_option("-C", "--convertFormat", dest="convertFormat",
                      help = "Format to convert fasta file to; clustal, fastq-illumina, fastq-sanger, fastq-solexa, phylip, phd, tab, stockholm")
    parser.add_option("-r", "--generateReads", dest="generateReads",
                      help = "Specify the method to use for read generation; 454, exact, illumina-single, or illumina-paired")
    parser.add_option("-l", "--readLength", dest="readLength",
                      help = "The length of the read")
    parser.add_option("-c", "--readCoverage", dest="readCoverage",
                      help = "Specify the amount of coverage")
    parser.add_option("-f", "--referenceFile", dest="referenceFile",
                      help = "The path of an existing reference file in FastA format.  Using this option will not generate a FastA reference file.")
                
if __name__ == '__main__':
    parser = OptionParser()
    addOptions(parser)
    (options, args) = parser.parse_args()
    optionsDictionary = vars(options)
    
    baseList = ['G', 'C', 'T', 'A']
    numberOfBases = 1000
    readLength = 0
    readCoverage = 0
    heterozygosity = 0
    filename = None
    baseFilename = None
    generateReads = None
    
    if 'bases' in optionsDictionary and optionsDictionary['bases'] is not None:
        numberOfBases = int(optionsDictionary['bases'])    
    
    if 'output' in optionsDictionary and optionsDictionary['output'] is not None:
        filename = optionsDictionary['output']
    else:
        filename = 'originalsequence-' + str(numberOfBases) + '.fasta'
        
    if 'referenceFile' in optionsDictionary and optionsDictionary['referenceFile'] is not None:
        filename = optionsDictionary['referenceFile']
    
    baseFilename = basename(filename).rsplit(".")[0]  
    
    if 'heterozygosity' in optionsDictionary and optionsDictionary['heterozygosity'] is not None:
        try:
            heterozygosity = float(optionsDictionary['heterozygosity'])
        except ValueError:
            print "Heterozygosity must be numeric."    
        
    if 'generateReads' in optionsDictionary and optionsDictionary['generateReads'] is not None:
        generateReads = optionsDictionary['generateReads']
        
        if generateReads != "454" and generateReads != "exact" and generateReads != "illumina-single" and generateReads != "illumina-paired":
            print "Only 454, exact, illumina-single, and illumina-paired methods are available to generate reads."
            sys.exit()
            
        if generateReads == "illumina-single":
            readLength = 150
        elif generateReads == "illumina-paired":
            readLength = 100
        elif generateReads == "454":
            readLength = 400
        
    if 'readLength' in optionsDictionary and optionsDictionary['readLength'] is not None:
        try:
            readLength = int(optionsDictionary['readLength'])
        except ValueError:
            print "Read length must be an int, i.e. 1000."
            sys.exit()
        
    if 'readCoverage' in optionsDictionary and optionsDictionary['readCoverage'] is not None:
        try:
            readCoverage = int(optionsDictionary['readCoverage'])
            
            if heterozygosity > 0:
                readCoverage /= 2
        except ValueError:
            print "Read coverage must be an int, i.e. 25."
            sys.exit()
        
    if generateReads is not None and readCoverage == 0:
        print "Read coverage must be specified to generate reads."
        sys.exit()
        
    if generateReads is not None and generateReads == "exact" and readLength == 0:
        print "Read length must be specified to generate reads for exact.  There is no default value."
        sys.exit()
    
    if 'referenceFile' in optionsDictionary and optionsDictionary['referenceFile'] is None:
        f = FastaGenerator(baseList, optionsDictionary, filename, baseFilename, numberOfBases, heterozygosity)
    
    if generateReads is not None:
        if generateReads == "454":
            g = FourFiveFour(baseList, filename, baseFilename, readCoverage, readLength)
        elif generateReads == "illumina-single":
            g = IlluminaSingle(baseList, filename, baseFilename, readCoverage, readLength)
        elif generateReads == "illumina-paired":
            g = IlluminaPaired(baseList, filename, baseFilename, readCoverage, readLength)
        elif generateReads == "exact":
            g = Exact(baseList, filename, baseFilename, readCoverage, readLength)
        
    sys.exit()