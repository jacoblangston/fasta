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
from fastagenerator import *
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
    parser.add_option("-r", "--errorModel", dest="errorModel",
                      help = "Specify the error model to use for read generation; 454, exact, illumina-single, or illumina-paired")
    parser.add_option("-l", "--readLength", dest="readLength",
                      help = "The length of the read")
    parser.add_option("-c", "--readCoverage", dest="readCoverage",
                      help = "Specify the amount of coverage")
    parser.add_option("-f", "--referenceFile", dest="referenceFile",
                      help = "The path of an existing reference file in FastA format.  Using this option will not generate a FastA reference file.")
                
if __name__ == '__main__':
    parser = OptionParser() # Creates the parser object
    addOptions(parser) 
    (options, args) = parser.parse_args() # Parses the arguments that the user enters
    optionsDictionary = vars(options) # Assigns the arguments to the parser destination
    
    baseList = ['G', 'C', 'T', 'A'] # Array of characters
    numberOfBases = 1000
    readLength = 0
    readCoverage = 0
    heterozygosity = 0
    filename = None
    baseFilename = None # Filename without the extension
    errorModel = None
    
    """
    If the key exists in the optionsDictionary and the user has entered in a value, 
    then assigns the entered value to its appropriate key in the dictionary.
    
    The key is the parser dest for an option.
    """
    if 'bases' in optionsDictionary and optionsDictionary['bases'] is not None:
        try:
            numberOfBases = int(optionsDictionary['bases'])
        except ValueError:
            print "Number of bases must be a number."
            sys.exit()        
    
    """
    If the key "output" is in the optionsDictionary and the user has entered in a value for it,
    then that value will be assigned to the filename variable.
    
    If the key "output" does not exist in the optionsDictionary or if the user does not enter
    in a value for the key "output", then the filename variable will be 
    "originalsequence-<numberOfBases>.fasta".
    """
    if 'output' in optionsDictionary and optionsDictionary['output'] is not None:
        filename = optionsDictionary['output']
    else:
        filename = 'originalsequence-' + str(numberOfBases) + '.fasta'
        
    if 'referenceFile' in optionsDictionary and optionsDictionary['referenceFile'] is not None:
        filename = optionsDictionary['referenceFile']
    
    """
    Getting the file name without the extention. Defines the variable baseFilename.
    Takes the filename, searches from right to left for a period, splits at each
    period it finds, and selects the filename. Only one period is assumed.
    """
    baseFilename = basename(filename).rsplit(".")[0]  
    
    if 'heterozygosity' in optionsDictionary and optionsDictionary['heterozygosity'] is not None:
        try:
            heterozygosity = float(optionsDictionary['heterozygosity'])
        except ValueError:
            print "Heterozygosity must be numeric."
            sys.exit()
        
    if 'errorModel' in optionsDictionary and optionsDictionary['errorModel'] is not None:
        errorModel = optionsDictionary['errorModel']
        
        if errorModel != "454" and errorModel != "exact" and errorModel != "illumina-single" and errorModel != "illumina-paired":
            print "Only 454, exact, illumina-single, and illumina-paired error models are available to generate reads."
            sys.exit()
            
        if errorModel == "illumina-single":
            readLength = 150
        elif errorModel == "illumina-paired":
            readLength = 100
        elif errorModel == "454":
            readLength = 400
        
    if 'readLength' in optionsDictionary and optionsDictionary['readLength'] is not None:
        try:
            readLength = int(optionsDictionary['readLength'])
        except ValueError:
            print "Read length must be a number, i.e. 1000."
            sys.exit()
        
    if 'readCoverage' in optionsDictionary and optionsDictionary['readCoverage'] is not None:
        try:
            readCoverage = int(optionsDictionary['readCoverage'])
            
            if heterozygosity > 0:
                readCoverage /= 2
        except ValueError:
            print "Read coverage must be a number, i.e. 25."
            sys.exit()
        
    if errorModel is not None and readCoverage <= 0:
        print "Read coverage must be specified to generate reads."
        sys.exit()
        
    if errorModel is not None and errorModel == "exact" and readLength <= 0:
        print "Read length must be specified to generate reads for exact.  There is no default value."
        sys.exit()
    
    """
    Generates the fasta file if the option -f is not specified.
    """
    if 'referenceFile' in optionsDictionary and optionsDictionary['referenceFile'] is None:
        FastaGenerator(baseList, optionsDictionary, filename, baseFilename, numberOfBases, heterozygosity)
    
    """
    If the user provides a value for the variable errorModel and that value is "454",
    then the program will generate a 454 reads file. The listed variables are
    passed in as parameters to the constructor of the class.
    """
    if errorModel is not None:
        if errorModel == "454":
            FourFiveFour(baseList, filename, baseFilename, readCoverage, readLength)
        elif errorModel == "illumina-single":
            IlluminaSingle(baseList, filename, baseFilename, readCoverage, readLength)
        elif errorModel == "illumina-paired":
            IlluminaPaired(baseList, filename, baseFilename, readCoverage, readLength)
        elif errorModel == "exact":
            Exact(baseList, filename, baseFilename, readCoverage, readLength)
        
    sys.exit()
