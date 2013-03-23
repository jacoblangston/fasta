# Author: Kelsey and Jacob Langston
# Url: https://github.com/jacoblangston/fasta
#
# Generates fasta and qual files using the format specified at:
# http://en.wikipedia.org/wiki/FASTA_format
# http://bioperl.org/wiki/Qual_sequence_format
#
# How to convert fasta to fastq using BioPython:
# http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec:SeqIO-fastq-conversion

import random
import sys
from optparse import OptionParser
from os.path import basename
from Bio import SeqIO
from Bio.SeqIO.QualityIO import PairedFastaQualIterator
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
    
def changeBase(baseList, currentBase):
    newBase = currentBase
    while newBase == currentBase:
        index = random.randint(0, len(baseList) - 1)
        newBase = baseList[index]
        
    return newBase

def convert(convertFormat, filename, qualFilename):
    print "\nFastA file information:"
    for seq_record in SeqIO.parse(filename, "fasta"):
        print str(seq_record)  
    
    try:
        with open(baseFilename + ".fastq", "w") as q:
            records = PairedFastaQualIterator(open(filename), open(qualFilename))
            count = SeqIO.write(records, q, convertFormat)
            print "Converted %i records" % count
    except ValueError, e:
        print "Encountered error converting file: " + str(e)
    
def generateQual(baseFilename, numberOfBases, minQuality, maxQuality):
    qualFilename = baseFilename + ".qual"
    quality = ">" + baseFilename + "\n"
    
    i = 0
    while i < numberOfBases:
        quality += str(getQuality(minQuality, maxQuality)) + " "
        if i > 0 and i % 79 == 0:
            quality += "\n"
        i += 1
        
    with open(qualFilename, "w") as q:
        q.write(quality)
        
    return qualFilename
    
def getQuality(minQuality, maxQuality):
    return random.randint(minQuality, maxQuality)
            
if __name__ == '__main__':
    parser = OptionParser()
    addOptions(parser)
    (options, args) = parser.parse_args()
    optionsDictionary = vars(options)
    
    baseList = ['G', 'C', 'T', 'A']
    filename = None
    baseFilename = None
    numberOfBases = 1000
    gcContentPercentage = .41
    heterozygosity = 0
    insertionPercent = 0
    deletionPercent = 0
    numberOfRepeatCharacters = 0
    readLength = 0
    readCoverage = 0
    minQuality = 11
    maxQuality = 41
    repeatBase = []
    repeatRegionCount = []
    segmentCount= []
    repeatBaseLength = []
    segment = []
    repeatSegment = None
    convertFormat = None
    generateReads = None
    
    if 'gc' in optionsDictionary and optionsDictionary['gc'] is not None:
        gcContentPercentage = float(optionsDictionary['gc'])
        
    if 'bases' in optionsDictionary and optionsDictionary['bases'] is not None:
        numberOfBases = int(optionsDictionary['bases'])
            
    if 'output' in optionsDictionary and optionsDictionary['output'] is not None:
        filename = optionsDictionary['output']
    else:
        filename = 'originalsequence-' + str(numberOfBases) + '.fasta'
    
    baseFilename = basename(filename).rsplit(".")[0]  
    
    if 'repeatBase' in optionsDictionary and optionsDictionary['repeatBase'] is not None:
        repeatBase = optionsDictionary['repeatBase'].split(',')
    
    if 'repeatRegionCount' in optionsDictionary and optionsDictionary['repeatRegionCount'] is not None:
        repeatRegionCount = optionsDictionary['repeatRegionCount'].split(',')
    
    if 'repeatBaseLength' in optionsDictionary and optionsDictionary['repeatBaseLength'] is not None:
        repeatBaseLength = optionsDictionary['repeatBaseLength'].split(',')
        
    if 'segment' in optionsDictionary and optionsDictionary['segment'] is not None:
        segment = optionsDictionary['segment'].split(',')
        
    if 'segmentCount' in optionsDictionary and optionsDictionary['segmentCount'] is not None:
        segmentCount = optionsDictionary['segmentCount'].split(',')
        
    if 'heterozygosity' in optionsDictionary and optionsDictionary['heterozygosity'] is not None:
        heterozygosity = float(optionsDictionary['heterozygosity'])
        
    if 'insertionPercent' in optionsDictionary and optionsDictionary['insertionPercent'] is not None:
        insertionPercent = float(optionsDictionary['insertionPercent'])
        
    if 'deletionPercent' in optionsDictionary and optionsDictionary['deletionPercent'] is not None:
        deletionPercent = float(optionsDictionary['deletionPercent'])
        
    if 'convertFormat' in optionsDictionary and optionsDictionary['convertFormat'] is not None:
        convertFormat = optionsDictionary['convertFormat']
        
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
        readLength = int(optionsDictionary['readLength'])
        
    if 'readCoverage' in optionsDictionary and optionsDictionary['readCoverage'] is not None:
        readCoverage = int(optionsDictionary['readCoverage'])
    
    if len(repeatBase) != len(repeatRegionCount) or len(repeatBase) != len(repeatBaseLength) or len(repeatRegionCount) != len(repeatBaseLength):
        print "Repeat bases, repeat region count, and repeat base length must have the same number of items specified."
        sys.exit()
        
    if len(segment) != len(segmentCount):
        print "Segment and segment count must have the same number of items specified."
        sys.exit()
        
    if heterozygosity == 0 and (insertionPercent > 0 or deletionPercent > 0):
        print "Heterozygosity must be specified for insertions and deletions."
        sys.exit()

    i = 0
    baseCharacterCountTotal = 0
    while i < len(repeatBase):
        baseCharacterCountTotal += (len(repeatBase[i]) * int(repeatBaseLength[i]) * int(repeatRegionCount[i]))
        i += 1
        
    i = 0
    segmentCharacterCountTotal = 0
    while i < len(segment):
        segmentCharacterCountTotal += (len(segment[i]) * int(segmentCount[i]))
        i += 1
            
    numberOfRepeatCharacters = baseCharacterCountTotal + segmentCharacterCountTotal

    output = '>' + baseFilename + '\n'
    bases = []
    
    i = 0
    while i < (numberOfBases - numberOfRepeatCharacters):
        if i <= gcContentPercentage * numberOfBases:
            bases.append(baseList[random.randint(0,1)])
        else:
            bases.append(baseList[random.randint(2,3)])
        i += 1
        
    random.shuffle(bases)
    
    i = len(repeatRegionCount)
    while i > 0:
        k = 0
        while k < int(repeatRegionCount[i - 1]):
            index = random.randint(0, len(bases))
            l = 0
            while l < int(repeatBaseLength[i - 1]):
                bases.insert(index + l, repeatBase[i - 1].strip())
                l += 1
            k += 1
        i -= 1
    
    i = len(segmentCount)
    while i > 0:
        x = 0
        while x < int(segmentCount[i - 1]):
            index = random.randint(0, len(bases))
            l = 0
            while l < len(segment[i - 1]):
                bases.insert(index + l, segment[i - 1][l].strip())
                l += 1
            x += 1
        i -= 1
    
    sortedBases = []
    i = 0
    while i < len(bases):
        sortedBases.append(bases[i])
        i += 1
            
    if heterozygosity > 0:
        sortedBases += sortedBases
        
        i = 0
        changedBases = []
        numberOfBasesToChange = int((len(sortedBases) - 1) * heterozygosity)
        while i < numberOfBasesToChange:
            index = random.randint(numberOfBases, len(sortedBases) - 1)
            while changedBases.count(index) > 0:
                index = random.randint(numberOfBases, len(sortedBases) - 1)
            
            currentBase = sortedBases[index]
            sortedBases[index] = changeBase(baseList, currentBase)
            changedBases.append(index)
            i += 1
            
        i = 0
        numberOfDeletions = int((len(sortedBases) - 1) * deletionPercent)
        while i < numberOfDeletions:
            index = random.randint(numberOfBases, len(sortedBases) - 1)
            base = baseList[random.randint(0, len(baseList) - 1)]
            sortedBases.pop(index)
            sortedBases.append(base)
            i += 1        
        
        i = 0
        numberOfInsertions = int((len(sortedBases) - 1) * insertionPercent)
        while i < numberOfInsertions:
            index = random.randint(numberOfBases, len(sortedBases) - 1)
            base = baseList[random.randint(0, len(baseList) - 1)]
            sortedBases.insert(index, base)
            sortedBases.pop()
            i += 1
         
    totalBases = len(sortedBases)
    
    i = 0
    while i < len(sortedBases):
        if i > 0 and i % 79 == 0:
            sortedBases.insert(i - 1, "\n")
        i += 1
    
    output += ''.join(sortedBases)
    output += "\n"
    
    with open(filename, 'w') as fa:
        fa.write(output)
        print "Wrote " + filename
    
    qualFilename = generateQual(baseFilename, totalBases, minQuality, maxQuality)
    
    if convertFormat is not None:
        convert(convertFormat, filename, qualFilename)
        
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