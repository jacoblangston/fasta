# Author: Kelsey and Jacob Langston
# Url: https://github.com/jacoblangston/fasta
#
# Generates fasta files using the format specified at
# http://en.wikipedia.org/wiki/FASTA_format.

import random
import sys
from optparse import OptionParser

def addOptions(parser):
    parser.add_option("-n", "--number", dest="bases", 
                      help = "the number of base pairs")
    parser.add_option("-o", "--output", dest="output",
                      help = "the name of the output file")
    parser.add_option("-g", "--gcContentPercentage", dest="gc",
                      help = "the percentage of the simulated reference genome that consists of guanines and cytosines")
    parser.add_option("-b", "--repeatBase", dest="repeatBase",
                      help = "the base to repeat in a repeat region containing a single base")
    parser.add_option("-r", "--repeatRegionCount", dest="repeatRegionCount",
                      help = "the number of times to repeat a single base region")
    parser.add_option("-l", "--repeatBaseLength", dest="repeatBaseLength",
                      help = "the length of a single base region")
    parser.add_option("-s", "--segment", dest="segment",
                      help = "a specific region to repeat")
    parser.add_option("-t", "--segmentCount", dest="segmentCount",
                       help = "the number of times to repeat a segment")    

if __name__ == '__main__':
    parser = OptionParser()
    addOptions(parser)
    (options, args) = parser.parse_args()
    optionsDictionary = vars(options)
    
    baseList = ['G', 'C', 'T', 'A']
    filename = None
    numberOfBases = 1027
    gcContentPercentage = .41
    repeatBase = []
    repeatRegionCount = []
    segmentCount= []
    repeatBaseLength = []
    numberOfRepeatCharacters = 0
    segment = []
    repeatSegment = ''
    
    if 'gc' in optionsDictionary and optionsDictionary['gc'] is not None:
        gcContentPercentage = float(optionsDictionary['gc'])
        
    if 'bases' in optionsDictionary and optionsDictionary['bases'] is not None:
        numberOfBases = int(optionsDictionary['bases'])
            
    if 'output' in optionsDictionary and optionsDictionary['output'] is not None:
        filename = optionsDictionary['output']
    else:
        filename = 'originalsequence-' + str(numberOfBases) + '.fa'
    
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
    
    if len(repeatBase) != len(repeatRegionCount) or len(repeatBase) != len(repeatBaseLength) or len(repeatRegionCount) != len(repeatBaseLength):
        print "Repeat bases, repeat region count, and repeat base length must have the same number of items specified."
        sys.exit()
        
    if len(segment) != len(segmentCount):
        print "Segment and segment count must have the same number of items specified."
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

    output = '>' + filename + '\n'
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
    j = 0
    while j < len(bases):
        sortedBases.append(bases[j])
        j += 1
        if j > 0 and j % 79 == 0:
            sortedBases.append("\n")
            print j
        
    output += ''.join(sortedBases)
    output += '\n'
    
    file = open(filename, 'w')
    file.write(output)
    
    print "Wrote " + filename
    sys.exit()