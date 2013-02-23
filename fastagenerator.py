# Author: Kelsey Langston
# Date: 1/28/2013
#
# Generates fasta files using the format specified at
# http://en.wikipedia.org/wiki/FASTA_format.

import random
import sys
from optparse import OptionParser

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-n", "--number", dest="bases", 
                      help = "the number of base pairs")
    parser.add_option("-o", "--output", dest="output",
                      help = "the name of the output file")
    
    (options, args) = parser.parse_args()
    optionsDictionary = vars(options)
    
    list = ['G', 'T', 'C', 'A']
    filename = 'originalSequence.txt'
    numberOfBases = 0
    
    if 'bases' in optionsDictionary and optionsDictionary['bases'] is not None:        
        numberOfBases = int(optionsDictionary['bases'])
    else:
        print "Number of bases must be specified.  Type --help for help."
        sys.exit(1)
            
    if 'output' in optionsDictionary:
        filename = optionsDictionary['output']
    else:
        filename = 'originalsequence-' + str(numberOfBases) + '.fa'
        
    output = '>' + filename + '\n'
    
    a = 0
    while a < numberOfBases:
        output += list[random.randint(0,3)]
        a += 1
        
        if a % 79 == 0:
            output += "\n"    
    
    file = open(filename, 'w')
    file.write(output)
    sys.exit()