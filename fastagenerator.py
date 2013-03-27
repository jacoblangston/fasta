# Author: Kelsey and Jacob Langston
# Url: https://github.com/jacoblangston/fasta
import random
from optparse import OptionParser
from Bio import SeqIO
from Bio.SeqIO.QualityIO import PairedFastaQualIterator

class FastaGenerator(object):
    baseList = []
    filename = None
    baseFilename = None
    numberOfBases = 1000
    gcContentPercentage = .41
    heterozygosity = 0
    insertionPercent = 0
    deletionPercent = 0
    numberOfRepeatCharacters = 0
    minQuality = 11
    maxQuality = 41
    repeatBase = []
    repeatRegionCount = []
    segmentCount= []
    repeatBaseLength = []
    segment = []
    repeatSegment = None
    convertFormat = None    

    def __init__(self, baseList, optionsDictionary, filename, baseFilename, numberOfBases):
        self.baseList = baseList
        self.optionsDictionary = optionsDictionary
        self.filename = filename
        self.baseFilename = baseFilename
        self.numberOfBases = numberOfBases
        self.parseOptions()
        self.generateFasta()
    
    def parseOptions(self):
        if 'gc' in self.optionsDictionary and self.optionsDictionary['gc'] is not None:
            self.gcContentPercentage = float(self.optionsDictionary['gc'])
                
        if 'repeatBase' in self.optionsDictionary and self.optionsDictionary['repeatBase'] is not None:
            self.repeatBase = self.optionsDictionary['repeatBase'].split(',')
        
        if 'repeatRegionCount' in self.optionsDictionary and self.optionsDictionary['repeatRegionCount'] is not None:
            self.repeatRegionCount = self.optionsDictionary['repeatRegionCount'].split(',')
        
        if 'repeatBaseLength' in self.optionsDictionary and self.optionsDictionary['repeatBaseLength'] is not None:
            self.repeatBaseLength = self.optionsDictionary['repeatBaseLength'].split(',')
            
        if 'segment' in self.optionsDictionary and self.optionsDictionary['segment'] is not None:
            self.segment = self.optionsDictionary['segment'].split(',')
            
        if 'segmentCount' in self.optionsDictionary and self.optionsDictionary['segmentCount'] is not None:
            self.segmentCount = self.optionsDictionary['segmentCount'].split(',')
            
        if 'heterozygosity' in self.optionsDictionary and self.optionsDictionary['heterozygosity'] is not None:
            try:
                self.heterozygosity = float(self.optionsDictionary['heterozygosity'])
            except ValueError:
                print "Heterozygosity must be numeric."
            
        if 'insertionPercent' in self.optionsDictionary and self.optionsDictionary['insertionPercent'] is not None:
            try:
                self.insertionPercent = float(self.optionsDictionary['insertionPercent'])
            except ValueError:
                print "Insertion percent must be numeric."
            
        if 'deletionPercent' in self.optionsDictionary and self.optionsDictionary['deletionPercent'] is not None:
            try:
                self.deletionPercent = float(self.optionsDictionary['deletionPercent'])
            except ValueError:
                print "Deletion percent must be numeric."
            
        if 'convertFormat' in self.optionsDictionary and self.optionsDictionary['convertFormat'] is not None:
            self.convertFormat = self.optionsDictionary['convertFormat']        
            
        if len(self.repeatBase) != len(self.repeatRegionCount) or len(self.repeatBase) != len(self.repeatBaseLength) or len(self.repeatRegionCount) != len(self.repeatBaseLength):
            print "Repeat bases, repeat region count, and repeat base length must have the same number of items specified."
            sys.exit()
                
        if len(self.segment) != len(self.segmentCount):
            print "Segment and segment count must have the same number of items specified."
            sys.exit()
            
        if self.heterozygosity == 0 and (self.insertionPercent > 0 or self.deletionPercent > 0):
            print "Heterozygosity must be specified for insertions and deletions."
            sys.exit()        
        
    def generateFasta(self):
        i = 0
        baseCharacterCountTotal = 0
        while i < len(self.repeatBase):
            baseCharacterCountTotal += (len(self.repeatBase[i]) * int(self.repeatBaseLength[i]) * int(self.repeatRegionCount[i]))
            i += 1
            
        i = 0
        segmentCharacterCountTotal = 0
        while i < len(self.segment):
            segmentCharacterCountTotal += (len(self.segment[i]) * int(self.segmentCount[i]))
            i += 1
                
        numberOfRepeatCharacters = baseCharacterCountTotal + segmentCharacterCountTotal
    
        output = '>' + self.baseFilename + '\n'
        bases = []
        
        i = 0
        while i < (self.numberOfBases - numberOfRepeatCharacters):
            if i <= self.gcContentPercentage * self.numberOfBases:
                bases.append(self.baseList[random.randint(0,1)])
            else:
                bases.append(self.baseList[random.randint(2,3)])
            i += 1
            
        random.shuffle(bases)
        
        i = len(self.repeatRegionCount)
        while i > 0:
            k = 0
            while k < int(self.repeatRegionCount[i - 1]):
                index = random.randint(0, len(bases))
                l = 0
                while l < int(self.repeatBaseLength[i - 1]):
                    bases.insert(index + l, self.repeatBase[i - 1].strip())
                    l += 1
                k += 1
            i -= 1
        
        i = len(self.segmentCount)
        while i > 0:
            x = 0
            while x < int(self.segmentCount[i - 1]):
                index = random.randint(0, len(bases))
                l = 0
                while l < len(self.segment[i - 1]):
                    bases.insert(index + l, self.segment[i - 1][l].strip())
                    l += 1
                x += 1
            i -= 1
        
        sortedBases = []
        i = 0
        while i < len(bases):
            sortedBases.append(bases[i])
            i += 1
                
        if self.heterozygosity > 0:
            sortedBases += sortedBases
            
            i = 0
            changedBases = []
            numberOfBasesToChange = int((len(sortedBases) - 1) * self.heterozygosity)
            while i < numberOfBasesToChange:
                index = random.randint(self.numberOfBases, len(sortedBases) - 1)
                while changedBases.count(index) > 0:
                    index = random.randint(self.numberOfBases, len(sortedBases) - 1)
                
                currentBase = sortedBases[index]
                sortedBases[index] = changeBase(currentBase)
                changedBases.append(index)
                i += 1
                
            i = 0
            numberOfDeletions = int((len(sortedBases) - 1) * deletionPercent)
            while i < numberOfDeletions:
                index = random.randint(self.numberOfBases, len(sortedBases) - 1)
                base = baseList[random.randint(0, len(baseList) - 1)]
                sortedBases.pop(index)
                sortedBases.append(base)
                i += 1        
            
            i = 0
            numberOfInsertions = int((len(sortedBases) - 1) * insertionPercent)
            while i < numberOfInsertions:
                index = random.randint(self.numberOfBases, len(sortedBases) - 1)
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
        
        with open(self.filename, 'w') as fa:
            fa.write(output)
            print "Wrote " + self.filename
        
        qualFilename = self.generateQual(totalBases)
        
        if self.convertFormat is not None:
            convert(convertFormat, qualFilename)        
    
    def changeBase(self, currentBase):
        newBase = currentBase
        while newBase == currentBase:
            index = random.randint(0, len(self.baseList) - 1)
            newBase = self.baseList[index]
            
        return newBase
    
    def convert(self, qualFilename):
        print "\nFastA file information:"
        for seq_record in SeqIO.parse(self.filename, "fasta"):
            print str(seq_record)  
        
        try:
            with open(baseFilename + ".fastq", "w") as q:
                records = PairedFastaQualIterator(open(self.filename), open(qualFilename))
                count = SeqIO.write(records, q, self.convertFormat)
                print "Converted %i records" % count
        except ValueError, e:
            print "Encountered error converting file: " + str(e)
        
    def generateQual(self, numberOfBases):
        qualFilename = self.baseFilename + ".qual"
        quality = ">" + self.baseFilename + "\n"
        
        i = 0
        while i < numberOfBases:
            quality += str(self.getQuality()) + " "
            if i > 0 and i % 79 == 0:
                quality += "\n"
            i += 1
            
        with open(qualFilename, "w") as q:
            q.write(quality)
            print "Wrote " + qualFilename
            
        return qualFilename
        
    def getQuality(self):
        return random.randint(self.minQuality, self.maxQuality)
    