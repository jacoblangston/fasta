# This is based on the ReadMaker perl script by Nozomu Okuda.
# Author: Jacob Langston.
# Url: https://github.com/jacoblangston/fasta
import random
from Bio import SeqIO
from Bio import Seq
from Bio import SeqRecord

class Exact(object):
    baseList = []
    baseFilename = None
    coverage = 0
    readLength = 0
    error = 0
    numberOfErrors = 0
    
    def __init__(self, baseList, fastaFilename, baseFilename, coverage, readLength):
	self.baseList = baseList
	self.baseFilename = baseFilename
        self.coverage = coverage
        self.readLength = readLength
	self.error = 0
                
	# Only one sequence is expected.  An exception will be thrown if more than one is detected.
	# Switch to SeqIO.parse if more than one sequence ever needs to be returned.  This will
	# require calling generateReads for each sequence.
	reference = SeqIO.read(fastaFilename, "fasta")

	print "Generating reads for %s" % fastaFilename
	self.generateReads(reference)
	print "Finished generating reads for %s" % fastaFilename
        
    def generateReads(self, sequence):
	filename = str(self.baseFilename) + ".fa"
	sequenceLength = len(sequence)
	randomLimit = sequenceLength - self.readLength
	cutoff = self.coverage * sequenceLength/self.readLength
	complete = 0
	records = []
        i = 0
        while i < cutoff:
            start = int(random.randint(0, randomLimit)) + 1;
	    tempSequence = sequence[start:start + self.readLength]

	    if random.random() < 0.5:
		tempSequence = tempSequence.reverse_complement()
    
	    records.append(self.generateErrors(tempSequence))

	    if (i * 100)/cutoff > complete:
		complete += 10
		print "%d%% of reads generated" % complete

            i += 1
	    
	SeqIO.write(records, filename, "fasta")
	print "Wrote " + filename

    def generateErrors(self, read):
	return read
    
class FourFiveFour(Exact):	
    def generateErrors(self, read):
	self.error = .0021
	output = ""
	self.numberOfErrors = 0
	numberOfCharactersToAdd = 0
	numberOfCharactersToDelete = 0	
	x = self.readLength - (self.error * self.readLength)
	
	i = 0
	while i < self.readLength:
	    char = read[i:i + 1]

	    if i > x:
		self.numberOfErrors += 1
		errorType = 0
		y = random.randint(0, 100)
		if y < 81:
		    errorType = 0
		else:
		    z = random.randint(0, 1)
		    if z == 0:
			errorType = 1
		    else:
			errorType = 2	    
		
		if errorType == 0: #substitution
		    output += self.baseList[random.randint(0, 3)]
		elif errorType == 1: #insertion
		    output += self.baseList[random.randint(0, 3)]
		    output += char
		    numberOfCharactersToDelete += 1
		elif errorType == 2: #deletion
		    numberOfCharactersToDelete += 1
		    numberOfCharactersToAdd += 1
	    else:
		output += char
		    
	    i += 1
	
	if numberOfCharactersToDelete > 0:
	    output = output[0: len(output) - numberOfCharactersToDelete]
	    
	if numberOfCharactersToAdd > 0:
	    i = 0
	    while i < numberOfCharactersToAdd:
		output += self.baseList[random.randint(0,3)]
		i += 1
	
	return output    
	
class IlluminaSingle(Exact):	
    def generateErrors(self, read):
	self.error = .0226
	output = ""
	self.numberOfErrors = 0
	numberOfCharactersToAdd = 0
	numberOfCharactersToDelete = 0	
	x = self.readLength - (self.error * self.readLength)
	
	i = 0
	while i < self.readLength:
	    char = read[i:i + 1]

	    if i > x:
		self.numberOfErrors += 1
		errorType = 0
	       
		y = random.randint(0, 100)
		if y < 98:
		    errorType = 0
		else:
		    z = random.randint(0, 1)
		    if z == 0:
			errorType = 1
		    else:
			errorType = 2
		
		if errorType == 0: #substitution
		    output += self.baseList[random.randint(0, 3)]
		elif errorType == 1: #insertion
		    output += self.baseList[random.randint(0, 3)]
		    output += char
		    numberOfCharactersToDelete += 1
		elif errorType == 2: #deletion
		    numberOfCharactersToDelete += 1
		    numberOfCharactersToAdd += 1
	    else:
		output += char
		    
	    i += 1
	
	if numberOfCharactersToDelete > 0:
	    output = output[0: len(output) - numberOfCharactersToDelete]
	    
	if numberOfCharactersToAdd > 0:
	    i = 0
	    while i < numberOfCharactersToAdd:
		output += self.baseList[random.randint(0,3)]
		i += 1
	
	return output
    
class IlluminaPaired(Exact):
    def generateReads(self, sequence):
	self.error = .0335
	filename = str(self.baseFilename) + ".fa"
	sequenceLength = len(sequence)
	records = []
	cutoff = self.coverage * sequenceLength/self.readLength
	randomLimit = sequenceLength - self.readLength
	complete = 0
	
	i = 0
	while i < cutoff:
	    start = int(random.randint(0, randomLimit)) + 1;
	    tempSequence = sequence[start:start + 500]
	    
	    firstSection = tempSequence[i:self.readLength]
	    records.append(self.generateErrors(firstSection, .0226))

	    lastSection = tempSequence[sequenceLength - self.readLength: self.readLength]
	    lastSection = lastSection.reverse_complement()
	    records.append(self.generateErrors(lastSection, self.error))
	    
	    if (i * 100)/cutoff > complete:
		complete += 10
		print "%d%% of reads generated" % complete	    
	    
	    i += 1
		
	SeqIO.write(records, filename, "fasta")
	print "Wrote " + filename    
	
    def generateErrors(self, read, error):
	output = ""
	self.numberOfErrors = 0
	numberOfCharactersToAdd = 0
	numberOfCharactersToDelete = 0	
	x = self.readLength - (error * self.readLength)
	
	i = 0
	while i < self.readLength:
	    char = read[i:i + 1]

	    if i > x:
		self.numberOfErrors += 1
		errorType = 0
	       
		y = random.randint(0, 100)
		if y < 98:
		    errorType = 0
		else:
		    z = random.randint(0, 1)
		    if z == 0:
			errorType = 1
		    else:
			errorType = 2
		
		if errorType == 0: #substitution
		    output += self.baseList[random.randint(0, 3)]
		elif errorType == 1: #insertion
		    output += self.baseList[random.randint(0, 3)]
		    output += char
		    numberOfCharactersToDelete += 1
		elif errorType == 2: #deletion
		    numberOfCharactersToDelete += 1
		    numberOfCharactersToAdd += 1
	    else:
		output += char
		    
	    i += 1
	
	if numberOfCharactersToDelete > 0:
	    output = output[0: len(output) - numberOfCharactersToDelete]
	    
	if numberOfCharactersToAdd > 0:
	    i = 0
	    while i < numberOfCharactersToAdd:
		output += self.baseList[random.randint(0,3)]
		i += 1
	
	return output    
