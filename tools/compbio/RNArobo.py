#!/usr/bin/python -w

# usage : python RNArobo.py <FASTA file> <output file>

import sys
import os
import time

#time.sleep(50)

fasta = sys.argv[1]
descriptor = sys.argv[2]
options = sys.argv[3].replace(',', '')
nratio = sys.argv[4]
output = sys.argv[5]


class ComputeRNArobo():
	def __init__(self, fasta, descriptor, options, nratio, output):
		self.fasta = fasta
		self.descriptor = descriptor
		self.options = options
		self.nratio = nratio
		self.output = output
		self.globalPath, self.localPath = os.getcwd().split("database")


	def run(self):
		params = self.getParameters()

		os.chdir(self.globalPath + "tools/compbio")

		os.chdir('RNArobo')
		os.system("./rnarobo " + params + " > " + self.output)
		os.chdir("..")

		if "f" in self.options:		# if format is .fasta than cut ---- Search Done ---- part of file 
			occurrencesFile = open(self.output, "r")
			occurrences = occurrencesFile.read()
			occurrencesFile.close()
			occurrences = self.cutDownSearchPart(occurrences)
			occurrencesFile = open(self.output, "w")
			occurrencesFile.write(occurrences)
			occurrencesFile.close()



	def getParameters(self):
		result = ""
		if self.options != "None":
			result += "-" + self.options + " "
		result += "--nratio " + self.nratio + " "
		result += self.descriptor + " " + self.fasta
		return result


	def cutDownSearchPart(self, occurrences):
		content = occurrences.split("----- SEARCH DONE -----")
		return content[0]



ComputeRNArobo(fasta, descriptor, options, nratio, output).run()

