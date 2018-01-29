#!/usr/bin/python -w

# usage : python RNArobo.py <FASTA file> <output file>

import sys
import os



fasta = sys.argv[1]
descriptor = sys.argv[2]
options = sys.argv[3].replace(',', '')
nratio = sys.argv[4]
output = open(sys.argv[5], 'w')


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
		os.system("./rnarobo " + params + " > ../occurrences.txt")
		os.chdir("..")

		occurrences = open("occurrences.txt")
		output.write(occurrences.read())
		occurrences.close()

		os.chdir(self.globalPath + "database" + self.localPath)

		output.close()

	def getParameters(self):
		result = ""
		if self.options != "None":
			result += "-" + self.options + " "
		result += "--nratio " + self.nratio + " "
		result += descriptor + " " + fasta
		return result



ComputeRNArobo(fasta, descriptor, options, nratio, output).run()

