import sys
import os
import optparse


class JoinBedWithFilter():
    def __init__(self):
        parser = optparse.OptionParser()
        (options, args) = parser.parse_args()
        bed_file, FoldFilter_outup_file, output_file = args
        bed = self.getBEDdictionary(bed_file)
        scores = self.getFilterArray(FoldFilter_outup_file)
        result = ""
        for score in scores:
            bed[score[0]][4] = score[1]
            result += "\t".join(bed[score[0]]) + "\n"
        self.writeOutput(output_file, result)


    def getBEDdictionary(self, bed_file):
        bed_file = open(bed_file, "r")
        bed_text = bed_file.read().strip()
        bed_file.close()
        lines = bed_text.split("\n")
        lines = [line.split("\t") for line in lines]
        res = dict()
        for line in lines:
            res[line[3]] = line
        return res

    def getFilterArray(self, FoldFilter_outup_file):
        bed_file = open(FoldFilter_outup_file, "r")
        bed_text = bed_file.read().strip()
        bed_file.close()
        lines = bed_text.split("\n")
        lines = [line.split("\t") for line in lines[1:]]
        return lines

    def writeOutput(self, output_file, result):
        output = open(output_file, "w")
        output.write(str(result))
        output.close() 


j = JoinBedWithFilter()

