import optparse
import os
from motif import Motif

class FoldFilter():
    
    def __init__(self):
        self.motifMap = []
        self.filter = []
        self.data = []
        self.finalResult = []
        self.globalPath, self.localPath = os.getcwd().split("database")

    def parseMotifMap(self, desFile):
        file = open(desFile, "r")
        line = file.readline()
        while line:
            line = line.strip()
            if line == "":
                pass
            elif line[0] == "#":
                pass
            else:
                self.motifMap = line.split()
                break
            line = file.readline()
        file.close()

    def parseFilterConfig(self, ftrFile):
        file = open(ftrFile, "r")
        line = file.readline()
        while line:
            line = line.strip()
            if line == "":
                pass
            elif line[0] == "#":
                pass
            else:
                self.filter.append(Motif(line))
            line = file.readline()
        file.close()

    def recalculateWeights(self):
        wsum = 0
        for submotif in self.filter:
            wsum += submotif.getWeight()
        for submotif in self.filter:
            w = submotif.getWeight() * (float(10) / float(wsum))
            submotif.setWeight(w)

    def parseFastaFile(self, faFile):
        f = open(faFile)
        lines = f.read().split("\n")
        f.close()
        for i in range(len(lines)):
            if lines[i] == "":
                pass
            elif lines[i][0] == "|":
                sequence = lines[i][1:-1].replace("|", " ")
                self.data.append((lines[i-1], sequence))


    def computeScores(self):
        for data in self.data:
            sequence = data[1].split()
            seqScore = 0
            subMotifsHash = dict()
            for subMotif in self.filter:
                subSeq = []
                map = []
                for i in range(len(self.motifMap)):
                    if self.motifMap[i] in subMotif.elements  or  self.motifMap[i][:-1] in subMotif.elements:
                        if len(sequence) - 1 < i:
                            print "sequence: ", sequence
                            print "motifMap: ", self.motifMap
                            print 'Descriptor is not consistent with RNArobo output.'
                            raise Exception('Descriptor is not consistent with RNArobo output.')
                        subSeq.append(sequence[i])
                        map.append(self.motifMap[i])
                nested = self.isNested(map)
                pseudoknotted = self.isPseudoknotted(map)
                if nested:
                    tmpSeq = "".join(subSeq[:len(subSeq)/2]) + "&" + "".join(subSeq[len(subSeq)/2:])
                    structure, min_en = self.runRNAcofold(tmpSeq)
                elif pseudoknotted:
                    tmpSeq = "".join(subSeq)
                    structure, min_en = self.runDotKnot(tmpSeq)
                else:
                    tmpSeq = "".join(subSeq)
                    structure, min_en = self.runRNAfold(tmpSeq)

                subMotif.unw_score = -1 * min_en
                if subMotif.per_nucleotide:
                    subMotif.unw_score = float(subMotif.unw_score) / (0.1 * len(structure))
                if subMotif.tail: # tail
                    pairPos = [None] * len(structure)
                    stack1 = []
                    stack2 = []
                    for ix in range(len(structure)):
                        if structure[ix] == "(":
                            stack1.append(ix)
                        elif structure[ix] == "[":
                            stack2.append(ix)
                        elif structure[ix] == ")":
                            top = stack1.pop()
                            pairPos[top] = ix
                            pairPos[ix] = top
                        elif structure[ix] == "]":
                            top = stack2.pop()
                            pairPos[top] = ix
                            pairPos[ix] = top
                    count = 0
                    if subMotif.tail > structure:
                        bound = len(structure) / 2
                    else:
                        bound = subMotif.tail

                    for ix in range(bound):
                        if pairPos[ix] and pairPos[ix] >= len(structure) - bound:
                            count += 1
                    subMotif.unw_score = (float(count) / float(bound)) * 100
                subMotif.score = subMotif.unw_score * subMotif.weight
                subMotifHash = dict()
                subMotifHash["unw_score"] = subMotif.unw_score
                subMotifHash["score"] = subMotif.score
                subMotifHash["min_en"] = min_en
                subMotifsHash[subMotif.name] = subMotifHash
                seqScore += subMotif.score
            self.finalResult.append((seqScore, subMotifsHash))


    def getTabularFormat(self):
        names = [filter.name for filter in self.filter]
        result = "#sequence\ttotal score\t" + "\t".join(names) + "\n"
        for i, res in enumerate(self.finalResult):
            seqName = self.data[i][0].split()
            seqName = seqName[2] + ":" + seqName[0] + "-" + seqName[1]
            result += seqName + "\t"
            result +=  str(res[0]) + "\t"
            scores = []
            for name in names:
                score = res[1][name]["score"]
                if score == 0:
                    score = int(score)
                scores.append(str(score))
            result += "\t".join(scores) + "\n"
        return result


    def isNested(self, map):
        nested = False
        if len(map)%2 == 0:
            nested = True
            for i in range(len(map)/2):
                if map[i]+"'" != map[len(map)-i-1]:
                    nested = False
        return nested

    def isPseudoknotted(self, map):
        pseudoknotted = False
        stack = []
        for el in map:
            if el[0] in ["r", "h"]  and  el[-1] != "'":
                stack.append(el+"'")
            elif el[0] in ["r", "h"]:
                if stack.pop() != el:
                    pseudoknotted = True
        return pseudoknotted

    def runRNAcofold(self, seq):
        tmpSeq = seq.split("&")
        if tmpSeq[0] != tmpSeq[1] and tmpSeq[0] != tmpSeq[1][:-1]:
            result = os.popen("echo '" + seq + "'|RNAcofold --noPS")
            result = result.read().split()
            structure = result[1].replace("&", "")
            min_en = result[-1][:-1].replace("(", "")
            min_en = float(min_en)
        else:
            # RNAcofold creates warnings, which is considered as error in galaxy system
            # input seq cases: GGTT&GGTT   or   GGTT&GGTTT
            # therefore we use python print function
            print "WARNING: Both input strands are identical, thus inducing rotationally symmetry! Symmetry correction might be required to compute actual MF"
            structure, min_en = '.'*len(seq), 0
        return structure, min_en

    def runRNAfold(self, seq):
        result = os.popen("echo '" + seq + "'|RNAfold --noPS")
        result = result.read().split()
        structure = result[1]
        min_en = result[-1][:-1].replace("(", "")
        min_en = float(min_en)
        return structure, min_en

    def runDotKnot(self, seq):
        os.chdir(self.globalPath + "tools/compbio")
        
        os.chdir("dotknot")
        file = open("file.txt", "w")
        file.write(">tmp\n" + seq + "\n")
        file.close()
        resFile = os.popen("python dotknot.py file.txt -klg")
        os.chdir('..')
        result = resFile.read()
        resFile.close()

        os.chdir(self.globalPath + "database" + self.localPath)
        
        resultSize = result.split("\n")
        if len(resultSize) > 25:
            lines = result.split("energy:")[1].strip()
            lines = lines.split("\n")
            start, end, min_en = lines[0].split()
            start, end, min_en = int(start), int(end), float(min_en)
            tmpStruct = lines[2]
            structure = (start - 1) * "."
            structure += tmpStruct
            structure += (len(tmpStruct) - end) * "."
        else:
            structure, min_en = self.runRNAfold(seq)
        return structure, min_en



def main():
    parser = optparse.OptionParser()
    foldFilter = FoldFilter()

    (options, args) = parser.parse_args()
    roboFile, desFile, ftrFile, output_file = args

    foldFilter.parseMotifMap(desFile)
    foldFilter.parseFilterConfig(ftrFile)
    foldFilter.recalculateWeights()
    foldFilter.parseFastaFile(roboFile)
    foldFilter.computeScores()
    result = foldFilter.getTabularFormat()
    
    out = open(output_file, "w")
    out.write(str(result))
    out.close() 


if __name__ == "__main__":
    main()
    print "10" 
