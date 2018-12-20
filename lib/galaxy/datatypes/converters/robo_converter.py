import optparse

class RoboConverter():
    def __init__(self):
        pass

    # ------ parse robo format to fasta
    def roboToFasta(self, roboFile):
        rf = open(roboFile);
        text = rf.read().split("\n")
        result = ""
        separators = 0
        for line in text:
            if separators != 3:
                if line == "":
                    result += "\n"
                else:
                    result += ";" + line + "\n"
            else:
                if line == "":
                    result += "\n"
                elif line[0] == "|"  or  line[0] in "ACGTU":
                    result += line.replace("|", "") + "\n"
                elif line[0] != "-":
                    positionsAndName = line.split()
                    startPos = str(int(positionsAndName[0]) + 1)
                    endPos = str(int(positionsAndName[1]) + 1)
                    seqName =  "_".join(positionsAndName[2:])
                    startWith = ""
                    if line[0] != ">":
                        startWith = ">"
                    result += startWith + seqName + ":" + startPos + "-" + endPos + "\n"

            if line != "":
                if line[0] == "-":
                    separators += 1
        rf.close()      
        return result


    # ------ parse robo format to BED
    def roboToBed(self, roboFile):
        rf = open(roboFile);
        text = rf.read().split("\n")
        result = ""
        separators = 0
        for line in text:
            if separators == 3:
                if line == "":
                    pass
                elif line[0] == "|"  or  line[0] in "ACGTU":
                    pass
                elif line[0] != "-":
                    result += self.createBadLine(line)

            if line != "":
                if line[0] == "-":
                    separators += 1
        rf.close()      
        return result


    def createBadLine(self, line):
        parameters = line.split()
        startPos = int(parameters[0])
        endPos = int(parameters[0])
        strand = "+"
        chromosone = " ".join(parameters[2:])
        if startPos > endPos:
            strand = "-"
            chromosone = " ".join(parameters[2:-1])
            startPos, endPos = endPos, startPos
        name = chromosone + ":" + str(startPos) + "-" + str(endPos)
        res = chromosone + "\t" + str(startPos) + "\t" + str(endPos) + "\t"
        res += name + "\t0\t" + strand + "\n"
        return res



def main():
    parser = optparse.OptionParser()
    converter = RoboConverter()

    (options, args) = parser.parse_args()
    input_file, output_file, option = args


    if option == "f":
        result = converter.roboToFasta(input_file)
    else:
        result = converter.roboToBed(input_file)

    out = open(output_file, "w")
    out.write(result)
    out.close()      


if __name__ == "__main__":
    main()