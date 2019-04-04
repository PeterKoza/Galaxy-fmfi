import optparse

class RoboConverter():
    def __init__(self):
        pass

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
        endPos = int(parameters[1])
        strand = "+"
        chromosone = " ".join(parameters[2:])
        if startPos > endPos:
            strand = "-"
            chromosone = " ".join(parameters[2:-1])
            startPos, endPos = endPos, startPos
        chrom = chromosone.split(" ")[0]
        name = chrom + ":" + str(startPos) + "-" + str(endPos)
        res = chrom + "\t" + str(startPos) + "\t" + str(endPos) + "\t"
        res += name + "\t0\t" + strand + "\n"
        return res



def main():
    parser = optparse.OptionParser()
    converter = RoboConverter()

    (options, args) = parser.parse_args()
    input, output = args

    result = converter.roboToBed(input)

    out = open(output, "w")
    out.write(result)
    out.close()      


if __name__ == "__main__":
    main()