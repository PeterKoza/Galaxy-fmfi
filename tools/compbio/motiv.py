class Motif():
    def __init__(self, line):
        # TODO: check motifmap
        self.name = ""
        self.weight = 0
        self.score = 0
        self.per_nucleotide = 0
        self.tail = 0
        self.unw_score = 0
        self.elements = []
        self.parseLine(line)


    def parseLine(self, line):
        line = line.split()
        self.setName(line.pop(0))
        self.setElements(line.pop())
        for param in line:
            if param[-1] == "w":
                self.setWeight(param)

            if param[-1] == "t":
                self.setTail(param)
            elif param == "pn":
                self.per_nucleotide = 1
            else:
                #TODO raise error
                pass


    def setName(self, name):
        self.name = name[:-1]

    def setElements(self, elements):
        self.elements = elements[1:-1].split(',')

    def setWeight(self, w):
        if type(w) is int  or  type(w) is float:
            self.weight = w
        else:
            self.weight = int(w[:-1])

    def getWeight(self):
        return self.weight

    def setTail(self, t):
        if type(t) is int  or type(t) is float:
            self.tail = t
        else:
            self.tail = int(t[:-1])

    def setUnweightScore(self, us):
        self.unw_score = us