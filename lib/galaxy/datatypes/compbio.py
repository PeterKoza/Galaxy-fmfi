import logging, os, sys, time, sets, tempfile, shutil
import data
from galaxy import util
from galaxy.datatypes.sniff import *
from cgi import escape
import urllib
from bx.intervals.io import *
from galaxy.datatypes.metadata import MetadataElement

log = logging.getLogger(__name__)

class Descriptor(data.Text):
    """
    delimited data in descriptor format
    http://compbio.fmph.uniba.sk/rnarobo/rnarobo-readme.pdf
    """
    
    file_ext = "des"
    
    def sniff( self, filename ):
        try:
            fh = open(filename)
            self.header = []
            while True:
                line = fh.readline()
                if not line:
                    break # EOF
                line = line.strip()
                if line == "":
                    pass # empty line
                elif line[0] == "#":
                    pass # comments
                elif len(self.header) == 0:
                    self.header = line.split()
                    if self.checkDuplicity() == False:
                        return False
                else:
                    if self.checkLine(line) == False:
                        return False
            fh.close()
            if len(self.header) != 0:
                return False
            return True
        except Exception:
            pass
        return True

    def set_peek( self, dataset, is_multi_byte=False ):
        pass

    def set_meta(self, dataset, **kwd):
        pass

    def checkLine(self, line):
        """Checks if the first letter of line is {s,h,r}
           And if header consists element of the line.
           Uses functions: checkNumberOfNuclotitsHR(line), checkNumberOfNuclotitsS(line) 
            Args:
                line(str): one line of the descriptor
            Returns:
                bool: The return value. True for success, False otherwise.
        """
        if len(line) < 2:
            return False
        line = line.split()
        if line[0][0] == "h" or line[0][0] == "r":
            if line[0] in self.header and line[0]+"\'" in self.header:
                self.header.remove(line[0])
                self.header.remove(line[0]+"\'")
                if self.checkNumberOfNuclotitsHR(line) == False:
                    return False
            else:
                return False
                
        elif line[0][0] == "s":
            if line[0] in self.header:
                self.header.remove(line[0])
                if self.checkNumberOfNuclotitsS(line) == False:
                    return False
            else:
                return False
        else:
            return False
        return True



    def checkDuplicity(self):
        """Checks if every item of header is uniqe.
            Returns:
                bool: The return value. True for success, False otherwise.
        """
        for i in range(len(self.header)):
            for j in range(i+1, len(self.header)):
                if self.header[i] == self.header[j]:
                    return False
        return True



    def checkNumberOfNuclotitsHR(self, line):
        """Checks right format of element, exactly: "h" and "r" 
           Right formats:  h1 0:0:2 NNN**CC:GG**NNN:A
                           h1 0:0   NNN**CC:GG**NNN
                           r1 0:0:2 NNN**CC:GG**NNN:A   TGCA
                           r1 0:0   NNN**CC:GG**NNN     TGCA
           uses functions: tryInt(char), checkUIPAC(letters)
            Args:
                line(str): one line of the descriptor
            Returns:
                bool: The return value. True for success, False otherwise.
        """
        numbers = line[1].split(":")
        letters = line[2].split(":")

        if len(letters[0]) != len(letters[1]):
            return False

        if len(letters) > 3 or len(numbers) > 3:
            return False

        if len(letters) != len(numbers):
            return False

        if len(letters) == 3:
            if len(letters[2]) != 1:
                return False

        firstNum, secondNum = self.tryInt(numbers[0]), self.tryInt(numbers[1])
        thirdNum = 0
        if len(numbers) == 3:
            thirdNum = self.tryInt(numbers[2])
        if firstNum == None  or  secondNum == None  or  thirdNum == None: 
            return False
        if firstNum > len(letters[0])  or  secondNum > len(letters[1]):
            return False
        if line[0][0] == "h":
            if len(line) != 3:
                return False
        elif line[0][0] == "r":
            if len(line) != 4:
                return False
            if len(line[3]) != 4  or  "[" in line[3]  or  "]" in line[3]:
                return False
            letters.append(line[3])
        if self.checkUIPAC(letters) == False:
            return False
        return True



    def checkNumberOfNuclotitsS(self, line):
        """Checks right format of element, exactly: "h" and "r" 
           Right formats:  s1 0:0:2 NNN**CC:A
                           s1 0:0   NNN**CC
           uses functions: tryInt(char), checkUIPAC(letters)
            Args:
                line(str): one line of the descriptor
            Returns:
                bool: The return value. True for success, False otherwise.
        """
        numbers = line[1].split(":")
        letters = line[2].split(":")
        
        if len(letters) > 2 or len(numbers) > 2:
            return False

        if len(letters) != len(numbers):
            return False
        
        if len(letters) == 2:
            if len(letters[1]) != 1:
                return False
        firstNum = self.tryInt(numbers[0])
        secondNum = 0
        if len(numbers) == 2:
            secondNum = self.tryInt(numbers[1])
        if firstNum == None  or  secondNum == None:
            return False
        if firstNum > len(letters[0]):
            return False
        if self.checkUIPAC(letters) == False:
            return False
        return True

        

    def checkUIPAC(self, letters):
        """Checks right format of sequences 
           Right formats: NNN**CC, NNNNNN, [10]ATG
           uses function: getNextBracket(letter, j)
            Args:
                letters(array[][]): array of sequences
            Returns:
                bool: The return value. True for success, False otherwise.
        """
        for i in range(len(letters)):
            j = 0
            while j < len(letters[i]):
                if letters[i][j] == "[":
                    index = self.getNextBracket(letters[i], j)
                    if index == 0:
                        return False
                    j = index+1
                elif letters[i][j] not in UIPAC:
                    return False
                j += 1
        return True



    def getNextBracket(self, letters, j):
        """Finds end of number in sequence
           EG: [10]ATG -> 3
               NNN[3]N -> 5
            Args:
                letters(string): one sequence
                j(int): index of "[" char in sequence
            Returns:
                bool: The return value. True for success, False otherwise.
        """
        firstBracket = j
        secondBracket = 0
        while j < len(letters):
            if letters[j] == "]":
                if self.tryInt(letters[firstBracket+1: j]) == None or letters[firstBracket+1] == "0":
                    return secondBracket
                secondBracket = j
            j += 1
        return secondBracket



    def tryInt(self, num):
        """Tryes to convert char to int
            Args:
                num(string): number in string format
            Returns:
                int: number if it is possible on None value
        """
        try:
            return int(num)
        except:
            return None