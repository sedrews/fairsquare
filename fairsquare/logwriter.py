import sys
from functools import reduce

def cs(ell):
    return reduce(lambda x,y : x + "," + y, ell)

class LogWriter:
    def __init__(self, filename, csvcols):

        self.csv = open(filename + ".csv", "w")
        self.log = open(filename + ".log", "w")

        self.csv.write(cs(csvcols) + "\n")

        self.log.write(reduce(lambda x,y : x + " " + y, sys.argv) + "\n")
        self.log.write("writing data to " + filename + ".csv\n")
    
    def data(self, cols):
        self.csv.write(cs(cols) + "\n")

    def message(self, message):
        self.log.write(message + "\n")

    def close(self):
        self.csv.close()
        self.log.close()
