#!/usr/bin/env python2.7
import random

class names(object):
    def __init__(self):
        self.genus = []
        self.species = []
        f = "names.txt"
        infile = open(f, 'r')
        lines = infile.readlines()
        for line in lines:
            l = line.strip().split()
            self.genus.append(l[1])
            self.species.append(l[2])
    def get_name(self):
        g = random.choice(self.genus)
        s = random.choice(self.species)
        return g+" "+s

if __name__ == "__main__":
    n = names()
    for i in range(10):
        print(n.get_name())
    
