from z3 import *
from z3extra import *
from vol import VComp
from fractions import Fraction
from abstract import Interval

class ProbComp(object):
    def __init__(self, name, phi, vdist, finmax, randomize, infmax, z3qe, adapt, rotate, numHists = 5, histBound = 3, verbose=False,rot_digits=None):

        self.name = name

        args = [vdist, finmax, randomize, infmax, z3qe, adapt, rotate, numHists, histBound, verbose, rot_digits]

        self._lvol = VComp(phi, *args)
        self._uvol = VComp(Not(phi), *args)
        
        self._lg = self._lvol.getUV()
        self._ug = self._uvol.getUV()

        self._l = Fraction(0)
        self._u = Fraction(0)

    def qtime(self):
        return self._lvol.qtime() + self._uvol.qtime()

    def exact(self):
        return not (self._l + self._u < 1)

    def refine(self):
        done = True
        try:
            self._l = next(self._lg)
            done = False
        except StopIteration:
            pass
        try:
            self._u = next(self._ug)
            done = False
        except StopIteration:
            pass
        if done:
            raise StopIteration
        assert self._l + self._u <= 1

    def numSamples(self):
        return self._lvol.numSamples + self._uvol.numSamples

    @property
    def value(self):
        return Interval((self._l, 1-self._u))

    def __str__(self):
        return str(self.name) + " : [" + str(float(self.under())) + ", " + str(float(self.over())) + "]"
