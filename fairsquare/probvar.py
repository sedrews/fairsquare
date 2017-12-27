from abc import ABC, abstractmethod
from scipy.stats import norm
from z3 import *
from z3extra import *

class ProbVar(ABC):
    @abstractmethod
    def evalCDF(self, l, u):
        ...
    @abstractmethod
    def makeApproxDist(self, numHists, histBound):
        ...
    def __init__(self):
        self.adist = None


class GaussVar(ProbVar):
    def __init__(self, mean, variance):
        super().__init__()
        self.mean = mean
        self.variance = variance

    def makeApproxDist(self, numHists, histBound):
        # NOTE: numHists is a parameter, number of histograms
        # NOTE: histBound defines the farthest histogram, at var*histBound, e.g.,
        # if mean = 0 variance = 2 and histBound =3, then the farthest histogram on the right,
        # its upperbound will be at 0 + 6.
        
        #if the value is too large, when RealVal calls str() on the float argument,
        #it gets a representation with scientific notation that causes it to crash
        try:
            u = RealVal(self.mean + ((self.variance**0.5) * histBound))
            l = RealVal(self.mean - ((self.variance**0.5) * histBound))
        except Z3Exception:
            u = RealVal('{:.100f}'.format(self.mean + ((self.variance**0.5) * histBound)))
            l = RealVal('{:.100f}'.format(self.mean - ((self.variance**0.5) * histBound)))
        
        width = simplify((u - l) / numHists)

        bars = []

        while numHists > 0:
            # last bar, the one in the middle of the Gaussian
            if numHists == 1:

                assert((u.numerator_as_long() * l.denominator_as_long() \
                       - l.numerator_as_long() * u.denominator_as_long()) \
                       * width.denominator_as_long() \
                       == width.numerator_as_long()*u.denominator_as_long()*l.denominator_as_long())

                height = norm.pdf(z3ToFloat(u), loc=self.mean, scale=(self.variance**0.5))
                bars.append((l, u, height))
                break

            # get bar on the right
            uleft = simplify(u - width)
            rheight = norm.pdf(z3ToFloat(u),loc=self.mean, scale=(self.variance**0.5))
            #print rheight
            bars.append((uleft, u, rheight))

            # get bar on the left
            lright = simplify(l + width)
            lheight = norm.pdf(z3ToFloat(l), loc=self.mean, scale=(self.variance**0.5))
            bars.append((l, lright, lheight))

            # move right bound and left bound to the center
            u = simplify(u - width)
            l = simplify(l + width)

            numHists = numHists - 2

        self.adist = bars

    def evalCDF(self, l, u):
        if u is None:
            vu = 1
        else:
            vu = norm.cdf(z3ToFloat(u), loc=self.mean, scale=(self.variance**0.5))
            vu = self.under10(vu) #underapproximate the upper bound

        if l is None:
            vl = 0
        else:
            vl = norm.cdf(z3ToFloat(l), loc=self.mean, scale=(self.variance**0.5))
            vl = self.over10(vl) #overapproximate the lower bound

        #we use these approx's to soundly underapproximate the
        #weighted volume of the hyperrectangles with respect to
        #scipy's gaussian CDF precision
        #TODO ideally, we would underapproximate the floating point truncation
        #of the interval itself that we pass to norm.cdf, but
        #it has to take a floating point number as input, so it is less
        #clear how to be sound without a better understanding of python float

        assert vu > vl, "vu: " + str(vu) + " <= vl: " + str(vl)

        return vu - vl #these are python Fraction s

    #returns the rational number equal to
    #the floating point representation truncated to 1E-10 relative precision
    def under10(self, v):
        if v == 0.:
            return 0
        #we get the number of zeros that offset the significant digits
        shift = math.ceil(-math.log10(v)) - 1
        #we care about 10 significant digits
        shift = 10 + shift
        exp = 10 ** shift
        return Fraction(math.floor(v * exp), exp)

    def over10(self, v):
        if v == 0.:
            return 0
        shift = math.ceil(-math.log10(v)) - 1
        shift = 10 + shift
        exp = 10 ** shift
        return Fraction(math.floor(v * exp) + 1, exp)




class StepVar(ProbVar):
    def __init__(self, bars):
        super().__init__()
        self.bars = bars

    def makeApproxDist(self, numHists, histBound):
        #ignores the numHists and histBound parameters
        self.adist = self.bars

    def evalCDF(self, l, u):
        if u is None:
            u = self.stepWiseUpperBound()
        else:
            u = z3ToFrac(u)

        if l is None:
            l = self.stepWiseLowerBound()
        else:
            l = z3ToFrac(l)

        vol = 0

        for (c1, c2, p) in self.bars:
            overlap = self.evalOverlap((l, u), (Fraction(c1), Fraction(c2)))
            overlap = overlap / (c2 - c1)
            vol = vol + (overlap * p)

        return vol


    """ size of overlap of two pairs of lower/upper intervals """
    def evalOverlap(self, int1, int2):
        (l, u) = int1
        (c1, c2) = int2
        left = max(l, c1)
        right = min(c2, u)
        if right - left <= 0: #XXX is this check due to floating point error?
            return 0
        return right - left


    """ get upperbound on stepwise-continuous distribution """
    def stepWiseUpperBound(self):
        # take max of all upper bounds of all bars
        res = max([b[1] for b in self.bars])

        return res

    """ get lowerbound on stepwise-continuous distribution """
    def stepWiseLowerBound(self):
        # take minimum of all lower bounds of all bars
        res = min([b[0] for b in self.bars])

        return res

