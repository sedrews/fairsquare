from z3 import *
from z3extra import *
from scipy.stats import norm
import numpy as np
import copy
import subprocess
from collections import deque
from redlogInterface import *
import rotationmath as rm
from numericalPhi import numericalPhi
import time
from functools import reduce
from probvar import GaussVar, StepVar

def wait():
    input("Press something...")


class VComp:

    """ Volume computation class """

    def __init__(self, phi, vdist, finmax, randomize, infmax, z3qe, adapt, rotate, numHists = 5, histBound = 3, verbose=False,rot_digits=None):
        # formula phi
        self.phi = close(phi)
        #print "phi", phi
        #print "CLOSED phi", self.phi

        # map from phi variables to distributions
        self.vdist = {}
        for v in vdist:
            if vdist[v][0] == 'S':
                self.vdist[v] = StepVar(*vdist[v][1:])
            elif vdist[v][0] == 'G':
                self.vdist[v] = GaussVar(*vdist[v][1:])
            else:
                assert(False)
        #self.vdist = vdist
        #print "VDIST,", vdist

        # collect all variables in phi
        self.phiVars = get_vars(self.phi)
        #print "PHIVARS,", self.phiVars

        # the parameters for approxGauss
        self.numHists = numHists
        self.histBound = histBound

        # approximate distributions for optimizations
        self.makeApproxDists()

        # each dimension is associated with variable that approximates
        # its density
        self.approxVars = self.getApproxVars()

        # Gaussian approximation constraints
        self.approxConsts = []

        # Sampling solver
        self.sampler = Solver()

        # whether or not to maximize finite bounds on samples
        self.finmax = finmax
        self.randomize = randomize
        if finmax:
            np.random.seed(randomize)

        # whether or not to drop infinite bounds on samples
        self.infmax = infmax

        # whether to use z3 quantifier elimination instead of redlog
        self.z3qe = z3qe

        # how much time was spent in quantifier elimination
        self.qelimtime = 0

        # variables for minimum samples size
        self.samplelb = 1.0
        self.decayrate = 0.5

        # variables for dynamically adjusting histBound
        self.adapt = adapt
        self.pastvol = -1

        # Verbose or not?
        self.verb = verbose
        
        # Do we rotate? Is the rotation something
        # other than the identity?
        self.rotate = rotate
        self.not_identity = False
        if self.rotate == 'forced':
            self.rel_gvars = self.getGaussians()
            self.rmat = np.identity(len(self.rel_gvars))

        # Set a digit limit on the rational
        # vectors.
        if rot_digits is None:
            self.rot_digits = min(max(8 - len(self.getGaussians()), 3), 5)
        else:
            self.rot_digits = rot_digits

        #keep track of the number of samples used for volume computation so far
        self.numSamples = 0


    def qtime(self):
        return self.qelimtime

    """ variables used to measure density in every dimension """

    def getApproxVars(self):
        d = {}

        for v in self.phiVars:
            d[v] = Real("d" + str(v))

        return d


    def makeApproxDists(self):
        for v in self.vdist:
            self.vdist[v].makeApproxDist(self.numHists, self.histBound)

    """ create approximation constraints """

    def getApproxConstraints(self, lbs, ubs):

        if self.approxConsts != []:
            return self.approxConsts

        consts = []
        for v in self.phiVars:

            # some vars won't have distributions
            if v not in self.vdist:
                continue

            (l, u) = (lbs[v], ubs[v])

            dist = self.vdist[v].adist

            cv = self.approxVars[v]

            weights = []

            i = 0
            for (lb, ub, pb) in dist:
                cvb = Real(str(cv) + str(i))
                i = i + 1

                # compute overlap between bounds of v and the bin
                lower = z3max(lb, l)
                upper = z3min(ub, u)

                cvb = If(lower >= upper, 0, (upper - lower) * pb)

                if self.verb:
                    print("CVB", cvb)
                weights.append(cvb)

            vconst = cv == Sum(*weights)
            consts.append(vconst)

        self.approxConsts = consts
        #print consts
        #exit(1)
        return consts

    def adaptApprox(self, svol):
        if self.pastvol == -1:
            self.pastvol = svol
        else:
            if svol < self.pastvol / 10:
                self.histBound = 2 * self.histBound
                self.makeApproxDists()
                self.pastvol = self.pastvol / 10


    """ maps each variable to a fresh upper bound variable """

    def getUbs(self):
        ubs = {}

        for v in self.phiVars:
            ubs[v] = Real('u_' + str(v))

        return ubs

    """ maps each variable to a fresh lower bound variable """

    def getLbs(self):
        lbs = {}

        for v in self.phiVars:
            lbs[v] = Real('l_' + str(v))

        return lbs

    """ use mathematica to integrate over form """

    def integrate(self, form):
        #XXX: this code probably doesn't work anymore
        mfile = open('mfile.mat', 'w', )
        alldists = []
        mfile.write("Put[ \"math.out\" ]\n")
        for v in self.phiVars:
            (t, mean, var) = self.vdist[v]
            distvstr = "dist" + str(v)
            alldists.append(distvstr)

            if t == 'G':
                distvDecl = distvstr + " = " + \
                    str(v) + " \\[Distributed] NormalDistribution[" + \
                    str(mean) + "," + str(var) + "]\n"
                if self.verb:
                    print(distvDecl)
                mfile.write(distvDecl)

        alldists = reduce(lambda v, v1: v + ", " + v1, alldists)
        alldists = "{" + alldists + "}"

        #options = ", AccuracyGoal -> Infinity"
        options = ""
        prob = "N[Probability[" + \
            z3ToMath(form) + ", " + alldists + options + "]]>>>\"math.out\"\n"
        mfile.write(prob)
        mfile.write("Quit[]")

        if self.verb:
            print("===========================")
            print("Calling Mathematica")

            mfile.close()

            subprocess.call(
                '/Applications/Mathematica.app/Contents/MacOS/MathKernel -noprompt -run \"<<mfile.mat\"', shell=True)

            print("DONE")
            print("===========================")

        resfile = open('math.out', 'r')
        res = resfile.readline().rstrip()
        return float(res)

    def getRectPhi(self, lbs, ubs):
        # conjunction of formulas: AND_v  l_v <= v <= u_v
        lhs = []
        for v in self.phiVars:
            lhs.append(And(v >= lbs[v], v <= ubs[v]))

        lhs = And(*lhs)

        # qelim tactic
        #t = Tactic('qe')

        # rectangles are FORALL V . lhs => phi, where V are variables in phi
        # performa quelim

        startqelim = time.time()

        if self.z3qe:
            if self.verb:
                print("\nCALLING QELIM Z3")
                print(self.phi)
            res = qelimForallDNF(self.phiVars, Implies(lhs, self.phi))
        else:
            if self.verb:
                print("\nCALLING QELIM REDLOG")
                print(self.phi)
            res = qelimRedlog(self.phiVars, Implies(lhs, self.phi))
        if self.verb:
            print(res)
        self.qelimtime = time.time() - startqelim

        # NOTE: sanity check to ensure result of qelimRedlog
        # is contained in quantified formula
        s = Solver()
        s.add(res)
        s.add(Not(Implies(lhs, self.phi)))
        assert(s.check() == unsat)

        # exit(1)

        # qelimForallDNF(self.phiVars, Implies(lhs, self.phi))
        rectsQelim = res

        # FIXME: Remove consistency constraints -- not needed anymore
        # consistency constraints: lower bounds less than upper bounds
        consistencyConst = [lbs[v] < ubs[v] for v in self.phiVars]
        consistencyConst = bigAnd(consistencyConst)

        return And(rectsQelim, consistencyConst)

    """ formula representation of sample """

    def sampleToFormula(self, sample):
        f = True
        for v in self.phiVars:
            f = And(f, v <= sample[v][1], v >= sample[v][0])

        return f

    """ checks sample included in phi and does not overlap with other rects """

    def checkSample(self, sample, rects, lbs, ubs):
        # check sample contained in phi -- sample => phi
        s = Solver()
        s.add(Not(phi))

        s.add(self.sampleToFormula(sample))

        if s.check() == sat:
            return False

        # check sample does not overlap with existing rectangles
        s = Solver()
        s.add(rects)
        for v in self.phiVars:
            s.add(lbs[v] == sample[v][0])
            s.add(ubs[v] == sample[v][1])

        if s.check() == sat:
            return True

        return False


    """ drop constraints from sample -- i.e., make bounds infinite """

    def maximizeSample(self, sample, rects, lbs, ubs):
        s = Solver()

        s.add(Not(rects))

        umap = {}
        soft = []

        # check that sample => rects, and remove constraints from unsat core
        for v in self.phiVars:
            # temporary Boolean var for unsat core
            umap[v] = (Bool('corel_' + str(v)), Bool('coreu_' + str(v)))

            # each Boolean var implies a lower/upper bound
            s.add(Implies(umap[v][1], ubs[v] <= sample[v][1]))
            s.add(Implies(umap[v][0], lbs[v] >= sample[v][0]))

            s.add(ubs[v] > lbs[v])

            soft.append(umap[v][0])
            soft.append(umap[v][1])

        if self.verb:
            print(s)
            print("Maximizing sample")

        s.check(soft)
        if s.check(soft) == sat:
            if self.verb:
                print(s.model())
        assert (s.check(soft) == unsat)

        # all constraints not appearing in unsat_core can be safely dropped
        if self.verb:
            print(s.unsat_core())
        toDrop = [v for v in soft if v not in s.unsat_core()]

        # drop all constraints in toDrop
        for c in toDrop:
            if self.verb:
                print("CORE: ", c)

            # finds variable denoted by constraint c and whether it's
            # upper or lower bound is to be removed
            for v_ in self.phiVars:
                if umap[v_][0] is c:
                    v = v_
                    upper = False
                    break
                if umap[v_][1] is c:
                    v = v_
                    upper = True
                    break

            # None indicates pos/neg infinity
            if upper:
                sample[v] = (sample[v][0], None)
            else:
                sample[v] = (None, sample[v][1])

        return sample

    """ try to optimize all finite bounds """
    def extendSample(self, sample, rects, lbs, ubs):
        sanity = Solver()
        sanity.add(Not(rects))
        sanity.check()
        M = Optimize()
        M.add(And(*self.sampler.assertions()))
        updated = False
        varordering = self.phiVars
        if self.randomize is not None:
            varordering = np.random.permutation(self.phiVars)
        for v in varordering:
            M.push()
            sanity.push()
            for v_ in [x for x in self.phiVars if x is not v]:
                sanity.add(lbs[v_] < ubs[v_])
                M.add(lbs[v_] == sample[v_][0])
                sanity.add(lbs[v_] == sample[v_][0])
                M.add(ubs[v_] == sample[v_][1])
                sanity.add(ubs[v_] == sample[v_][1])
            #lower bound
            M.push()
            sanity.push()
            M.add(ubs[v] == sample[v][1])
            if sample[v][1] is not None:
                sanity.add(ubs[v] == sample[v][1])
            M.minimize(lbs[v])
            assert(M.check() == sat)
            if z3ToFrac(M.model()[lbs[v]]) < z3ToFrac(sample[v][0]):
                sanity.add(lbs[v] == M.model()[lbs[v]])
                if sanity.check() == unsat:
                    if self.verb:
                        print(lbs[v], "updated from", z3ToFloat(sample[v][0]), "to", z3ToFloat(M.model()[lbs[v]]))
                    updated = True
                    sample[v] = (M.model()[lbs[v]], sample[v][1])
            M.pop()
            sanity.pop()
            #upper bound
            M.add(lbs[v] == sample[v][0])
            sanity.add(lbs[v] == sample[v][0])
            M.maximize(ubs[v])
            assert(M.check() == sat)
            if z3ToFrac(M.model()[ubs[v]]) > z3ToFrac(sample[v][1]):
                sanity.add(ubs[v] == M.model()[ubs[v]])
                if sanity.check() == unsat:
                    if self.verb:
                        print(ubs[v], "updated from", z3ToFloat(sample[v][1]), "to", z3ToFloat(M.model()[ubs[v]]))
                    updated = True
                    sample[v] = (sample[v][0], M.model()[ubs[v]])
            M.pop()
            sanity.pop()

        if self.verb:
            print("maximal sample:", sample)
        #if updated:
        #    raw_input()

        return sample

    """ Gets a sample rectangle by heuristically picking models and expanding """

    def getSampleHeuristic(self, rects, lbs, ubs):
        # we want to get a model of phi
        s = Solver()
        s.add(phi)

        s.add(rects)

        # to ensure we do not overlap with rects
        for v in self.phiVars:
            s.add(v == lbs[v])
            s.add(v == ubs[v])

        # no more samples
        if s.check() == unsat:
            return None

        m = s.model()

        # sample is a map from variables to upper and lower bounds
        sample = {}

        for v in self.phiVars:
            sample[v] = (m[lbs[v]], m[ubs[v]])

        if self.verb:
            print("hsample before: ", sample)
        sample = self.maximizeSampleHeuristic(sample, rects, lbs, ubs)

        return sample

    """ Gets a sample rectangle """

    def getSample(self, rects, lbs, ubs, obj=None, s=True):
        sampler = self.sampler

        if sampler.check() == unsat:
            return None

        c = 0
        sampler.push()
        c = c + 1
        # add optimization objective function
        # maximize SUM_v upper_v - lower_v
        # i.e., maximize perimiter of rectangle

        # assert approximation
        approxConsts = self.getApproxConstraints(lbs, ubs)
        sampler.add(approxConsts)

        if self.verb:
            print("Getting a sample")

        # NOTE: Change parameters
        lowerbound = self.samplelb
        decay = self.decayrate

        sampler.push()
        c = c + 1

        for v in self.phiVars:
            sampler.add(self.approxVars[v] >= lowerbound)

        # FIXME: make parameter
        cutoff = 1e-10
        while sampler.check() == unsat:
            sampler.pop()
            c = c - 1
            if lowerbound < cutoff:
                #raw_input("it happened")
                lowerbound = 0
            sampler.push()
            c = c + 1
            lowerbound = lowerbound * decay

            for v in self.phiVars:
                sampler.add(self.approxVars[v] >= lowerbound)

        self.samplelb = lowerbound

        sampler.pop()
        c = c - 1
        sampler.pop()
        c = c - 1

        # check that we did not forget to pop all the way
        assert (c==0)

        m = sampler.model()
        if self.verb:
            print("model from sampler", m)

        # sample is a map from variables to upper and lower bounds
        sample = {}

        for v in self.phiVars:
            sample[v] = (m[lbs[v]], m[ubs[v]])

        #we don't want to try to optimize parameters that have no bound
        if self.verb:
            print("before: ", sample)
        if self.finmax:
            sample = self.extendSample(sample, rects, lbs, ubs)
        if self.infmax:
            sample = self.maximizeSample(sample, rects, lbs, ubs)
        #sanity check our sample
        s = Solver()
        s.add(Not(rects))
        for v in self.phiVars:
            s.add(lbs[v] < ubs[v]) #this has to be a strict inequality
            if sample[v][0] is not None:
                s.add(lbs[v] == sample[v][0])
            if sample[v][1] is not None:
                s.add(ubs[v] == sample[v][1])
        assert(s.check() == unsat)


        return sample

    """ negates the rectangle defined by the sample rect """

    def negateRect(self, rect, lbs, ubs):
        overlapRect = []
        sameRect = []

        # to negate the rectangle, for one of the dimensions v,
        # we need to ensure the upper/lower bounds do not overlap
        # with the current upper and lower bounds.
        # we also need to ensure we do not get the same rectangle

        for v in self.phiVars:

            # skip variables that do not have distributions
            # these are dataflow variables
            if v not in self.vdist:
                continue

            # lower and upper bounds of dimension v in rect
            l = rect[v][0]
            u = rect[v][1]

            # FIXME: NEED TO FIX constraints > --> >=
            # constraints that specify all rectangles that do not
            # overlap with rect
            # None specifies upper bound is infty

            if u is None:
                aboveU = False
            else:
                aboveU = And(lbs[v] >= u, ubs[v] >= u)

            if l is None:
                underL = False
            else:
                underL = And(lbs[v] <= l, ubs[v] <= l)

            overlapRect.append(Or(aboveU, underL))

            # formula representation of rectangle

            if l is None:
                vl = True
            else:
                vl = lbs[v] == l

            if u is None:
                vu = True
            else:
                vu = ubs[v] == u

            vconst = And(vl, vu)
            sameRect.append(vconst)

        # new rectangles cannot overlap, unless on boundaries
        # and should not be equal to an already sampled rectangle

        overlapRect = Or(*overlapRect)
        sameRect = Not(And(*sameRect))

        return And(overlapRect, sameRect)

    def computeSampleVol(self, sample, rmat=None):
        vol = 1

        for v in self.phiVars:
            if self.verb:
                print("VDIST", self.vdist)

            if v not in self.vdist:
                continue

            (l, u) = sample[v]
            vol = vol * self.vdist[v].evalCDF(l, u)

        return vol

    """ Compute an underapproximating volume """

    def getUV(self):

        # map from variables to their upperbounds
        ubs = self.getUbs()

        # map from variables to their lowerbounds
        lbs = self.getLbs()

        # rotate the formula
        mat, gv = self.getRotation(self.rotate)
        #print("ROTATION MATRIX:\n", mat)
        #raw_input("press enter")
        self.phi = self.rotateFormula(mat, gv)

        #print "MADE A ROTATION MATRIX..."
        #wait()
        # Determine whether the rotation was anything
        # other than the identity.
        if np.array_equal(np.abs(1.0*mat),np.identity(mat.shape[0])):
            self.not_identity = False
        else:
            self.not_identity = True

        consts = []
        for v in self.phiVars:
            if v not in self.vdist:
                continue
            if type(self.vdist[v]) is StepVar:
                consts.append(lbs[v] >= min([lb for (lb,ub,pb) in self.vdist[v].bars]))
                consts.append(ubs[v] <= max([ub for (lb,ub,pb) in self.vdist[v].bars]))

        self.phi = And(self.phi, *consts)           
        self.getApproxConstraints(lbs, ubs)

        #wait()
        rects = self.getRectPhi(lbs, ubs)
        #print "GOT RECTANGLES!!!"
        #wait()
        self.sampler.add(rects)

        # sampling and add volume loop
        sampleCount = 0
        volume = 0
        sampleSet = []

        while (True):

            sample = self.getSample(rects, lbs, ubs)
            #print "GOT A SAMPLE!!!"
            #wait()
            # no sample is found, we've found exact volume
            if sample == None:
                #print "NO SAMPLE FOUND!!!!"
                #wait()
                if self.verb:
                    print("computed exact volume")
                break

            self.numSamples += 1

            if self.verb:
                print("\n==> Sample #", sampleCount)
                print("\t", sample)

            # add sample to sample set
            sampleSet.append(sample)

            # negate sample and add it to rects
            negrect = self.negateRect(sample, lbs, ubs)
            if self.verb:
                print("> negrect:", simplify(negrect))

            rects = And(rects, negrect)
            self.sampler.add(negrect)

            if self.verb:
                print("Previous volume", volume)

            sampleVolume = self.computeSampleVol(sample)
            volume = volume + sampleVolume

            sampleCount = sampleCount + 1
            if self.verb:
                print("Sample volume: ", sampleVolume)
                print("Current volume: ", volume)

            if self.adapt:
                self.adaptApprox(sampleVolume)
            #print "SAMPLE VOLUME:", sampleVolume
            #wait()
            yield volume

        if self.verb:
            print("Volume = ", volume)
            print("sampleCount = ", sampleCount)
        # return volume
 
    
    def getDnfUV(self):
        """
        Iteratively compute the volume of a formula by 
        decomposing it into DNF, and separately computing the 
        volume of each disjunct.
        """

        phiprime = self.phi
        dnf_phi = exclusiveToDNF(phiprime, maxlen=50)

        # Some formulae are just ill-conditioned
        # for DNF ("exponential blowup").
        if dnf_phi == []:
            print("DNF PHI:", dnf_phi)
            print("FAILED TO GET REASONABLE DNF. TRYING")
            print("ORDINARY GETUV.")
            gen = self.getUV()
            finished = False
            while True:
                newvol = next(gen)
                yield newvol

            if finished:
                return

        print("DNF:\n", dnf_phi)
        print(len(dnf_phi), " DISJUNCTS")
        dnf_vcomps = [VComp(dj, self.vdist, self.finmax, self.randomize, self.infmax, self.z3qe, 
                            self.adapt, self.rotate, numHists=self.numHists, histBound=self.histBound, 
                            verbose=self.verb,rot_digits=self.rot_digits) for dj in dnf_phi]
        genlist = list([x.getUV() for x in dnf_vcomps])
        volumes = [0 for gen in genlist]
        count = 0
        while True:

            if count == 0:
                for vc in dnf_vcomps:
                    if vc.not_identity == True:
                        self.not_identity = True

            newvols = []
            for i, gen in enumerate(genlist):
                try:
                    newv = next(gen)
                    volumes[i] = newv
                except:
                    continue
                #print "DISJUNCT VOLUME:", volumes[i]

            total_volume = sum(volumes)
            print(volumes, total_volume)
            yield total_volume


    """Get Gaussian Variables"""
    def getGaussians(self):
       
        gvars = np.array(list(self.vdist.keys()))
        gvars = np.array([gv for gv in gvars if type(self.vdist[gv]) is GaussVar])
        gvarnames = [str(gv) for gv in gvars]
        gvars = list(gvars[np.argsort(gvarnames)])

        if self.verb:
            print("GVARS:", gvars)

        return gvars


    def get_atom_coeffs(self, atom):
        """
        Given an atom, extract the coefficients of
        its variables in a dictionary indexed by 
        variable names. 
        """

        atom = simplify(atom)
       
        # Assume that we're given an atomic *inequality*.
        # We want to have "<="; if we have a ">=", we
        # will multiply everything by -1.
        sgn = 1
        if is_le(atom):
            sgn = 1
        elif is_ge(atom):
            sgn = -1
        else:
            #print "ATOM MUST BE A <= or >="
            return {}
        
        # Return a dictionary whose keys are gaussian variable names
        # and values are numbers --- the variables' coefficients.
        coeff_d = {}

        gvars = self.getGaussians()
        gvarnames = [str(var) for var in gvars]

        lhs_coeffs = self.get_lin_coeffs(atom.children()[0])
        rhs_coeffs = self.get_lin_coeffs(atom.children()[1])
        key_union = set(lhs_coeffs.keys()).union(list(rhs_coeffs.keys()))

        for varname in key_union:
            coeff_d[varname] = 0
            if varname in gvarnames:
                if varname in list(lhs_coeffs.keys()):
                    coeff_d[varname] += sgn*lhs_coeffs[varname]
                if varname in list(rhs_coeffs.keys()):
                    coeff_d[varname] += -1*sgn*rhs_coeffs[varname]
            elif varname=='const':
                if varname in list(lhs_coeffs.keys()):
                    coeff_d[varname] += -1*sgn*lhs_coeffs[varname]
                if varname in list(rhs_coeffs.keys()):
                    coeff_d[varname] += sgn*rhs_coeffs[varname]
            else:
                #print "{}\nNOT A USABLE ATOM.".format(atom)
                return {}

        return coeff_d 
      

    """Get coefficients from a Z3 linear combination"""
    def get_lin_coeffs(self,lincom):

        result_dict = {}
        
        # Some code that gets reused a lot:
        def extract_from_product(prod):
            var_inds = []
            coeff = 1
            for i, factor in enumerate(prod.children()):
                if is_rational_value(factor) or is_algebraic_value(factor):
                    coeff = coeff * factor
                else:
                    var_inds.append(i)
            assert len(var_inds) <= 1, "{} must have only one variable.".format(prod)
            value = simplify(coeff)
            value = 1.0*value.numerator_as_long()/value.denominator_as_long()
            return value, prod.children()[var_inds[0]] 

        result_dict['const'] = 0

        #print "LINEAR COMBINATION", lincom
       
        if is_rational_value(lincom):
            value = simplify(lincom)
            value = 1.0*value.numerator_as_long()/value.denominator_as_long()
            result_dict['const'] += value 
        elif is_mul(lincom):
            coeff, var = extract_from_product(lincom)
            result_dict[str(var)] = coeff
        elif is_add(lincom):
            for child in lincom.children():
                ch = simplify(child)
                if is_mul(child):
                    coeff, var = extract_from_product(child)
                    result_dict[str(var)] = coeff
                elif is_rational_value(child) or is_algebraic_value(child):
                    result_dict['const'] += child
                else:
                    result_dict[str(child)] = 1
        elif lincom.decl().kind() == Z3_OP_UNINTERPRETED:
            result_dict[str(lincom)] = 1

        return result_dict


    """Get a rotation matrix that is well-suited for this formula"""
    def getRotation(self,heuristic='identity'):
       
        if heuristic not in ['nonzero-entries',
                             'surface-integral',
                             'nearest',
                             'identity',
                             'first-sample',
                             'forced']:
            heuristic = 'identity'

        if heuristic=='identity':
            rel_gvars = self.getGaussians()
            rmat = np.identity(len(rel_gvars),dtype=int)
            return rmat, rel_gvars 

        #preds = list(getPreds(self.phi)) 
        #if self.verb:
        #    print "PREDS\n", preds
        #atomcoeffs = [self.getAtomCoeffs(pred) for pred in preds]
      
        # This encodes the linear constraints in
        # a more numerically palatable fashion.
        numphi = numericalPhi(self.phi,valid_vars = self.getGaussians()) 

        # Select the face we wish to align with
        # by one of the available heuristics.
        if heuristic=='nonzero-entries':
            rel_gvars, rel_coeffs = self.selectByNonzero(numphi)
            rel_coeffs = [rel_coeffs]

        elif heuristic=='surface-integral':
            rel_gvars, rel_coeffs = self.selectByIntegral(numphi)

        elif heuristic=='nearest':
            rel_gvars, rel_coeffs = self.selectByNonzero(numphi)
            rel_coeffs = [rel_coeffs]

        elif heuristic=='forced':
            rel_gvars, rel_coeffs = self.rel_gvars, self.rel_coeffs

        elif heuristic=='first-sample':
            rel_gvars, rel_coeffs = self.selectByFirstSample()

        if len(rel_gvars) == 0:
            rel_gvars = self.getGaussians()
            rmat = np.identity(len(rel_gvars),dtype=int)
            return rmat, rel_gvars 

        #print "RELEVANT COEFFICIENTS:", rel_coeffs
        #print "RELEVANT GAUSSIAN VARS:", rel_gvars

        # Return a rotation matrix. Let it be the identity
        # if we aren't really rotating.
        #rmat = rm.rh_align(rel_coeffs)
        #rmat = rm.pairwise_rh_align(rel_coeffs)
        rmat = rm.full_rational_align(rel_coeffs,precision=self.rot_digits)

        if self.verb:
            print("ROTATED VARS:", rel_gvars)
            print("COEFFS:", rel_coeffs)
            print("RMAT\n", 1.0*rmat)

        return rmat, rel_gvars 



    def selectByNonzero(self,numphi):
        """
        Given a list of atom data (each is a dict mapping
        variables to coefficients), select one of the 
        atoms and return its variables and coefficients
        in a canonical order.
        This method selects a dictionary based on how many 
        nonzero entries it has.
        """
        A = numphi.constraint_matrix
        if self.verb:
            print("PHI CONSTRAINTS:", A) 
        atomlens = [len(A[row,A[row,:]!=0]) for row in range(A.shape[0])]
        biggest = A[np.argmax(atomlens),:]
        if self.verb:
            print("BIGGEST", biggest)

        rel_gvars = self.getGaussians() 

        return rel_gvars, rel_coeffs


    def selectByFirstSample(self):
        """
        Try each possible rotation for one iteration, then
        stick with the one that does best
        on that one iteration.
        """
        ordered_gvars = self.getGaussians()

        preds = list(getPreds(self.phi)) 
        #print "PREDS\n", preds
        atomcoeffs = [self.get_atom_coeffs(pred) for pred in preds]
        to_try = []
        for d in atomcoeffs:
            if len(list(d.keys())) > 2:
                coeffs = [d[str(gv)] for gv in ordered_gvars if str(gv) in list(d.keys())]
                redundant = False
                for already in to_try:
                    already_coeffs = [already[str(gv)] for gv in ordered_gvars if str(gv) in list(already.keys())]
                    if len(already_coeffs)==len(coeffs)\
                        and set(already.keys()) == set(d.keys())\
                        and abs(np.dot(already_coeffs,coeffs)) < 1e-9:
                        redundant = True
                if not redundant:
                    to_try.append(d)

        first_samples = []
        rel_gvars_list = [[ordered_gvars[0]]]
        rel_coeffs_list = [[1]]
        #print "{} TO TRY;".format(len(to_try))
        for i, d in enumerate(to_try):

            rel_gvars = [gv for gv in ordered_gvars if str(gv) in list(d.keys())]
            rel_coeffs = [d[str(gv)] for gv in rel_gvars]

            #print "TRYING {} OF {}".format(i+1, len(to_try))
            subvcomp = VComp(self.phi, self.vdist, self.finmax, self.randomize, self.infmax, 
                             self.z3qe, self.adapt, 'forced', numHists = self.numHists, 
                             histBound=self.histBound, verbose=self.verb, rot_digits=self.rot_digits)

            subvcomp.rel_gvars = rel_gvars
            subvcomp.rel_coeffs = [rel_coeffs]

            sub_uv = subvcomp.getUV()
            first_samples.append(next(sub_uv))
            rel_gvars_list.append(rel_gvars)
            rel_coeffs_list.append(rel_coeffs)

        best_ind = int(np.argmax(first_samples))
        #print first_samples
        #print "{} WAS BEST".format(best_ind)
        #wait()

        return  rel_gvars_list[best_ind], [rel_coeffs_list[best_ind]]


    def selectByIntegral(self, numphi):
        """
        Given a list of atom data (each is a dict mapping
        variables to coefficients), select one of the 
        atoms and return its variables and coefficients
        in a canonical order.
        Each atom corresponds to a face of the region; this
        method computes a surface integral of the joint pdf
        over each of the faces and returns the one with largest
        computed integral.
        """
        ordered_gvars = self.getGaussians()

        A = numphi.constraint_matrix
        b = numphi.constraint_b
        d = numphi.disjuncts

        if A.shape[0] == 0:
            return [], []

        rankings = rm.mc_gauss_surface_integrals(A, b, d,
                                                 mc_sample=10000)

        #print "RANKINGS", rankings
        ordered_inds = np.argsort(rankings)[::-1][:len(rankings)]
        ordered_coeffs = A[ordered_inds,:]

        sums = np.sum(np.abs(ordered_coeffs),axis=0)
        rel_coeffs = ordered_coeffs[:,sums != 0]
        rel_gvars = [ogv for i, ogv in enumerate(ordered_gvars) if sums[i] != 0]
            
        assert rel_coeffs.shape[1] == len(rel_gvars)

        return rel_gvars, rel_coeffs


    """Return a version of phi whose (gaussian) free variables have been rotated."""
    def rotateFormula(self, rmat, rvars):

        if len(rvars) == 0:
            return self.phi 
            
        assert (len(rmat.shape)==2 and rmat.shape[0]==len(rvars) \
            and rmat.shape[1]==len(rvars)),\
            "rmat must be a square matrix with \
             size corresponding to gaussian variables."

        subslist = []
        # Do the transformation: we will replace
        # each gaussian variable with a linear combination of them;
        # the inner product of the variables with a row of the rotation matrix.
        for i, gv1 in enumerate(rvars):
            expression = 0
            for k, gv2 in enumerate(rvars):
                expression = expression + rmat[i,k]*gv2
            #print "RESULTING COEFFICIENT:", simplify(expression)
            subslist.append((gv1,simplify(expression)))

        result = substitute(self.phi, *subslist)
        result = simplify(result)

        return result

# x = Real('x')
# y = Real('y')
# dists = {}
# dists[x] = ('G', 0, 1)
#dists[y] = ('G', 12, 5)
# phi = And(x >= 0, y == x)
# #phi = And(x**2 > y)
# #phi = And(x > y)
# phi = True
# prev = None
# for i in range(0, 2):
#     x = Real('x' + str(i))
#     dists[x] = ('G', i, i + 1)
#     if prev is not None:
#         phi = And(phi,  x + prev <= 10)
#     else:
#         phi = And(phi, x >= i)
#     prev = x
# phi = simplify(phi)
# print phi
# vc = VComp(phi, dists, square=False, hsample=False,
#            pointRect=False, optimize=True)
#
# v = vc.getUV()
# for i in range(1,2):
#     print v.next()
# vc.getUV(200)
# print vc.integrate(phi)
