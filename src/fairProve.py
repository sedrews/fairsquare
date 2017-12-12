from parse import Encoder
import argparse
import ast
from vol import VComp
from z3 import *
from z3extra import *
from threading import Thread
from redlogInterface import *
from simulate import *
import time
from probcomp import ProbComp
from logwriter import LogWriter
from plotter import Plotter
from functools import reduce

#a global variable
globalqelimtime = 0

def parseArgs():

    parser = argparse.ArgumentParser(
        description='''\u25B6\u25B6\u25B6 The Glorious Madison Algorithmic Fairness Prover.''')

    parser.add_argument("-f", "--file", required=True,
                        help="Input python bench")

    parser.add_argument("-e", "--epsilon", required=False, default=.1,
                        help="Epsilon fairness guarantee")

    parser.add_argument("-s", "--simulate", required=False, action='store_true', default=False,
                        help="Simulate program to check property")

    parser.add_argument("-sc", "--simTimes", required=False, default=1000,
                        help="Number of simulation times")

    parser.add_argument("-mf", "--finiteMaximize", required=False, action='store_true', default=False,
                        help="attempt to maximize the finite bounds of samples")

    parser.add_argument("-r", "--randomize", required=False, default=None,
                        help="If using finiteMaximize, randomly permutes the order of the variables to try maximizing using the given random seed")

    parser.add_argument("-mi", "--infiniteMaximize", required=False, action='store_true', default=False,
                        help="attempt to drop bounds for unbounded sides of samples")

    parser.add_argument("-p", "--plot", required=False, action='store_true', default=False,
                        help="Plots progress of sampling")

    parser.add_argument("-z", "--z3qe", required=False, action='store_true', default=False,
                        help="Use z3 quantifier elimination instead of redlog")

    parser.add_argument("-nh", "--numHists", required=False, default=5,
                        help="Number of buckets to use in histogram approximations (preferably odd)")

    parser.add_argument("-hb", "--histBound", required=False, default=3,
                        help="Sets the histogram bound to +/- variance times this number")

    parser.add_argument("-o", "--output", required=False, default=None,
                        help="File to write results")

    parser.add_argument("-t", "--timeout", required=False, default=None,
                        help="Use a timeout period in seconds after which to terminate (not including time for qelim) instead of an epsilon") 

    parser.add_argument("-a", "--adapt", required=False, action='store_true', default=False,
                        help="Dynamically increase the histogram bound for approximations when samples are small (in a rudimentary way")

    parser.add_argument("-ro", "--rotate", required=False, default='identity',
                        help="Enable rotations of the formula for improved convergence in some circumstances")

    parser.add_argument("-v", "--verbose", required=False, action='store_true', default=False,
                        help="Output a lot")

    return parser.parse_args()


def projectNonProbVars(phi, pvars, z3qe):

    allvars = get_vars(phi)

    # nonProbVars = allvars \cap \not pvars
    nonProbVars = list(set(allvars) - set(pvars))

    # all vars probabilistic
    if len(nonProbVars) == 0:
        return phi

    qelimstart = time.time()
    if z3qe:
        #res = qelimExists(nonProbVars, phi)
        res = Not(qelimForallDNF(nonProbVars, Not(phi)))
    else:
        res = qelimRedlog(nonProbVars, phi, exists=True)
    global globalqelimtime
    globalqelimtime += time.time() - qelimstart

    return res


def proveFairness(p, output, epsilon, finmax, randomize, infmax, plot, z3qe, numHists, histBound, timeout, adapt, rotate, verbose):

    start = time.time()

    noqual = True
    if p.qualified is not None:
        s = Solver()
        s.add(Not(p.qualified))
        if s.check() == sat: #if qual is tautologically true, don't bother
            noqual = False
    
    if noqual:
        M = And(p.model, p.sensitiveAttribute)
        MH = And(p.model, p.program, p.sensitiveAttribute, p.fairnessTarget)
        notMH = And(p.model, p.program, Not(p.sensitiveAttribute), p.fairnessTarget)
        phis = [M, MH, notMH]
        phinames = ["m", "mh", "nmh"]
    else:
        MQ = And(p.model, p.sensitiveAttribute, p.qualified)
        notMQ = And(p.model, Not(p.sensitiveAttribute), p.qualified)
        MQH = And(p.model, p.program, p.sensitiveAttribute, p.fairnessTarget, p.qualified)
        notMQH = And(p.model, p.program, Not(p.sensitiveAttribute), p.fairnessTarget, p.qualified)
        phis = [MQ, notMQ, MQH, notMQH]
        phinames = ["mq", "nmq", "mqh", "nmqh"]
    
    # all probabilistic vars
    pvars = [x for x in p.vdist]

    # existentially eliminate all non-probabilistic (e.g., SSA) variables
    print("running quantifier elimination (this may take some time) ...\n") # univ qe happens later though
    phis = [projectNonProbVars(phi, pvars, z3qe) for phi in phis]

    vcargs = [p.vdist, finmax, randomize, infmax, z3qe, adapt, rotate, numHists, histBound, verbose]

    probs = [ProbComp(name, phi, *vcargs) for name, phi in zip(phinames,phis)]

    def gratiol(PrM, PrMH, PrnotMH):
        #underapproximate the numerator, P(H|M) = P(HM) / P(M)
        phgm = PrMH.under() / PrM.over()
        #overapproximate the denominator, P(H|!M) = P(H!M) / P(!M)
        phgnm = PrnotMH.over() / (1 - PrM.over())
        return phgm / phgnm
    def gratiou(PrM, PrMH, PrnotMH): 
        #overapproximate the numerator
        phgm = PrMH.over() / PrM.under()
        #underapproximate the denominator
        phgnm = PrnotMH.under() / (1 - PrM.under())
        return phgm / phgnm
    def gpost(PrM, PrMH, PrnotMH, epsilon):
        try: l = gratiol(PrM, PrMH, PrnotMH)
        except ZeroDivisionError: l = None
        try: u = gratiou(PrM, PrMH, PrnotMH)
        except ZeroDivisionError: l = None
        print("ratio bounds:", float(l) if l is not None else "None", float(u) if u is not None else "None")
        if u is not None and u < 1 - epsilon:
            return False
        elif l is not None and l > 1 - epsilon:
            return True
        else:
            return None

    def gqratiol(PrMQ, PrnotMQ, PrMQH, PrnotMQH):
        #underapproximate the numerator, P(H|MQ) = P(HMQ) / P(MQ)
        phgmq = PrMQH.under() / PrMQ.over()
        #overapproximate the denominator, P(H|!MQ) = P(H!MQ) / P(!MQ)
        phgnmq = PrnotMQH.over() / PrnotMQ.under()
        return phgmq / phgnmq
    def gqratiou(PrMQ, PrnotMQ, PrMQH, PrnotMQH):
        #overapproximate the numerator
        phgmq = PrMQH.over() / PrMQ.under()
        #underapproximate the denominator
        phgnmq = PrnotMQH.under() / PrnotMQ.over()
        return phgmq / phgnmq
    def gqpost(PrMQ, PrnotMQ, PrMQH, PrnotMQH, epsilon):
        try: l = gqratiol(PrMQ, PrnotMQ, PrMQH, PrnotMQH)
        except ZeroDivisionError: l = None
        try: u = gqratiou(PrMQ, PrnotMQ, PrMQH, PrnotMQH)
        except ZeroDivisionError: u = None
        print("ratio bounds:", float(l), float(u))
        if u is not None and u < 1 - epsilon:
            return False
        elif l is not None and l > 1 - epsilon:
            return True
        else:
            return None

    if noqual:
        ratiol = lambda : gratiol(*probs)
        ratiou = lambda : gratiou(*probs)
        post = lambda : gpost(*probs, epsilon)
    else:
        ratiol = lambda : gqratiol(*probs)
        ratiou = lambda : gqratiou(*probs)
        post = lambda : gqpost(*probs, epsilon)

    proof = False


    qelimtime = globalqelimtime + reduce(lambda x,y : x + y, [t.qtime() for t in probs])
    
    if plot:
        plotter = Plotter(epsilon)

    if output is not None:
        def csvf(time, *probs):
            return [str(x) for x in [time, ratiol(), ratiou(), *reduce(lambda y,z : y + z , [[pr.under(), pr.over()] for pr in probs])]]
        log = LogWriter(output, ["time","ratiol","ratiou", *reduce(lambda x,y: x + y, [[pr.name + "l", pr.name + "u"] for pr in probs])])
        log.message("using " + ("z3qe" if z3qe else "rlqe"))
        log.message("computing probs: " + reduce(lambda x,y : x + " " + y, [pr.name for pr in probs]))

    if timeout is not None:
        loop = lambda e,t: t
    else:
        loop = lambda e,t: e

    alldone = False

    while loop(not proof, time.time() - start - qelimtime < timeout if timeout is not None else -1):

        alldone = True
        
        for Pr in probs:
            if not Pr.exact():
                try:
                    Pr.refine()
                    alldone = False
                except StopIteration:
                    pass

        qelimtime = globalqelimtime + reduce(lambda x,y : x + y, [t.qtime() for t in probs])
        
        res = post()
        if res is not None:
            proof = True
            racist = not res

        if alldone:
            break

        if plot:
            plotter.draw(ratiol(), ratiou())
        if output is not None:
            log.data(csvf(time.time() - start - qelimtime, *probs))

    end = time.time()

    msg = ""
    if timeout is None or alldone:
        if racist:
            msg += "postcondition does not hold: UNFAIR\n"
        else:
            msg += "postcondition holds: FAIR\n"
    msg += "samples computed: " + str(reduce(lambda x,y: x + y, [p.numSamples() for p in probs])) + "\n"
    msg += "sampling run time: " + str(end - start - qelimtime) + "\n"
    msg += "qelim time: " + str(qelimtime)

    print(msg)

    if output is not None:
        log.message(msg)
        log.close()

def main():

    args = parseArgs()

    if args.simulate:
        simulate(args.file, args.simulate, args.simTimes)
        return

    f = open(args.file, "r")
    node = ast.parse(f.read())

    e = Encoder()
    e.visit(node)

    print("\n\n== Population and program encoding == ")
    print("Population model: ", e.model)
    print("Program: ", e.program)

    print("Probability dists: ", e.vdist)
    print("fairness target", e.fairnessTarget)
    print("sensitiveAttribute", e.sensitiveAttribute)
    print("\n\n")


    #print("ROTATE COMMAND LINE ARGUMENT:\n", args.rotate)

    timeoutarg = int(args.timeout) if args.timeout is not None else None
    randarg = int(args.randomize) if args.randomize is not None else None
    proveFairness(e, args.output, float(args.epsilon), args.finiteMaximize, randarg, args.infiniteMaximize, 
            args.plot, args.z3qe, int(args.numHists), int(args.histBound), timeoutarg, args.adapt,
            args.rotate, args.verbose)

if __name__ == "__main__":
    main()
