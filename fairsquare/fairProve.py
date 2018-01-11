from .parse import parse_fr
import ast
from .vol import VComp
from z3 import *
from .z3extra import *
from threading import Thread
from .redlogInterface import *
from .simulate import *
import time
from .probcomp import ProbComp
from .logwriter import LogWriter
from .plotter import Plotter
from functools import reduce

#a global variable
globalqelimtime = 0

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

def proveFairness(source, output, epsilon, finmax, randomize, infmax, plot, z3qe, numHists, histBound, timeout, adapt, rotate, verbose):
#def proveFairness(source, ADF_params, vol_params):

    p = parse_fr(source)

    print("\n\n== Population and program encoding == ")
    print("Population model: ", p.model)
    print("Program: ", p.program)

    print("Probability dists: ", p.vdist)
    print("fairness target", p.fairnessTarget)
    print("sensitiveAttribute", p.sensitiveAttribute)
    print("qualified", p.qualified)

    #print("ROTATE COMMAND LINE ARGUMENT:\n", args.rotate)

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

    def gpost(PrM, PrMH, PrnotMH, epsilon):
        ratio = None
        res = None
        try:
            PrHgM = PrMH / PrM
            PrHgnotM = PrnotMH / (1 - PrM)
            ratio = PrHgM / PrHgnotM
            print("ratio bounds:", str(float(ratio.lower)) + " " + str(float(ratio.upper)) if ratio is not None else "None None")
            res = bool(ratio > 1-epsilon) # Will raise ValueError if unknown
        except (ZeroDivisionError, ValueError):
            pass
        return res
    def gqpost(PrMQ, PrnotMQ, PrMQH, PrnotMQH, epsilon):
        ratio = None
        res = None
        try:
            PrHgMQ = PrMQH / PrMQ
            PrHgnotMQ = PrnotMQH / PrnotMQ
            ratio = PrHgMQ / PrHgnotMQ
            print("ratio bounds:", str(float(ratio.lower)) + " " + str(float(ratio.upper)) if ratio is not None else "None None")
            res = bool(ratio > 1-epsilon)
        except (ZeroDivisionError, ValueError):
            pass
        return res

    if noqual:
        post = lambda : gpost(*[p.value for p in probs], epsilon)
    else:
        post = lambda : gqpost(*[p.value for p in probs], epsilon)

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

