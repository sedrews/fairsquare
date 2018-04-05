import ast
from asteval import Interpreter
import numpy.random as np
import random


def gaussian(m, v):
    return np.normal(m, v ** 0.5)

def uniform(l, u):
    assert( u > l )
    return np.uniform(l, u)

def step_(massDist):
    randRoll = random.random()  # in [0,1)
    sum = 0
    result = 1
    for mass in massDist:
        sum += mass
        if randRoll < sum:
            return result
        result += 1

    return result


def step(l):
    ps = [p[2] for p in l]
    i = step_(ps)
    res = float(l[i - 1][1]) - \
        ((float(l[i - 1][1]) - float(l[i - 1][0])) / float(2))
    return res

m = False
h = False


def sensitiveAttribute(b):
    global m
    m = b


def fairnessTarget(b):
    global h
    h = b


def simulate(fn, check, times):
    f = open(fn, "r")
    p = ""

    for line in f:

        if ("def " in line) and ("popModel" not in line):
            continue

        if "return" in line:
            continue

        p += line

    p += "    return\n"
    p += "res = popModel()\n"
    print(p)

    tree = ast.parse(p)

    mcount = 0
    hcount = 0
    mhcount = 0
    notmhcount = 0
    for counter in range(1, int(times)):
        global m
        global h
        exec(compile(tree, filename="<ast>", mode="exec"))
        if m:
            mcount += 1
        if h:
            hcount += 1
        if m and h:
            mhcount += 1
        if (not m) and  h:
            notmhcount += 1

        m = False
        h = False

    print("\n\nSIMULATION REPORT:\n")
    print("P(M) = ", float(mcount) / float(times))
    print("P(MH) = ", float(mhcount) / float(times))
    print("P(H) = ", float(hcount) / float(times))
    print("P(not M and H) = ", float(notmhcount) / float(times))
