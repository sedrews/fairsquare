from z3 import *
from z3extra import *
from subprocess import Popen, PIPE
import ast
from parse import Encoder
import platform
from itertools import tee, zip_longest
from functools import reduce

def pairwise(iterable):
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def z3toRedlog(vars, phi, exists):
    # map from z3 to redlog
    m = {}

    def infix(e, op):
        ops = list(map(z3toRedlog_, e.children()))

        res = "(" + reduce(lambda v, v1: v + " " + op + " " + v1, ops) + ")"

        m[e] = res

        return res

    def z3toRedlog_(e):
        if e in m:
            return m[e]

        if is_not(e):
            assert(len(e.children()) == 1)
            arg = e.children()[0]
            res = "(not (" + z3toRedlog_(arg) + "))"
            m[e] = res

            return res

        if is_and(e):
            return infix(e, "and")

        if is_or(e):
            return infix(e, "or")

        if is_add(e):
            return infix(e, "+")

        if is_sub(e):
            return infix(e, "-")

        if is_mul(e):
            return infix(e, "*")

        # not in Z3's Python API
        if is_app_of(e, Z3_OP_POWER):
            return infix(e, "^")

        #if is_div(e):
        #    return infix(e, "/")

        if is_gt(e):
            return infix(e, ">")

        if is_ge(e):
            return infix(e, ">=")

        if is_lt(e):
            return infix(e, "<")

        if is_le(e):
            return infix(e, "<=")

        if is_eq(e):
            return infix(e, "=")

        if is_rational_value(e):
            return str(e)

        if is_const(e):
            return str(e)

        assert(False)

    phinnf = (Tactic('nnf')(phi)).as_expr()
    rvars = []

    for v in vars:
        m[v] = str(v)
        rvars.append(m[v])

    rvarStr = "{" + reduce(lambda v, v1: v + ", " + v1, rvars) + "}"

    if exists:
        redForm = "ex(" + rvarStr + ", " + z3toRedlog_(phinnf) + ")"
    else:
        redForm = "all(" + rvarStr + ", " + z3toRedlog_(phinnf) + ")"

    return redForm

# fixes the way powers are weirdly printed by redlog
def fixPowers(s):

    def merge(a,b):
        if a == " ":
            return b
        assert(b==" ")
        return "**" + a

    ls = s.split('\n')
    res = ""
    lastPower = False
    for l1,l2 in pairwise(ls):
        print("l1: ", l1)
        print("l2: ", l2)
        if lastPower:
            lastPower = False
            continue
        #remove white space and check if digit
        if ("".join(l1.split())).isdigit():
            z = list(zip_longest(l1,l2,fillvalue=" "))
            merged = "".join([merge(a_b[0],a_b[1]) for a_b in z])
            res += merged
            print("merged", merged)
            lastPower = True
        else:
            res += l1
            lastPower = False

    print(lastPower)
    if not lastPower: res += l2

    return res

def qelimRedlog(vars, phi, exists=False):
    if len(vars) == 0:
        return phi

    s = ""
    s += "rlset reals;\n"
    s += "phi := " + z3toRedlog(vars, phi, exists) + ";\n"
    s += "out \"" + "qelim.res" + "\";\n"
    s += "rlqe phi;\n"
    s += "shut \"" + "qelim.res" + "\";\n"
    s += "quit;\n"

    # remove temp.red
    try:
        os.remove("temp.red")
    except OSError:
        pass

    # remove qelim.res
    try:
        os.remove("qelim.res")
    except OSError:
        pass

    # create and write to problem file
    f = open("temp.red", "w+")
    f.write(s)
    f.close()

    # create result file
    resfile = open("qelim.res", "w+")
    resfile.close()

    p = Popen(['./tools/reduce', "temp.red"], stdout=PIPE)
    p.wait()

    # open and read results file
    resfile = open("qelim.res", "r")
    resstr = resfile.read()

    # HACK: fixPowers assumes powers < 10
    # uncomment at your own risk
    # resstr = fixPowers(resstr)

    resstr = resstr.replace("\n", " ")
    resstr = resstr.lstrip()
    resstr = resstr.replace(" = ", " == ")
    resstr = resstr.replace(" <> ", " != ")
    resstr = resstr.replace("true", "True")
    resstr = resstr.replace("false", "False")
    resstr = resstr.replace("^", "**")


    resfile.close()

    return redlogToZ3(resstr)


def redlogToZ3(s):
    e = Encoder()
    past = ast.parse(s).body
    assert len(past) == 1

    z3exp = e.exprToZ3(past[0].value)

    return z3exp


# x = Real('x')
# y = Real('y')
# phi = Implies(y > x, x == 1)
# res = qelimRedlog([y], phi)
# print res
