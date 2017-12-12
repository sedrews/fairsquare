import ast
import codegen
from z3 import *
from z3extra import *
import rotationmath as rm
# helpers


def name(node):
    return ast.dump(node).split("(")[0]

# HACK:


def isTrue(node):
    return "True" in ast.dump(node)
# HACK:


def isFalse(node):
    return "False" in ast.dump(node)


def isCall(node):
    return name(node) == 'Call'


def isAssign(node):
    return name(node) == 'Assign'


def isIf(node):
    return name(node) == 'If'


def isExpr(node):
    return name(node) == 'Expr'


def isReturn(node):
    return name(node) == 'Return'


def isBinOp(node):
    return name(node) == 'BinOp'


def isBoolOp(node):
    return name(node) == 'BoolOp'


def isUnaryOp(node):
    return name(node) == 'UnaryOp'


def isCompareOp(node):
    return name(node) == 'Compare'


def isAdd(node):
    return name(node) == 'Add'


def isPow(node):
    return name(node) == 'Pow'


def isAnd(node):
    return name(node) == 'And'


def isOr(node):
    return name(node) == 'Or'


def isUAdd(node):
    return name(node) == 'UAdd'


def isUSub(node):
    return name(node) == 'USub'


def isMult(node):
    return name(node) == 'Mult'


def isDiv(node):
    return name(node) == 'Div'

def isSub(node):
    return name(node) == 'Sub'


def isNot(node):
    return name(node) == 'Not'


def isEq(node):
    return name(node) == 'Eq'


def isNotEq(node):
    return name(node) == 'NotEq'

#Lt | LtE | Gt | GtE


def isLt(node):
    return name(node) == 'Lt'


def isLtE(node):
    return name(node) == 'LtE'


def isGt(node):
    return name(node) == 'Gt'


def isGtE(node):
    return name(node) == 'GtE'


def isNum(node):
    return name(node) == 'Num'


def isName(node):
    return name(node) == 'Name'


def evalAST(node):
    node = codegen.to_source(node)
    return eval(node)


def makeBin(op, l, r):
    if isAdd(op):
        return l + r
    elif isSub(op):
        return l - r
    elif isMult(op):
        return l * r
    elif isDiv(op):
        return l / r
    elif isPow(op):
        return l**r
    else:
        assert(False and "Weird binary op")


def makeUnary(op, e):
    if isUSub(op):
        return -e
    elif isUAdd(op):
        return +e
    elif isNot(op):
        return Not(e)
    else:
        assert(False and "Weird unary op")


def makeBool(op, *args):
    if isAnd(op):
        return And(*args)
    elif isOr(op):
        return Or(*args)
    else:
        assert(False and "Weird bool op")


def makeCompare(op, l, r):
    if isGt(op):
        return l > r
    elif isGtE(op):
        return l >= r
    elif isLt(op):
        return l < r
    elif isLtE(op):
        return l <= r
    elif isEq(op):
        return l == r
    elif isNotEq(op):
        return l != r
    else:
        assert(False)


class Encoder(ast.NodeVisitor):

    def generic_visit(self, node):
        return

    def visit_Module(self, node):
        fpop = node.body[0]
        assert (fpop.name == 'popModel')

        fprog = node.body[1]

        d = initDict(node)

        # place holder for max values
        self.gd = d.copy()

        self.model = self.doSeq(fpop.body, d)
        self.model = simplify(self.model)

        self.program = self.doSeq(fprog.body, d)
        self.progam = simplify(self.program)

    def __init__(self):
        self.vdist = {}

        self.sensitiveAttribute = None
        self.qualified = None
        self.fairnessTarget = None

        self.model = None
        self.program = None

    def doSeq(self, seq, d):
        trans = []
        for s in seq:
            #print "ENCODING", ast.dump(s)[0:200]
            if isAssign(s):
                t = self.visit_Assign(s, d)
                trans.append(t)

            elif isIf(s):
                t = self.visit_If(s, d)
                trans.append(t)

            elif isExpr(s):
                # this is a call with no lhs
                t = self.visit_Call(s.value, d)
                trans.append(t)

            elif isReturn(s):
                # do nothing on return
                continue

            else:
                assert(False)

        res = simplify(And(*trans))

        # print "\n----------"
        # print "res: ", res
        # print "sres: ", simplify(res)

        return res

    def probAssign(self, lhs, rhs, d):

        self.refresh(d, lhs)
        zlhs = var(lhs, d[lhs])

        fname = rhs.func.id
        
        assert(zlhs not in self.vdist)

        if fname == 'gaussian':
            mean = evalAST(rhs.args[0])
            variance = evalAST(rhs.args[1])
            std_dev = RealVal(rm.sigfig_str(variance ** 0.5, 4))

            assert(variance >= 0)

            # We normalize all gaussian variables,
            # and let the program shift and scale them.
            #self.vdist[zlhs] = ('G', mean, variance)
          
            self.refresh(d, lhs)
            zlhs2 = var(lhs, d[lhs])
            self.refresh(d, lhs)
            zlhs3 = var(lhs, d[lhs])
            
            self.vdist[zlhs] = ('G', 0, 1)

            res = And(zlhs2 == zlhs*std_dev, zlhs3 == zlhs2+mean)

            return res


        elif fname == 'step':
            step = ('S', evalAST(rhs.args[0]))
            self.vdist[zlhs] = step

            # sum of probs == 1
            l = evalAST(rhs.args[0])
            s = sum([a[2] for a in l])
            
            assert(abs(s-1.0) <= 0.00001)

            lbounds = [a_b_c[0] for a_b_c in step[1]]
            ubounds = [a_b_c1[1] for a_b_c1 in step[1]]

            res = And(zlhs >= min(lbounds), zlhs <= max(ubounds))

            return res

        else:
            assert(False and "Supported distributions: gaussian")

        return True

    def exprToZ3(self, e, d=None):
        #print "we are here ", ast.dump(e)
        if isBinOp(e):
            op = e.op
            zlhs = self.exprToZ3(e.left, d)
            zrhs = self.exprToZ3(e.right, d)
            return makeBin(op, zlhs, zrhs)

        elif isUnaryOp(e):
            op = e.op
            zexp = self.exprToZ3(e.operand, d)
            return makeUnary(op, zexp)

        elif isBoolOp(e):
            op = e.op
            zexprs = [self.exprToZ3(v, d) for v in e.values]
            return makeBool(e.op, zexprs)

        elif isCompareOp(e):
            assert(len(e.ops) == 1)

            op = e.ops[0]
            zlhs = self.exprToZ3(e.left, d)
            zrhs = self.exprToZ3(e.comparators[0], d)
            return makeCompare(op, zlhs, zrhs)

        elif isNum(e):
            num = evalAST(e)
            return num

        elif isName(e):
            # HACK:
            if isTrue(e):
                return And(True)
            if isFalse(e):
                return And(False)
            
            if d is None:
                n = var(e.id)
            else:
                n = var(e.id, d[e.id])
            return n

        else:
            assert(False and ("Weird expression" + ast.dump(e)))

    def refresh(self, d, v):
        # print "--------->", d[v], self.gd[v]
        d[v] = max(d[v], self.gd[v]) + 1
        self.gd[v] = d[v]

    def visit_Assign(self, node, d):
        #print "Processing Assign"
        #print "Assign", ast.dump(node)[0:10]
        #print "POPULATION MODEL: ", self.model

        assert(len(node.targets) == 1)

        lhs = node.targets[0].id
        rhs = node.value

        # increment SSA ID

        if isCall(rhs):
            return self.probAssign(lhs, rhs, d)

        zrhs = self.exprToZ3(rhs, d)
        self.refresh(d, lhs)
        zlhs = var(lhs, d[lhs])
        res = zlhs == zrhs

        return res

    def visit_Call(self, node, d):
        # print "Processing Call"
        # print "Call", ast.dump(node)[0:10]

        fn = node.func.id
        assert(len(node.args) == 1)

        if fn == 'assume':
            res = self.exprToZ3(node.args[0], d)
            return res

        elif fn == 'sensitiveAttribute':
            self.sensitiveAttribute = self.exprToZ3(node.args[0], d)
            return True

        elif fn == 'qualified':
            self.qualified = self.exprToZ3(node.args[0], d)
            return True

        elif fn == 'fairnessTarget':
            self.fairnessTarget = self.exprToZ3(node.args[0], d)
            return True

        else:
            assert(False and "Unrecognizable function call")

    def createPhiNode(self, d, dt):
        # then branch constraints
        phiT = []

        # else branch constraints
        phiE = []

        for v in d:
            # check if vars have different subscripts
            # on the two branches, if so, refresh subscript
            # and add to phi nodes
            if d[v] != dt[v]:
                i = d[v]
                ev = var(v, d[v])

                it = dt[v]
                tv = var(v, dt[v])

                self.refresh(d, v)

                ni = d[v]
                nv = var(v, ni)

                phiT.append(nv == var(v, it))
                phiE.append(nv == var(v, i))

        # print phiT
        # print phiE
        return (And(*phiT), And(*phiE))

    def visit_If(self, node, d):
        #print "Processing If"
        #print "If", ast.dump(node)[0:10]
        
        zcond = self.exprToZ3(node.test, d)

        dt = d.copy()
        zthen = self.doSeq(node.body, dt)
        zelse = self.doSeq(node.orelse, d)

        (zphiT, zphiE) = self.createPhiNode(d, dt)

        #print "COND: ", zcond
        #print "ZT: ", zthen
        #print "ZELSE: ", zelse

        resT = Implies(zcond, And(zthen, zphiT))
        resE = Implies(Not(zcond), And(zelse, zphiE))
        res = And(resT, resE)

        #res = And(If(zcond, And(zthen, zphiT), And(zelse, zphiE)))

        return res


class VDict(ast.NodeVisitor):

    def __init__(self):
        self.vdict = {}

    def visit_Name(self, node):
        # print ast.dump(node)
        self.vdict[node.id] = 0


def initDict(node):
    v = VDict()
    v.visit(node)
    return v.vdict

if __name__=="__main__":


    from runTests import projectNonProbVars
    node = ast.parse('''
def popModel():
    m = gaussian(0,100)
    rank = gaussian(25,25)
    exp = gaussian(10,9)

    if m < 10:
        rank = rank + 5

    sensitiveAttribute(m > 10)

def F():
    t = 1.0 * m + 3.0 * rank + -1*exp
    fairnessTarget(rank < 5 or exp-rank > -5)
    return t
    ''')
    e = Encoder()
    e.visit(node)
    print("ENCODER:\n", e)
    print("MODEL:\n", e.model)
    print("PROGRAM:\n", e.program)
    phi = And(e.model, e.program, e.sensitiveAttribute, e.fairnessTarget)
    pvars = [x for x in e.vdist]
    print("VDIST:\n", e.vdist)
    print("PHI:\n", phi)
    print("NONPROBVARS PROJECTED:\n", projectNonProbVars(phi,pvars,False)) 

# print codegen.to_source(node)

