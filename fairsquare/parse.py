import ast
import codegen
from z3 import *
from collections import namedtuple


def parse_fr(filename):
    f = open(filename, "r")
    node = ast.parse(f.read())
    f.close()

    e = Z3Encoder()
    e.visit(SATransformer().visit(node))

    ParsedFR = namedtuple('ParsedFR', ['vdist', 'model', 'program', 'sensitiveAttribute', 'fairnessTarget', 'qualified'])

    return ParsedFR(e.vdist, e.model, e.program, e.sensitiveAttribute, e.fairnessTarget, e.qualified)


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
    if isAdd(op): return l + r
    elif isSub(op): return l - r
    elif isMult(op): return l * r
    elif isDiv(op): return l / r
    elif isPow(op): return l**r
    else: assert False, "Weird binary op"
def makeUnary(op, e):
    if isUSub(op): return -e
    elif isUAdd(op): return +e
    elif isNot(op): return Not(e)
    else: assert False, "Weird unary op"
def makeBool(op, *args):
    if isAnd(op): return And(*args)
    elif isOr(op): return Or(*args)
    else: assert False, "Weird bool op"
def makeCompare(op, l, r):
    if isGt(op): return l > r
    elif isGtE(op): return l >= r
    elif isLt(op): return l < r
    elif isLtE(op): return l <= r
    elif isEq(op): return l == r
    elif isNotEq(op): return l != r
    else: assert False, "Weird compare"


# A visitor class to gather all identifiers (used for SSA)
class VDict(ast.NodeVisitor):
    def __init__(self):
        self.vdict = {}
    def visit_Name(self, node):
        # print(ast.dump(node))
        self.vdict[node.id] = 0
def init_dict(node):
    v = VDict()
    v.visit(node)
    return v.vdict


# Convert the AST to (not-static) single-assignment form
# This facilitates a direct conversion to logical formulae
class SATransformer(ast.NodeTransformer):

    # This is specialized for our .fr files where we assume
    # all variables in all functions use the same namespace
    def visit_Module(self, node):
        d = VDict()
        d.visit(node)
        self.live_index = d.vdict
        
        for func in node.body:
            func.body = self.doSeq(func.body)

        return node

    def doSeq(self, seq):
        return [self.visit(s) for s in seq]

    def visit_If(self, node):

        self.visit(node.test)

        before = self.live_index.copy()
        node.body = self.doSeq(node.body)
        after_if = self.live_index.copy()

        self.live_index = before.copy()
        node.orelse = self.doSeq(node.orelse)
        after_else = self.live_index.copy()
        
        for name in self.live_index:
            # Updated indexes after the join
            self.live_index[name] = max(after_if[name], after_else[name])
            # Add the Phi functions
            if after_if[name] < self.live_index[name]:
                node.body.append(ast.Assign(targets=[ast.Name(id=self.append_index(name))],    
                                            value=ast.Name(id=name+"_"+str(after_if[name]))))
            if after_else[name] < self.live_index[name]:
                node.orelse.append(ast.Assign(targets=[ast.Name(id=self.append_index(name))],    
                                              value=ast.Name(id=name+"_"+str(after_else[name]))))

        return node

    def visit_Assign(self, node):
        assert len(node.targets) == 1
        value = self.visit(node.value)
        self.live_index[node.targets[0].id] += 1
        target = self.visit(node.targets[0])
        return ast.Assign(targets=[target], value=value)

    def visit_Call(self, node):
        for arg in node.args:
            self.visit(arg)
        return node

    def visit_Name(self, node):
        #return ast.copy_location(
        #        ast.Name(id=self.append_index(node.id), ctx=ast.Load()),
        #        node)
        return ast.Name(id=self.append_index(node.id))

    def append_index(self, name):
        if name not in self.live_index:
            self.live_index[name] = 0
        return name + "_" + str(self.live_index[name])


# Essentially syntactic conversion of an SA-form Python program
# to a z3 formula
class Z3Encoder(ast.NodeVisitor):

    def __init__(self):
        self.vdist = {}
        
        self.sensitiveAttribute = None
        self.qualified = None
        self.fairnessTarget = None

        self.model = None
        self.program = None

    def generic_visit(self, node):
        assert False, "generic visit on node " + ast.dump(node)

    # ast.parse returns an ast.Module object
    def visit_Module(self, node):
        # The body of node should be a sequence:
        # popModel, F
        assert len(node.body) == 2
        fpop = node.body[0]
        fprog = node.body[1]
        assert fpop.name == 'popModel' and fprog.name == 'F'

        self.model = self.doSeq(fpop.body)
        self.model = simplify(self.model)

        self.program = self.doSeq(fprog.body)
        self.program = simplify(self.program)

    def doSeq(self, seq):
        trans = [self.visit(s) for s in seq]
        return simplify(And(*trans))

    def visit_Assign(self, node):
        assert len(node.targets) == 1, "No unpacked tuple assignments allowed"

        lhs = node.targets[0].id
        rhs = node.value

        if isCall(rhs):
            # TODO Ensure that the call is to a supported prob dist
            return self.probAssign(lhs, rhs)

        zrhs = exprToZ3(rhs)
        zlhs = Real(lhs)
        return zlhs == zrhs

    # Store the probability distribution information
    def probAssign(self, name, rhs):
        fname = rhs.func.id
        
        if fname == 'gaussian':
            mean = evalAST(rhs.args[0])
            variance = evalAST(rhs.args[1])
            assert variance > 0
            self.vdist[Real(name)] = ('G', mean, variance)
            return True

        elif fname == 'step':
            step = ('S', evalAST(rhs.args[0]))
            self.vdist[name] = step

            # sum of probs == 1
            l = evalAST(rhs.args[0])
            s = sum([a[2] for a in l])
            assert abs(s-1.0) <= 0.00001

            lbounds = [a_b_c[0] for a_b_c in step[1]]
            ubounds = [a_b_c1[1] for a_b_c1 in step[1]]
            res = And(Real(name) >= min(lbounds), Real(name) <= max(ubounds))

            return res

        else:
            assert False, "Supported distributions: gaussian, step"

    def visit_If(self, node):
        zcond = exprToZ3(node.test)
        zthen = self.doSeq(node.body)
        zelse = self.doSeq(node.orelse)
        return And(Implies(zcond, zthen), Implies(Not(zcond), zelse))

    # We handle expressions in If and Assign statements in their visit_'s,
    # so these are always calls
    def visit_Expr(self, node):
        assert isinstance(node.value, ast.Call), "Unexpected expression"
        return self.visit(node.value)

    # We handle calls to prob vars in visit_Assign,
    # so these are always the program annotations
    def visit_Call(self, node):
        fn = node.func.id
        assert len(node.args) == 1

        if fn == 'sensitiveAttribute':
            self.sensitiveAttribute = exprToZ3(node.args[0])
        elif fn == 'qualified':
            self.qualified = exprToZ3(node.args[0])
        elif fn == 'fairnessTarget':
            self.fairnessTarget = exprToZ3(node.args[0])
        else:
            assert False, "Unrecognizable function call"

        return True

    def visit_Return(self, node):
        return True


# e is an ast expr node
def exprToZ3(e):
    if isBinOp(e):
        zlhs = exprToZ3(e.left)
        zrhs = exprToZ3(e.right)
        return makeBin(e.op, zlhs, zrhs)
    elif isUnaryOp(e):
        zexp = exprToZ3(e.operand)
        return makeUnary(e.op, zexp)
    elif isBoolOp(e):
        zexprs = [exprToZ3(v) for v in e.values]
        return makeBool(e.op, zexprs)
    elif isCompareOp(e):
        assert len(e.ops) == 1
        zlhs = exprToZ3(e.left)
        zrhs = exprToZ3(e.comparators[0])
        return makeCompare(e.ops[0], zlhs, zrhs)
    elif isNum(e):
        return evalAST(e)
    elif isName(e):
        if isTrue(e):
            return BoolVal(True)
        elif isFalse(e):
            return BoolVal(False)
        else:
            return Real(e.id)
    else:
        assert False, "Weird expression" + ast.dump(e)



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

