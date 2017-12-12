# "numericalPhi.py"
# 2017-02-07
# David Merrell
#
# A class that bridges the gap between Z3 and
# numpy. It encodes a formula in linear real arithmetic
# as a set of linear constraints
# Ax <= b.
# For now it just handles atoms containing linear
# real inequalities.
#
# Can be used to query whether an array of points
# satisfy the formula.
#
# Assumes variables in the Z3 formula are represented
# by Z3_OP_UNINTERPRETED 


from z3 import *
import numpy as np
import z3extra as ze

class numericalPhi:

    def __init__(self, phi, valid_vars=None):
        """
        Build an instance from a Z3 formula.
        Ignore predicates containing variables
        that are not valid.
        """

        print("PHI:\n", phi)
        phi = simplify(phi)
        print("SIMPLIFIED PHI:\n", phi)
        phi = ze.close(phi)
        print("CLOSED PHI:\n", phi)

        if valid_vars is None:
            self.valid_vars = list(ze.get_vars(phi))
        else:
            self.valid_vars = valid_vars

        # We want to encode the formula as a set
        # of linear constraints of the form 
        # Ax <= b.
        # We also want a vector d which maps
        # inequalities to disjuncts.
        self.constraint_matrix, self.constraint_b, self.disjuncts \
                = self.build_constraint_matrix(phi)


    def build_constraint_matrix(self, phi):
        """
        Build arrays A, b, and d, where
        Ax <= b
        encodes the formula, and d keeps track
        of which disjunct a given constraint
        belongs to.
        """
        A = []
        b = []
        d = []

        disjuncts = ze.toDNF(phi)
        print("DNF PHI:\n", disjuncts) 
       
        dj_num = 0
        for dj in disjuncts:
            print("DISJUNCT:\n", dj)

            if is_and(dj):
                for cnj in dj.children():
                  
                    print("CONJUNCT:\n", cnj)
                    row, bconst = self.extract_ineq(cnj)
                    if len(row) > 0 :
                        A.append(row)
                        b.append(bconst)
                        d.append(dj_num)

            elif is_le(dj) or is_ge(dj) or is_not(dj):
                row, bconst = self.extract_ineq(dj)
                if len(row) > 0 :
                    A.append(row)
                    b.append(bconst)
                    d.append(dj_num)

            dj_num += 1


        A = np.array(A)
        b = np.array(b)
        d = np.array(d)

        assert A.shape[0] == len(b) and len(d) == len(b)

        return A, b, d


    def extract_ineq(self, ineq):

        row = []
        vnames = self.get_ordered_varnames()

        # We only handle linear inequalities.
        if not( is_le(ineq) or is_ge(ineq) or is_not(ineq) ):
            print("{} is not an inequality.".format(ineq))
            return [], 0  

        # Handle the case of a negation.
        sgn = 1
        if is_not(ineq):
            assert len(ineq.children()) == 1,\
                    "{} MULTIPLE CHILDREN?".format(ineq)
            ineq = ineq.children()[0]
            sgn = -1

        # Extract coefficients from the atom.
        coeff_d = self.get_atom_coeffs(ineq)

        # Add coefficients to the row vector, in order.
        for vname in vnames:
            if vname in list(coeff_d.keys()):
                row.append(coeff_d[vname])
            else:
                row.append(0.0) 

        if len(coeff_d) > 0 and 'const' in list(coeff_d.keys()):
            return row, coeff_d['const']
        else:
            return [], 0



    def get_ordered_vars(self):
        """
        Get the variables in sorted (alphabetical) order.
        This is the canonical ordering used throughout
        this class.
        """
        varnames = np.array([str(v) for v in self.valid_vars])
        sorted_inds = np.argsort(varnames)
        sorted_vars = list(varnames[sorted_inds])

        return sorted_vars


    def get_ordered_varnames(self):
        """
        Get variable names in sorted order.
        """
        return [str(v) for v in self.get_ordered_vars()]


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
            print("ATOM MUST BE A <= or >=")
            return {}
        
        # Return a dictionary whose keys are gaussian variable names
        # and values are numbers --- the variables' coefficients.
        coeff_d = {}

        varnames = self.get_ordered_varnames() 

        lhs_coeffs = self.get_lin_coeffs(atom.children()[0])
        rhs_coeffs = self.get_lin_coeffs(atom.children()[1])
        key_union = set(lhs_coeffs.keys()).union(list(rhs_coeffs.keys()))

        for varname in key_union:
            coeff_d[varname] = 0
            if varname in varnames:
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
                print("{}\nNOT A USABLE ATOM.".format(atom))
                return {}

        return coeff_d 
      

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


if __name__=="__main__":

    x = Real('x')
    y = Real('y')

    phi = Or(And(2*x + 3*y < 100,x > 0, y > 0), And(-7*x + -4*y < 123, x < -2, y < -3)) 

    nphi = numericalPhi(phi)

    print("MATRIX:\n", nphi.constraint_matrix)
    print("RHS:\n", nphi.constraint_b)
    print("DISJUNCTS:\n", nphi.disjuncts)



