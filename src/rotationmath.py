# "rotationmath.py"
# David Merrell
# 2016-11-17
#
# Math related to rotations.
#

import numpy as np
import pandas as pd
#from z3 import *
import fractions as frc
import pickle as cpkl
from sklearn.neighbors import KDTree
from scipy.stats import norm
from functools import reduce

def qr_align(targetvec):
    """
    Given a 'target' vector, find the basis
    vector that projects onto it maximally.
    Then replace that basis vector with the 
    target vector, and obtain an orthonormal basis
    that is aligned with the target vector.

    We employ numpy's QR factorization method,
    observing that it generates an orthonormal 
    basis whose 1st vector is parallel to the 
    original basis' 1st vector.

    When the target vector is multiplied by the
    output matrix, the resulting vector has
    a single non-zero entry corresponding to the
    maximal basis vector described above.
    """

    d = len(targetvec)

    targetvec = np.array(targetvec)
    assert len(targetvec.shape) == 1, "targetvec should be 1-dimensional"

    original_basis = np.identity(d)
    #print "ORIGINAL BASIS\n", original_basis

    projections = np.absolute(np.dot(original_basis, targetvec))
    #print "PROJECTIONS\n", projections

    # Get the original basis vector 
    # which is the best match
    best_match_ind = np.argmax(projections)
    #print "BEST MATCH:\n", best_match_ind

    # Get the indices that weren't the best match
    not_best_inds = [i for i in range(d) if i != best_match_ind]

    # Compose a new matrix 
    newmat = np.concatenate((targetvec.reshape(d,1),
                             original_basis[:,not_best_inds]),axis=1)
    #print "NEWMAT:\n", newmat

    # QR factorize the new matrix. The resulting Q should
    # contain an orthonormal basis, with the first vector
    # parallel to the target vector.
    q, r = np.linalg.qr(newmat)

    #print "ORTHONORMAL BASIS\n", q

    #temp = np.copy(q[:,best_match_ind])
    #q[:,best_match_ind] = q[:,0]
    #q[:,0] = temp[:]

    #return np.transpose(q)
    return q



def b_generator(n):
    """
    Generate (n-1)-tuples of integers that enable 
    the computation of pythagorean n-tuples.
    """

    assert n >= 3

    ls = [0 for i in range(n-3)]
    ls = ls + [1,2]

    yield ls

    while True:
        for ind, num in enumerate(ls):
            if ind < n-3:
                if num < ls[ind+1]:
                    ls[ind] = num+1
                    ls[:ind] = [0 for i in range(ind)]
                    yield ls
                    break
                    
            elif ind == n-3:
                if num < ls[ind+1]-1:
                    ls[ind] = num+1
                    ls[:ind] = [0 for i in range(ind)]
                    yield ls
                    break

            elif ind > n-3:
                ls[ind] = num+1
                ls[:ind] = [0 for i in range(ind)]
                yield ls




def pythagorean_ntuples(n, lim):
    """
    Generate pythagorean n-tuples
    """

    assert(n > 2)

    gen = b_generator(n)
    tpl = next(gen)
    res_set = set()

    count = 0
    while True:

        result = [0 for i in range(n)]

        result[-2] = tpl[-1]*tpl[-1]
        for x in tpl[:-1]:
            result[-2] = result[-2] - x*x

        for ind, x in enumerate(tpl[:-1]):
            result[ind] = 2*tpl[-1]*x

        result[-1] = 0
        for ind, x in enumerate(tpl):
            result[-1] = result[-1] + x*x

        lhs = 0
        for a in result[:-1]:
            lhs += a*a
        rhs = result[-1]*result[-1]
        assert lhs == rhs

        tpl = next(gen)

        div = reduce(frc.gcd, result)
        result = [abs(res) / div for res in result]

        count += 1
        #print count, tuple(sorted(result))

        if max(result) >= lim:
            break

        res_set.add(tuple(sorted(result)))
        if len(res_set) % 10000 == 0:
            print(len(res_set))

    return res_set


def precompute_rational_vecs(dim,digit_lim=4,arr_filename="pythagorean/tuples"):
                             #tree_filename="kdtree"):
    """
    Generate (dim)-dimensional unit vectors
    whose entries are rational. Store them
    for later use. 
    """

    # Generate the pythagorean (dim+1)-tuples
    pyth_tuples = list(pythagorean_ntuples(dim+1,10**digit_lim))

    # Store the raw pythagorean n-tuples
    pyth_df = pd.DataFrame(pyth_tuples)
    pyth_df.to_csv(arr_filename+".csv")

    # Make a KD Tree object for each d <= dim.
    #rat_vecs = [[1.0*num/tpl[-1] for num in tpl[:-1]] for tpl in pyth_tuples]

    #for d in range(2,dim):
    #tree = KDTree(rat_vecs)
    #cpkl.dump(tree,open(tree_filename+"_{}_{}.csv".format(dim,digit_lim),"wb"))

    return


def nearest_rational_vec(in_vec,tuple_file="pythagorean/triples.csv"):
    """
    Given a vector of dimension d, return
    a unit vector with rational entries
    that is most aligned with it.
    """

    in_dim = len(in_vec)
   
    pyth_df = pd.read_csv(tuple_file,index_col="Unnamed: 0")
    stored_dim = len(pyth_df.columns) - 1

    # If the precomputed n-tuples do not cover 
    # enough dimensions, then we precompute a larger set
    # and try again.
    if stored_dim < in_dim:
        precompute_rational_vecs(in_dim,3.5)
        return nearest_rational_vec(in_vec)

    # Get the pythagorean (in_dim+1)-tuples.
    pyth_arr = np.array(pyth_df)
    sums = np.sum(pyth_arr[:,list(range(stored_dim-in_dim))],axis=1)
    relevant_arr = pyth_arr[sums == 0]
    comps = relevant_arr[:,list(range(stored_dim-in_dim,stored_dim))]
    norms = relevant_arr[:,stored_dim]
    vecs = np.apply_along_axis(lambda x: x * 1.0/ norms, arr=comps, axis=0)

    # positivize and sort in_vec for matching purposes.
    neg_inds = np.sign(in_vec)
    in_pos = in_vec * neg_inds
    sorted_inds = np.argsort(in_pos)
    in_pos_srt = in_pos[sorted_inds]

    # Find the best match
    best_ind = np.argmax(np.dot(vecs, in_pos_srt))
    best_pos_srt = relevant_arr[best_ind,list(range(stored_dim-in_dim,stored_dim))]

    # Un-sort and un-positivize the match
    best_pos = best_pos_srt[np.argsort(sorted_inds)]
    best = best_pos * neg_inds

    best = [int(b) for b in best]
    best = np.append(best,int(norms[best_ind]))

    return best


def build_householder(vec, k, vecnorm=None):
    """
    Construct a Householder reflection
    matrix that reflects the vector
    "vec" into alignment with the 
    kth standard basis vector.
    """

    vec = np.array(vec,dtype=int)
    I = np.identity(len(vec), dtype=int)

    if vecnorm is None:
        vecnorm = int(np.linalg.norm(vec))

    diff = vec - I[k,:]*vecnorm

    M = np.outer(diff,diff)
    denom = int(np.inner(diff,diff))

    if denom != 0:
        M = [[frc.Fraction(int(num),denom) for num in row] for row in M]
        M = np.array(M)

    return I - 2*M



def rh_align(target_vec,precision):
    """
    Given a "target" vector, we return
    a reflection matrix with rational entries,
    that approximately reflects that target vector
    into alignment with one of the standard basis vectors.
    """

    #match = nearest_rational_vec(target_vec,"pythagorean/tuples.csv")
    match = nearest_rational_vec2(target_vec,dig_lim=precision)
    h_mat = build_householder(match[:-1],
                              np.argmax(np.abs(match[:-1])),
                              vecnorm=match[-1])

    return h_mat



def mc_gauss_surface_integrals(constraint_A,constraint_b, dsjncts, mc_sample=1000):
    """
    Given a formula (i.e. region in d-space),
    approximate the "mass" associated with
    each face of the formula; "mass" refers to 
    the joint probability density integrated
    over a face. Used for selecting
    rotations.
    """

    #print constraint_A

    #print constraint_b

    # Suppose we get coefficient vectors for each face
    # and store them in k x d matrix A. We also get the
    # "offset" associated with each of these k faces.
    k = constraint_A.shape[0]
    assert len(constraint_b.shape) == 1
    assert k == constraint_b.shape[0]

    constraint_b = np.reshape(constraint_b,(k,1))

    masses = []

    # Each row corresponds to a face of 
    # a region; each gets a mass computed for it
    for row in range(k):
    
        dj = dsjncts[row]

        # We compute an affine transformation from
        # (d-1)-space to d-space (a parameterization).
        c_vec = constraint_A[row,:]
        d = len(c_vec)
        affine_A = qr_align(c_vec)
        affine_A = affine_A[:,1:]
        affine_b = np.reshape(c_vec * constraint_b[row] / np.dot(c_vec,c_vec),
                              (d,1))

        samples = np.random.normal(size=[d-1,mc_sample])
        transformed = np.dot(affine_A,samples) + affine_b

        sat_inds = np.logical_and(dsjncts == dj, np.arange(k) != row)
        satmat = constraint_A[sat_inds,:]
        satisfies = np.dot(constraint_A[sat_inds,:],transformed)\
                        <= constraint_b[sat_inds,:]

        satsums = np.sum(satisfies,axis=0)

        # We multiply the face's "conditional probability
        # mass" by this in order to get the face's "probability".
        #print "NORM OF C_VEC:", np.linalg.norm(c_vec)
        #print "SIZE OF CONSTRAINT b ENTRY:", constraint_b[row]
        factor = float(norm.pdf(constraint_b[row] / np.linalg.norm(c_vec)))
        masses.append(factor * len(satsums[satsums == satmat.shape[0]]) / mc_sample)

    return np.array(masses)


################################################################
def pairwise_rh_align(target_vec,flex_order=True,precision=5):
    """
    Given a "target" vector, we return
    a reflection matrix with rational entries
    that approximately reflects that target vector
    into alignment with one of the standard basis vectors.

    Differs from "rh_align"; we use the insight 
    that this can be done with a composition
    of householder reflections, each acting on 
    a pair of dimensions. So we only ever need to 
    use pythagorean triples (rather than n-tuples).
    """

    if np.linalg.norm(target_vec) == 0:
        return np.identity(len(target_vec), dtype=int)

    # Make sure target is normalized.
    target_vec = np.array(target_vec) / np.linalg.norm(target_vec)

    # Figure out the order in which we'll visit dimensions.
    # For now, we go in the order of decreasing component size
    # in the target vector.
    if flex_order:
        dim_order = np.argsort(np.abs(target_vec))[::-1][:len(target_vec)]
    else:
        dim_order = np.array(list(range(len(target_vec))))

    result = np.identity(len(target_vec),dtype=int)
    
    # Handle the case where the target vector is a standard basis vector.
    if abs(abs(np.sum(target_vec))-1)<1e-9 and abs(np.dot(target_vec,target_vec)-1)<1e-9:
        target_ind = np.argmax(np.abs(target_vec))
        swp = np.copy(result[target_ind,:])
        result[target_ind,:] = np.copy(result[dim_order[0],:])
        result[dim_order[0],:] = swp[:]
        return result

    # We proceed by starting with a standard basis 
    # vector and applying reflections to it until
    # it becomes the target vector.
    curr_vec = np.identity(len(target_vec))[:,dim_order[0]]
    for i, dim in enumerate(dim_order[:-1]):

        # Destination vector (numerical)
        comp_sgn = np.sign(target_vec[dim_order[i+1]])
        dest_vec = np.zeros((len(curr_vec),))
        dest_vec[dim_order[:(i+1)]] = target_vec[dim_order[:(i+1)]]
        subvec_size = min(1.0,np.dot(target_vec[dim_order[:(i+1)]],
                                     target_vec[dim_order[:(i+1)]]))

        dest_vec[dim_order[i+1]] = comp_sgn*(1.0 - subvec_size)**0.5
        #print "DEST VEC!!!", dest_vec
        #print "DEST VEC NORM!!!", np.linalg.norm(dest_vec)

        # Difference vector (numerical)
        diff = dest_vec - curr_vec
        if np.max(np.abs(diff)) == 0:
            return np.identity(len(curr_vec),dtype=int)
        diffnormed = diff / np.linalg.norm(diff)
        diff2 = [diffnormed[dim_order[i]],diffnormed[dim_order[i+1]]]

        # Difference vector (rational approximation)
        rdiff = nearest_rational_vec2(diff2,dig_lim=precision)
        rdiff_norm = rdiff[-1]
        rdiff = [frc.Fraction(num,rdiff_norm) for num in rdiff[:-1]]

        # Build a householder matrix (from rational approximation)
        hh2 = np.identity(2,dtype=int) - 2*np.outer(rdiff,rdiff)
        hh = np.identity(len(curr_vec),dtype=object)
        hh[dim,dim] = hh2[0,0]
        hh[dim_order[i+1],dim] = hh2[1,0]
        hh[dim,dim_order[i+1]] = hh2[0,1]
        hh[dim_order[i+1],dim_order[i+1]] = hh2[1,1]

        # Update the result matrix, update "current" vector
        result = np.dot(hh,result)
        curr_vec = dest_vec[:]

    return result


################################################################
def full_rational_align(target_vecs,precision=5,maxdig=3):
    """
    Assume the target vectors are given in order
    of decreasing priority; that is, we wish to be
    aligned most with the first row, next-most aligned
    with the second, etc. Then this function
    returns an orthogonal matrix with rational entries,
    constructed from composed householder reflections.
    The result is constrained in m dimensions, where

    m = min(number of target vectors, number of dimensions)
    """

    target_vecs = np.array(target_vecs)[:]
    #print "TARGET ALIGNMENT VECTORS:\n", target_vecs
    
    # Eliminate redundancies
    normed_vecs = target_vecs / np.reshape(np.linalg.norm(target_vecs,axis=1),
                                           (target_vecs.shape[0],1))
    absprods = np.abs(np.dot(normed_vecs,normed_vecs.T))
    is_one = np.abs(absprods - 1) < 1e-9
    keep_inds = []
    for row in range(is_one.shape[0]):
        not_copy = True
        for col in range(row):
            if is_one[row,col]:
                not_copy = False
                break
        if not_copy:
            keep_inds.append(row)
    target_vecs = target_vecs[keep_inds,:]
    print("NONREDUNDANT TARGET ALIGNMENT VECTORS:\n",target_vecs)
    d = target_vecs.shape[1]

    result = np.identity(d, dtype=object)
    num_res = np.identity(d,dtype=float)

    for i in range(d-1):

        # If we run out of target vecs, we're done.
        if i >= target_vecs.shape[0]:
            break
        # If we accumulate too many digits, we're done.
        few_digs = True
        for row in range(result.shape[0]):
            for col in range(result.shape[1]):
                if result[row,col].numerator > 10**maxdig\
                        or result[row,col].denominator > 10**maxdig:
                    few_digs = False
                    break
        if not few_digs:
            return result

        target = target_vecs[i,:]

        orthog_space_basis = 1.0 * result[:,i:] 
        projected_target = np.dot(target,orthog_space_basis)

        rot = np.identity(d,dtype=object)
        rot[i:,i:] = pairwise_rh_align(1.0*projected_target,flex_order=False,precision=precision)

        result = np.dot(result,rot)

    return result


################################################################
def sigfig_str(num,figs):
    lead_place = np.floor(np.log10(abs(num)))
    trunc_place = int(np.around(figs - lead_place -1))
    num = np.around(num,trunc_place)
    str_trunc_place = abs(max(trunc_place, 0))
    format_bed = "{:." + str(str_trunc_place) + "f}"
    return format_bed.format(num)
   

################################################################
def nearest_rational_vec2(target_vec, dig_lim=4):
    """
    Exploit the methods described by Shiu
    (http://www.jstor.org/stable/3617358?seq=1#page_scan_tab_contents)
    to generate a rational-valued 2D unit vector near
    to the target vector, without exceeding a certain
    number of digits.
    """

    # Put the target vector in the first quadrant,
    # and sort its indices into decreasing order
    # (want an angle in [0,pi/4]
    target_sgns = np.sign(target_vec)
    target_pos = target_vec*target_sgns
    target_inds = np.argsort(target_pos)[::-1][:2]
    target_pos_srt = target_pos[target_inds]

    # Obtain u
    angle = np.arctan2(target_pos_srt[1],target_pos_srt[0])
    if abs(angle - np.pi/2) < 1e-9:
        approx_pos_srt = np.array([1,0])
        approx_pos = approx_pos_srt[np.argsort(target_inds)]
        approx = approx_pos * target_sgns
        return [int(approx[0]),int(approx[1]),1] 

    u = np.tan(angle) + (1.0/np.cos(angle))

    # Iteratively converge on u
    rational = frc.Fraction(1,1)

    def cf_(f,i_list,prev_r):

        if abs(f - 0) < 1e-10:
            return contdfrac_to_frac(i_list)

        else:

            i_next = int(np.floor(1.0/f))
            f_next = (1.0/f) - i_next
            i_list.append(i_next)
            rational = contdfrac_to_frac(i_list)

            x = 2*rational.numerator*rational.denominator
            y = rational.numerator**2 - rational.denominator**2
            norm = rational.numerator**2 + rational.denominator**2
            div = reduce(frc.gcd, [x,y,norm])
            
            if not(x/div < 10**dig_lim and y/div < 10**dig_lim and norm/div < 10**dig_lim):
                return prev_r
            else:
                return cf_(f_next,i_list,rational)
        
    i0 = int(np.floor(u))
    f0 = u - i0
    rational = cf_(f0,[i0],rational)    
           
    approx_pos_srt = np.abs(np.array([2*rational.numerator*rational.denominator,
                               rational.numerator**2 - rational.denominator**2]))

    #if np.dot([1,0],target_pos_srt) == 1:
    #    approx_pos_srt = np.array([1,0])

    approx_pos = approx_pos_srt[np.argsort(target_inds)]
    approx = approx_pos * target_sgns

    norm = rational.numerator**2 + rational.denominator**2
    div = reduce(frc.gcd, [approx[0],approx[1],norm])
    approx = approx / div
    norm = norm / div
    return [int(approx[0]),int(approx[1]),int(norm)] 


################################################################
def float_to_contdfrac(x, k):
    """
    Given a positive number x, return a truncation
    of its continuing fractional form to a depth
    of k. For example, 
    (x, k) -> a0 + 1/(a1 + 1/(a2 + 1/(a3 + ... )...)
    We return it as a list of these a's: [a0,a1,...,ak].
    """
    i0 = int(np.floor(x))
    f0 = x - i0

    def cf_(f,i_list,k):
        if k <= 0 or abs(f - 0) < 1e-9:
            return i_list
        else:
            i_next = int(np.floor(1.0/f))
            f_next = (1.0/f) - i_next
            i_list.append(i_next)
            return  cf_(f_next,i_list,k-1)

    return cf_(f0,[i0],k) 

################################################################
def contdfrac_to_float(alist):

    assert len(alist) >= 1
    if len(alist) == 1:
        return alist[0]
    else:
        return alist[0] + 1.0/(contdfrac_to_float(alist[1:]))

################################################################
def contdfrac_to_frac(alist):

    assert len(alist) >= 1
    if len(alist) == 1:
        return alist[0]
    else:
        return alist[0] + frc.Fraction(1, (contdfrac_to_frac(alist[1:])))


if __name__=='__main__':
   

    # Test tuple generator
    #precompute_rational_vecs(8,digit_lim=3.5)
    #precompute_rational_vecs(2,digit_lim=5,arr_filename="pythagorean/triples_5")

    #query = [10,-3,4,-40]
    #print query
    #match = nearest_rational_vec(query)
    #print match

    #mat = build_householder(match[:-1],0)
    #print mat
    #print np.dot(mat,mat)

    # Tabulate some tuple counts, as an example.
    #import pandas as pd
    #ns = [3,4,5]
    #lims = [100,1000,10000]
    #df = pd.DataFrame(index=ns, columns=lims)
    #for n in ns:
    #    for lim in lims:
    #        df.loc[n,lim] = len(pythagorean_ntuples(n,lim))
    #print df


    # Selection by KD Tree
    # vs 
    # Selection by Inner Product

    # Test Pairwise Reflection alignment
    #target1 = np.array([1,1,1])/3**0.5
    #target2 = -1*np.array([1,1,1])/3**0.5
    #mat1 = pairwise_rh_align(target1,flex_order=True)
    #mat2 = pairwise_rh_align(target2,flex_order=True)
    #print "MAT1:\n", 1.0*mat1
    #print "MAT2:\n", 1.0*mat2
    #print "THIS SHOULD BE IDENTITY:\n", 1.0*np.dot(mat1,np.transpose(mat1))
    #print "PROD 1:\n", np.dot(target1/np.linalg.norm(target1),mat1)
    #print "PROD 2:\n", np.dot(target2/np.linalg.norm(target2),mat2)
    #print "GOOD MATRIX AND BAD VECTOR:\n", np.dot(target2/np.linalg.norm(target1),mat1)

    #targets = np.array([[1,1,1],[2,-1,-1],[0,1,-1],[9,10,11]])
    #targets = targets / np.reshape(np.linalg.norm(targets,axis=1),(4,1))
    #print "TARGETS:\n",targets

    #mat = full_rational_align(targets)

    #print "MAT:\n", mat
    #print "TEST 1st TARGET:", np.dot(mat,np.transpose(targets))

    
    target1 = [np.array([-1,-1,-1])/3**0.5]
    target2 = [np.array([1,1,1])/3**0.5]
    print("TARGET 1:", target1)
    print("TARGET 2:", target2)
   
    mat1 = full_rational_align(target1)
    mat2 = full_rational_align(target2)
    print("MAT 1:\n", 1.0*mat1)
    print(1.0*np.dot(mat1,mat1.T))
    print("MAT 2:\n", 1.0*mat2)
    print(1.0*np.dot(mat2,mat2.T))

    #print "PROD 1:", np.dot(target1, mat1)
    #print "PROD 2:", np.dot(target2, mat2)

    #print "TRY GOOD MATRIX WITH BAD VECTOR:\n", np.dot(target1, mat2)
    #print "SO THE MATRIX IS BAD."
    #num1 = np.array(1.0*mat1, dtype=float)
    #num2 = np.array(1.0*mat2, dtype=float)
    #print "DET 1:", np.linalg.det(1.0*num1)
    #print "DET 2:", np.linalg.det(1.0*num2)

    #print "SHOULD BE 12300:", sigfig_str(12345.123, 3)
    #print "SHOULD BE .12:", sigfig_str(.12345, 2)
    #print "SHOULD BE -0.01:", sigfig_str(-0.009876, 1)
    #print "SHOULD BE -987.2:", sigfig_str(-987.193, 4)
    
    #phi = (1.0+np.sqrt(5))/2.0
    #gr = float_to_contdfrac(phi, 20)
    #print "GOLDEN RATIO SHOULD BE ALL ONES:", gr
    #print "numerical diff:", contdfrac_to_float(gr) - phi
    #print "rational approx:", contdfrac_to_frac(gr)

    #gr = float_to_contdfrac(1.75, 20)
    #print "SHOULD END EARLY:", gr
    #print "FRACTION:", contdfrac_to_frac(gr)
    #print "numerical:", contdfrac_to_float(gr)

    #print "NEAREST RATIONAL UNIT VECTOR:", nearest_rational_vec2([-1,-1])
    #print "NEAREST RATIONAL UNIT VECTOR:", nearest_rational_vec2([7,24])
    #print "NEAREST RATIONAL UNIT VECTOR:", nearest_rational_vec2([-2,-25],dig_lim=3)
