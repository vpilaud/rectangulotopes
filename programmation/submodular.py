# submodular functions for rectangulotopes

from itertools import *

## SUBMODULAR FUNCTIONS

# simplex for an up arc I
# 1 if S intersects I, 0 otherwise
def f_loday(S, I):
    return int(bool(set(S) & set (I)))

# anti-simplex for a down arc I
# -1 if I is a subset of S, 0 otherwise
def f_antiloday(S, I):
    return - int(set(I).issubset(set(S)))

# simplex + antisimplex (cuboctahedron in 3d)
# 1 if S intersects I but I is not a subset of S, 0 otherwise
def f_cuboc(S, I):
    return int(bool(set(S) & set (I)) and not set(I).issubset(set(S)))

# weak rectangulotopes
# number of I such that S intersects I but I is not a subset of S
def f_weakrect(S, n):
    return sum([int(bool(set(S) & set (range(i, j+1))) and not set(range(i, j+1)).issubset(set(S))) for i in range(n-1) for j in range(i+1, n)])

# yin shard polytopes
# 1 if and only if S intersects I, and J is not a subset of S
def f_yin(S, I, J):
    return int(bool(set(S) & set(I)) and not (set(J).issubset(S)))

# yin rectangulotopes
# number of pairs I, J such that S intersects I, and J is not a subset of S
def f_yinrect(S, n):
    return sum([int (bool(set(S) & set(range(i, j+1))) and not (set(range(j+1, k+1)).issubset(S))) for i in range(n-1) for j in range(i, n-1) for k in range(j+1, n)])

# yang shard polytopes
# 1 if and only if S intersects J, and I is not a subset of S
def f_yang(S, I, J):
    return int(bool(set(S) & set(J)) and not (set(I).issubset(S)))

# yang rectangulotopes
# number of pairs I, J such that S intersects J, and I is not a subset of S
def f_yangrect(S, n):
    return sum([int (bool(set(S) & set(range(j+1, k+1))) and not (set(range(i, j+1)).issubset(S))) for i in range(n-1) for j in range(i, n-1) for k in range(j+1, n)])

# strong rectangulotopes
def f_strongrect(S, n):
    return f_yinrect(S, n) + f_yangrect(S, n)

## GENERIC FUNCTIONS

# powerset recipe
def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

# compute the polyhedron for a given submod function
# f takes two arguments : subset S and size of ground set n
def submod_polytope(f, n):
    all_ineq = [ (f(S, n),) + tuple(-1 if x in S else 0 for x in range(n)) for S in powerset(range(n)) ]
    return Polyhedron(ieqs = all_ineq, eqns = [all_ineq[-1]]) # the last ineq is an eq

## RECTANGULOTOPES

# associahedron as a Minkowski sum
def assoc(n):
    return sum([submod_polytope(lambda S, m : f_loday(S, range(i, j+1)), n) for i in range(n-1) for j in range(i+1, n)])

# weak rectangulation polytope (diagonal / Baxter)
# either sum all lodays and antilodays, or the cubocs
def wrp(n):
    return sum([submod_polytope(lambda S, m : f_cuboc(S, range(i, j+1)), n) for i in range(n-1) for j in range(i+1, n)])

# weak rectangulotope via the summed submodular function f_weakrect 
def wrp2(n):
    return submod_polytope(f_weakrect, n)

# yin rectangulotopes via yin shards
def yrp(n):
    return sum([submod_polytope(lambda S, m : f_yin(S, range(i, j+1), range(j+1, k+1)), n) for i in range(n-1) for j in range(i, n-1) for k in range(j+1, n)])

# yin rectangulotopes via the summed submod function
def yrp2(n):
    return submod_polytope(f_yinrect, n)

# yang polytopes, just the same
def yarp(n):
    return sum([submod_polytope(lambda S, m : f_yang(S, range(i, j+1), range(j+1, k+1)), n) for i in range(n-1) for j in range(i, n-1) for k in range(j+1, n)])

# strong rectangulotope
def srp(n):
    return submod_polytope(f_strongrect, n)

# examples
#WRP6 = wrp2(6)
#WRP6 == wrp(6)
#YRP4 = yrp(5)
#YRP4 == yrp2(5)







 
        
        


    






