# submodular functions for rectangulotopes

from itertools import *

## SUBMODULAR FUNCTIONS

# simplex for a Loday arc I
# 1 if S and I intersects, 0 otherwise
def f_loday(S, I):
    return int(bool(set(S) & set (I)))

# anti-simplex for an anti-Loday arc I
# -1 if I is a subset of S, 0 otherwise
def f_antiloday(S, I):
    return - int(set(I).issubset(set(S)))

# simplex + antisimplex (cuboctahedron in 3d)
# 1 if S and I intersects but I is not a subset of S, 0 otherwise
def f_cuboc(S, I):
    return int(bool(set(S) & set (I)) and not set(I).issubset(set(S)))

# weak rectangulation
# number of I such that S and I intersects but I is not a subset of S
def f_weakrect(S, n):
    return sum([int(bool(set(S) & set (range(i, j+1))) and not set(range(i, j+1)).issubset(set(S))) for i in range(n-1) for j in range(i+1, n)])

# yin shard polytopes
# 1 if and only if S intersects I but J is not a subset of S
def f_yin(S, I, J):
    return int(bool(set(S) & set(I)) and not (set(J).issubset(S)))

# yin rectangulations
# number of pairs I, J such that S intersects I but J is not a subset of S
def f_yinrect(S, n):
    return sum([int (bool(set(S) & set(range(i, j+1))) and not (set(range(j+1, k+1)).issubset(S))) for i in range(n-1) for j in range(i, n-1) for k in range(j+1, n)])

# yang shard polytopes
# 1 if and only if S intersects J but I is not a subset of S
def f_yang(S, I, J):
    return int(bool(set(S) & set(J)) and not (set(I).issubset(S)))

## GENERIC FUNCTIONS

# powerset recipe
def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

# compute the polyhedron for a given submod function
def submod_polytope(f, n):
    all_ineq = []
    for S in powerset(range(n)):
        all_ineq.append((f(S),) + tuple(-1 if x in S else 0 for x in range(n)))
    return Polyhedron(Polyhedron(ieqs=all_ineq).vertices()) # remove rays

## RECTANGULOTOPES

# associahedron as a Minkowski sum
def assoc(n):
    A = Polyhedron([(0,)*n])
    for i in range(n-1):
        for j in range(i+1, n):
            I = range(i, j+1)
            g_loday = lambda S : f_loday(S, I)
            A = A + submod_polytope(g_loday, n)
    return Polyhedron(A.vertices()) # remove rays

# weak rectangulation polytope (diagonal / Baxter)
# either sum all lodays and antilodays, or the cubocs
def wrp(n):
    A = Polyhedron([(0,)*n])
    for i in range(n-1):
        for j in range(i+1, n):
            I = range(i, j+1)
            # g_loday = lambda S : f_loday(S, range(i, j+1))
            # g_antiloday = lambda S : f_antiloday(S, range(i, j+1))
            # A = A + submod_polytope(g_loday, n) + submod_polytope(g_antiloday, n)
            g = lambda S : f_cuboc(S, I)
            A = A + submod_polytope(g, n)
    return Polyhedron(A.vertices()) # remove rays

# weak rectangulotope via the summed submodular function f_weakrect 
def wrp2(n):
    return Polyhedron(submod_polytope(lambda S : f_weakrect(S, n), n).vertices())

# yin rectangulotopes via yin shards
def yrp(n):
    A = sum([submod_polytope(lambda S : f_yin(S, range(i, j+1), range(j+1, k+1)), n) for i in range(n-1) for j in range(i, n-1) for k in range(j+1, n)])
    return Polyhedron(A.vertices()) # remove rays

# yin rectangulotopes via the summed submod function
def yrp2(n):
    return Polyhedron(submod_polytope(lambda S : f_yinrect(S, n), n).vertices())

# yang polytopes, just the same
def yarp(n):
    A = sum([submod_polytope(lambda S : f_yang(S, range(i, j+1), range(j+1, k+1)), n) for i in range(n-1) for j in range(i, n-1) for k in range(j+1, n)])
    return Polyhedron(A.vertices()) # remove rays


# examples
#WRP6 = wrp2(6)
#WRP6 == wrp(6)
#YRP4 = yrp(5)
#YRP4 == yrp2(5)







 
        
        


    






