# submodular functions for rectangulotopes

from itertools import *

# simplex for a Loday arc I
# 1 if S and I intersects, 0 otherwise
def f_loday(S, I):
    return int(bool(set(S) & set (I)))

# anti-simplex for an anti-Loday arc I
# -1 if I is a subset of S, 0 otherwise
def f_antiloday(S, I):
    return - int(bool(set(I).issubset(set(S))))

# simplex + antisimplex
# 1 if S and I intersects but I is not a subset of S, 0 otherwise
def f_cuboc(S, I):
    return int(bool(set(S) & set (I)) and not bool(set(I).issubset(set(S))))

# weak rectangulation
# number of I such that S and I intersects but I is not a subset of S
def f_weakrect(S, n):
    return sum([int(bool(set(S) & set (range(i, j+1))) and not bool(set(range(i, j+1)).issubset(set(S)))) for i in range(n-1) for j in range(i+1, n)])

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
    return Polyhedron(ieqs=all_ineq)

# test
n = 6
g = lambda S : f_antiloday(S, range(n-2))
print(submod_polytope(g, n))

# associahedron as a Minkowski sum
def assoc(n):
    A = Polyhedron([(0,)*n])
    for i in range(n-1):
        for j in range(i+1, n):
            I = range(i, j+1)
            g_loday = lambda S : f_loday(S, range(i, j+1))
            A = A + submod_polytope(g_loday, n)
    return Polyhedron(A.vertices()) # remove rays

# diagonal/Baxter/weak rectangulation polytope
# either sum all lodays and antilodays, or the cubocs
def wrp(n):
    A = Polyhedron([(0,)*n])
    for i in range(n-1):
        for j in range(i+1, n):
            I = range(i, j+1)
            # g_loday = lambda S : f_loday(S, range(i, j+1))
            # g_antiloday = lambda S : f_antiloday(S, range(i, j+1))
            # A = A + submod_polytope(g_loday, n) + submod_polytope(g_antiloday, n)
            g = lambda S : f_cuboc(S, range(i, j+1))
            A = A + submod_polytope(g, n)
    return Polyhedron(A.vertices()) # remove rays

# or use the summed submodular 
def wrp2(n):
    return Polyhedron(submod_polytope(lambda S : f_weakrect(S, n), n).vertices())




 
        
        


    






