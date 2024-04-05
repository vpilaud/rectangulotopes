### SOME EXPERIMENTS ON RECTANGULOTOPES

@cached_function
def diagonal_rectangulation_polytope(n, essential=False):
    return diagonal_rectangulations_congruence(n).quotientope(essential=essential)

@cached_function
def diagonal_rectangulation_polytope_vertices(n, essential=False):
    return diagonal_rectangulation_polytope(n).vertices()

def yinyangs(n):
    return Arcs(n, arcs=[arc for arc in Arcs(n, restriction='all') if len([i for i in range(arc._left + 1, arc._right) if arc._word[i] * arc._word[i+1] == -1]) == 1])

def yins(n):
    return Arcs(n, arcs=[arc for arc in yinyangs(n) if arc._word[arc._left+1] == 1])

def yangs(n):
    return Arcs(n, arcs=[arc for arc in yinyangs(n) if arc._word[arc._left+1] == -1])

def up_down_yins(n):
    return Arcs(n, arcs=Set(list(permutree_congruence([-1]*n)) + list(permutree_congruence([1]*n)) + list(yins(n))))

def up_down_yangs(n):
    return Arcs(n, arcs=Set(list(permutree_congruence([-1]*n)) + list(permutree_congruence([1]*n)) + list(yangs(n))))

def up_down_yinyangs(n):
    return Arcs(n, arcs=list(up_down_yins(n)) + list(up_down_yangs(n)))

def up_down_yins_quotientope(n, essential=False):
    return up_down_yins(n).quotientope(essential=essential)

def up_down_yangs_quotientope(n, essential=False):
    return up_down_yangs(n).quotientope(essential=essential)

def up_down_yinyangs_quotientope(n, essential=False):
    return up_down_yinyangs(n).quotientope(essential=essential)

def vertex_2_clumped_permutation(perm):
    n = len(perm)
    perm_inv = perm.inverse()
    res = vector(max(diagonal_rectangulation_polytope_vertices(n), key = lambda v: vector(v).dot_product(vector(perm_inv))))
    """
    for arc in yins(n):
        blacks = [arc._left] + [i for i in range(arc._left+1, arc._right) if arc._word[i] == 1]
        whites = [i for i in range(arc._left+1, arc._right) if arc._word[i] == -1] + [arc._right]
        v = max(blacks, key = lambda i: perm_inv(i+1)) + 1
        V = min(whites, key = lambda i: perm_inv(i+1)) + 1
        if perm_inv(v) > perm_inv(V):
            res = res + vector([0]*(v-1) + [1] + [0]*(V-v-1) + [-1] + [0]*(n-V))
    """
    for (v,V) in perm_inv.inversions():
        i = max([x for x in range(1,n+1) if x < v and perm_inv(x) > perm_inv(v)] + [0])
        k = min([x for x in range(1,n+1) if x > v and perm_inv(x) > perm_inv(v)] + [V])
        j = max([x for x in range(1,n+1) if x < V and perm_inv(x) < perm_inv(V)] + [v])
        l = min([x for x in range(1,n+1) if x > V and perm_inv(x) < perm_inv(V)] + [n+1])
        res = res + len([(x,y,z) for (x,y,z) in cartesian_product([range(i+1,v+1), range(j,k), range(V,l)]) if y-x >= 1 and z-y >= 2]) * vector([0]*(v-1) + [1] + [0]*(V-v-1) + [-1] + [0]*(n-V))
    return res

def up_down_yins_quotientope_bis(n):
    return Polyhedron([vertex_2_clumped_permutation(perm) for perm in up_down_yins(n).quotient()])

def vertex_2_clumped_permutation_bis(perm):
    perm_inv = perm.inverse()
    n = len(perm)
    S = perm.binary_search_tree() # careful, this has a left right inversion...
    T = perm.reverse().binary_search_tree()
    HS = [[]]; VS = [[]]; S.in_order_traversal(node_action=lambda node: VS.append([node.label()] + node[0].labels()) or HS.append([node.label()] + node[1].labels()))
    HT = [[]]; VT = [[]]; T.in_order_traversal(node_action=lambda node: HT.append([node.label()] + node[0].labels()) or VT.append([node.label()] + node[1].labels()))
    """
    Yins = [0] + [[0] + [max(0, len(HT[u]) * (min(max(VT[u]), v-1) - max(min(VS[v])-1, u) + 1) * len(HS[v]) - int((u+1) in VS[v]) * len(HS[v]) - int((v-1) in VT[u]) * len(HT[u]) + int(u == v-1)) for u in range(1, v)] for v in range(1, n+1)]
    Yangs = [0] + [[0] + [max(0, len(VS[u]) * (min(max(HS[u]), v-1) - max(min(HT[v])-1, u) + 1) * len(VT[v]) - int((u+1) in HT[v]) * len(VT[v]) - int((v-1) in HS[u]) * len(VS[u]) + int(u == v-1)) for u in range(1, v)] for v in range(1, n+1)]
    """
    Yins = [0] + [[0] + [max(0, len(HT[u]) * (min(max(VT[u]), v-1) - max(min(VS[v])-1, u) + 1) * len(HS[v])) for u in range(1, v)] for v in range(1, n+1)]
    Yangs = [0] + [[0] + [max(0, len(VS[u]) * (min(max(HS[u]), v-1) - max(min(HT[v])-1, u) + 1) * len(VT[v])) for u in range(1, v)] for v in range(1, n+1)]
    vertex = vector(range(binomial(n,2), -binomial(n,2)-1, -n)) # vector([len(HT[u]) * len(VT[u]) - len(HS[u]) * len(VS[u]) for u in range(1,n+1)]) # no need for Loday and antiLoday
    for v in range(1, n+1):
        for u in range(1, v):
            if perm_inv(u) > perm_inv(v):
                vertex = vertex + Yins[v][u] * vector([0]*(u-1) + [1] + [0]*(v-u-1) + [-1] + [0]*(n-v))
            else:
                vertex = vertex - Yangs[v][u] * vector([0]*(u-1) + [1] + [0]*(v-u-1) + [-1] + [0]*(n-v))
    return vertex

@cached_function
def up_down_yins_quotientope_ter(n):
    return Polyhedron([vertex_2_clumped_permutation_bis(perm) for perm in up_down_yins(n).quotient()])

@cached_function
def up_down_yangs_quotientope_ter(n):
    return Polyhedron([vertex_2_clumped_permutation_bis(perm) for perm in up_down_yangs(n).quotient()])

@cached_function
def up_down_yinyangs_quotientope_bis(n):
    return Polyhedron([vertex_2_clumped_permutation_bis(perm) for perm in rectangulations_congruence(n).quotient()])

### SOURCE AND TARGET TREES

@cached_function
def source_tree(perm):
    r"""
    Return the source tree of the permutation perm.
    Be careful: the vertical subtrees are left subtrees, the horizontal subtrees are right subtrees.

    EXAMPLES::
        sage: perm = Permutation([7, 5, 1, 14, 8, 6, 15, 11, 4, 2, 10, 9, 16, 13, 3, 12])
        sage: source_tree(perm)
        7[5[1[., 4[2[., 3[., .]], .]], 6[., .]], 14[8[., 11[10[9[., .], .], 13[12[., .], .]]], 15[., 16[., .]]]]
    """
    return perm.binary_search_tree()

@cached_function
def target_tree(perm):
    r"""
    Return the target tree of the permutation perm.
    Be careful: the horizontal subtrees are left subtrees, the vertical subtrees are right subtrees.

    EXAMPLES::
        sage: perm = Permutation([7, 5, 1, 14, 8, 6, 15, 11, 4, 2, 10, 9, 16, 13, 3, 12])
        sage: target_tree(perm)
        12[3[2[1[., .], .], 9[4[., 6[5[., .], 8[7[., .], .]]], 10[., 11[., .]]]], 13[., 16[15[14[., .], .], .]]]
    """
    return perm.reverse().binary_search_tree()

### BASIC OPERATIONS ON TREES

@cached_function
def tree_parents(tree):
    r"""
    Return the list whose ith entry is the parent of node i in the given tree.
    We put a 0 in front and at the root of the tree.
    """
    n = len(tree.labels())
    D = tree.as_digraph()
    return [0] + [(D.neighbors_in(i) + [0])[0] for i in range(1, n+1)]

@cached_function
def tree_children(tree):
    r"""
    Return the list whose ith entry is the pair of left and right children of node i in the given tree.
    We put a 0 in front and at the root of the tree.
    """
    res = [None]
    tree.in_order_traversal(node_action=lambda node: res.append((node[0].label(), node[1].label())))
    return res

@cached_function
def tree_child_side(tree):
    r"""
    Return the list whose ith entry is 'left' if i is the left child of its parent, 'right' if i is the left child of its parent, and 'root' if it is the root.
    """
    res = dict([(l, ['root']) for l in tree.labels()])
    tree.in_order_traversal(node_action=lambda node: (res[node[0].label()].append('left') if node[0] else False) or (res[node[1].label()].append('right') if node[1] else False))
    return dict([(l, res[l][-1]) for l in tree.labels()])

@cached_function
def tree_subtrees(tree):
    r"""
    Return the list whose ith entry is the subtree of the tree at node i.
    """
    res = [None]
    tree.in_order_traversal(node_action=lambda node: res.append(node))
    return res

@cached_function
def tree_left_corners(tree):
    r"""
    Return the left corners of the trees, that is the nodes with no left child, visited in pre-order.
    """
    res = []
    tree.in_order_traversal(node_action=lambda node: res.append(node.label()) if not node[0] else True)
    return res

@cached_function
def tree_right_corners(tree):
    r"""
    Return the right corners of the trees, that is the nodes with no right child, visited in post-order.
    """
    res = []
    tree.in_order_traversal(node_action=lambda node: res.append(node.label()) if not node[1] else True)
    return res[::-1]

@cached_function
def segments(perm):
    r"""
    Return the vertical and horizontal segments of the rectangulation of the permutation perm. This only depends on the weak class.

    EXAMPLES::
        sage: perm = Permutation([7, 5, 1, 14, 8, 6, 15, 11, 4, 2, 10, 9, 16, 13, 3, 12])
        sage: segments(perm)
        ([([1], [2, 4]),
          ([2], [3]),
          ([5], [6]),
          ([7], [8, 14]),
          ([4, 6, 8], [9, 10, 11]),
          ([3, 9, 10, 11], [12, 13]),
          ([14], [15]),
          ([15], [16])],
         [([2, 3], [4, 9]),
          ([1, 4], [5, 6]),
          ([5, 6], [7, 8]),
          ([9], [10]),
          ([10], [11]),
          ([12], [13]),
          ([8, 11, 13], [14, 15, 16])])
    """
    S = source_tree(perm).as_digraph()
    T = target_tree(perm).as_digraph()
    verticalSegments = []
    horizontalSegments = []
    for i in range(1, len(perm)):
        left = T.all_paths(i+1, i)
        right = S.all_paths(i, i+1)
        up = S.all_paths(i+1, i)
        down = T.all_paths(i, i+1)
        if left: verticalSegments.append((left[0][1:], right[0][1:][::-1]))
        if down: horizontalSegments.append((up[0][1:], down[0][1:][::-1]))
    return (verticalSegments, horizontalSegments)

@cached_function
def bounding_segments(perm):
    r"""
    Return the list whose ith entry gives the left, right, up, down segments bounding the rectangle i in the rectangulation of the permutation perm. This only depends on the weak class. The 0th entry is fixed to [0, 0, 0, 0].

    EXAMPLES::
        sage: perm = Permutation([7, 5, 1, 14, 8, 6, 15, 11, 4, 2, 10, 9, 16, 13, 3, 12])
        sage: bounding_segments(perm)
        [[0, 0, 0, 0],
         [0, ([1], [2, 4]), 0, ([1, 4], [5, 6])],
         [([1], [2, 4]), ([2], [3]), 0, ([2, 3], [4, 9])],
         [([2], [3]), ([3, 9, 10, 11], [12, 13]), 0, ([2, 3], [4, 9])],
         [([1], [2, 4]), ([4, 6, 8], [9, 10, 11]), ([2, 3], [4, 9]), ([1, 4], [5, 6])],
         [0, ([5], [6]), ([1, 4], [5, 6]), ([5, 6], [7, 8])],
         [([5], [6]), ([4, 6, 8], [9, 10, 11]), ([1, 4], [5, 6]), ([5, 6], [7, 8])],
         [0, ([7], [8, 14]), ([5, 6], [7, 8]), 0],
         [([7], [8, 14]),
          ([4, 6, 8], [9, 10, 11]),
          ([5, 6], [7, 8]),
          ([8, 11, 13], [14, 15, 16])],
         [([4, 6, 8], [9, 10, 11]),
          ([3, 9, 10, 11], [12, 13]),
          ([2, 3], [4, 9]),
          ([9], [10])],
         [([4, 6, 8], [9, 10, 11]),
          ([3, 9, 10, 11], [12, 13]),
          ([9], [10]),
          ([10], [11])],
         [([4, 6, 8], [9, 10, 11]),
          ([3, 9, 10, 11], [12, 13]),
          ([10], [11]),
          ([8, 11, 13], [14, 15, 16])],
         [([3, 9, 10, 11], [12, 13]), 0, 0, ([12], [13])],
         [([3, 9, 10, 11], [12, 13]), 0, ([12], [13]), ([8, 11, 13], [14, 15, 16])],
         [([7], [8, 14]), ([14], [15]), ([8, 11, 13], [14, 15, 16]), 0],
         [([14], [15]), ([15], [16]), ([8, 11, 13], [14, 15, 16]), 0],
         [([15], [16]), 0, ([8, 11, 13], [14, 15, 16]), 0]]
    """
    n = len(perm)
    res = [[0,0,0,0] for i in range(n+1)]
    for vs in segments(perm)[0]:
        for i in vs[1]:
            res[i][0] = vs
        for i in vs[0]:
            res[i][1] = vs
    for hs in segments(perm)[1]:
        for i in hs[1]:
            res[i][2] = hs
        for i in hs[0]:
            res[i][3] = hs
    return res

### CANONICAL EMBEDDING WEAK RECTANGULATION

def weak_insertion(perm):
    r"""
    Return the coordinates (left, right, down, up) of the rectangles of the weak insertion of perm.

        sage: perm = Permutation([7, 5, 1, 14, 8, 6, 15, 11, 4, 2, 10, 9, 16, 13, 3, 12])
        sage: weak_insertion(perm)
        ((0, 1, 12, 16),
         (1, 2, 13, 16),
         (2, 11, 13, 16),
         (1, 8, 12, 13),
         (0, 5, 10, 12),
         (5, 8, 10, 12),
         (0, 7, 0, 10),
         (7, 8, 3, 10),
         (8, 11, 7, 13),
         (8, 11, 6, 7),
         (8, 11, 3, 6),
         (11, 16, 4, 16),
         (11, 16, 3, 4),
         (7, 14, 0, 3),
         (14, 15, 0, 3),
         (15, 16, 0, 3))
    """
    n = len(perm)
    perm_inv = perm.inverse()
    return tuple([(max([0] + [j for j in range(1,i+1) if perm_inv(i) > perm_inv(j)]), min([n+1] + [j for j in range(i+1,n+1) if perm_inv(i) < perm_inv(j)])-1, n-min([n+1] + [j for j in range(i+1,n+1) if perm_inv(i) > perm_inv(j)])+1, n-max([0] + [j for j in range(1,i+1) if perm_inv(i) < perm_inv(j)])) for i in range(1,n+1)])

def check_weak_insertion(n):
    r"""
    Check consistency of the weak insertion of all permutations of size n.

    EXAMPLES::
        sage: check_strong_insertion(6)
        ok
    """
    rectangulations = []
    count = 0
    for cl in diagonal_rectangulations_congruence(n).congruence_classes().values():
        if len(Set([weak_insertion(perm) for perm in cl])) != 1:
            print("not well-defined", cl)
        rectangulation = weak_insertion(cl[0])
        if not all([r1 <= l2 or r2 <= l1 or t1 <= b2 or t2 <= b1 for ((l1,r1,b1,t1),(l2,r2,b2,t2)) in Subsets(rectangulation,2)]):
            print("not disjoint")
        if add((r-l) * (t-b) for (l,r,b,t) in rectangulation) != n*n:
            print("not covering")
        rectangulations.append(rectangulation)
        count = count + 1
    if len(Set(rectangulations)) != count:
        print("not injective")
    print("ok")

### CANONICAL EMBEDDING STRONG RECTANGULATION

# We first need to find the stones in all directions.

def up_stones(perm, i):
    perm_inv = perm.inverse()
    res = [i]
    for j in tree_left_corners(tree_subtrees(source_tree(perm))[i]):
        up_segment = bounding_segments(perm)[j][2]
        if not up_segment:
            return res
        x = [u for u in up_segment[0] if perm_inv(u) > perm_inv(j)]
        if len(x) > 1:
            return res + up_stones(perm, x[1])
    return res

def down_stones(perm, i):
    perm_inv = perm.inverse()
    res = [i]
    for j in tree_right_corners(tree_subtrees(target_tree(perm))[i]):
        down_segment = bounding_segments(perm)[j][3]
        if not down_segment:
            return res
        x = [d for d in down_segment[1] if perm_inv(d) < perm_inv(j)]
        if len(x) > 1:
            return res + down_stones(perm, x[-2])
    return res

def right_stones(perm, i):
    perm_inv = perm.inverse()
    res = [i]
    for j in tree_right_corners(tree_subtrees(source_tree(perm))[i]):
        right_segment = bounding_segments(perm)[j][1]
        if not right_segment:
            return res
        x = [r for r in right_segment[1] if perm_inv(r) > perm_inv(j)]
        if len(x) > 1:
            return res + right_stones(perm, x[-2])
    return res

def left_stones(perm, i):
    perm_inv = perm.inverse()
    res = [i]
    for j in tree_left_corners(tree_subtrees(target_tree(perm))[i]):
        left_segment = bounding_segments(perm)[j][0]
        if not left_segment:
            return res
        x = [l for l in left_segment[0] if perm_inv(l) < perm_inv(j)]
        if len(x) > 1:
            return res + left_stones(perm, x[1])
    return res

def strong_insertion(perm):
    r"""
    Return the coordinates (left, right, down, up) of the rectangles of the strong insertion of perm.

    EXAMPLES::
        sage: perm = Permutation([7, 5, 1, 14, 8, 6, 15, 11, 4, 2, 10, 9, 16, 13, 3, 12])
        sage: strong_insertion(perm)
        ((0, 1, 9, 16),
         (1, 2, 13, 16),
         (2, 25/2, 13, 16),
         (1, 8, 9, 13),
         (0, 5, 6, 9),
         (5, 8, 6, 9),
         (0, 7, 0, 6),
         (7, 8, 3, 6),
         (8, 25/2, 10, 13),
         (8, 25/2, 7, 10),
         (8, 25/2, 3, 7),
         (25/2, 16, 11, 16),
         (25/2, 16, 3, 11),
         (7, 23/2, 0, 3),
         (23/2, 15, 0, 3),
         (15, 16, 0, 3))
    """
    n = len(perm)
    res = [[0, n, 0, n] for i in range(n+1)]
    for vs in segments(perm)[0]:
        ds = down_stones(perm, vs[0][0])
        us = up_stones(perm, vs[1][-1])
        left = Set(range(1, vs[0][0]))
        right = Set(range(vs[1][-1]+1, n+1))
        x = 0
        for d in ds:
            diff = Set(tree_subtrees(target_tree(perm))[d].labels())
            left = left.union(diff)
            right = right.difference(diff)
            x = x + 1/2
        for u in us:
            diff = Set(tree_subtrees(source_tree(perm))[u].labels())
            left = left.difference(diff)
            right = right.union(diff)
            x = x - 1/2
        x = x + len(left)
        # print(perm, vs, ds, us, left, right, x)
        for i in vs[0]:
            res[i][1] = x
        for i in vs[1]:
            res[i][0] = x
    for hs in segments(perm)[1]:
        rs = right_stones(perm, hs[0][0])
        ls = left_stones(perm, hs[1][-1])
        up = Set(range(1, hs[0][0]))
        down = Set(range(hs[1][-1]+1, n+1))
        x = 0
        for r in rs:
            diff = Set(tree_subtrees(source_tree(perm))[r].labels())
            up = up.union(diff)
            down = down.difference(diff)
            x = x - 1/2
        for l in ls:
            diff = Set(tree_subtrees(target_tree(perm))[l].labels())
            up = up.difference(diff)
            down = down.union(diff)
            x = x + 1/2
        x = x + len(down)
        # print(perm, hs, rs, ls, up, down, x)
        for i in hs[0]:
            res[i][2] = x
        for i in hs[1]:
            res[i][3] = x
    return tuple([tuple(t) for t in res[1:]])

def check_strong_insertion(n):
    r"""
    Check consistency of the strong insertion of all permutations of size n.

    EXAMPLES::
        sage: check_strong_insertion(7)
        ok
    """
    rectangulations = []
    count = 0
    for cl in rectangulations_congruence(n).congruence_classes().values():
        # if len(Set([strong_insertion(perm) for perm in cl])) != 1:
        #     print("not well-defined", cl)
        rectangulation = strong_insertion(cl[0])
        if not all([r1 <= l2 or r2 <= l1 or t1 <= b2 or t2 <= b1 for ((l1,r1,b1,t1),(l2,r2,b2,t2)) in Subsets(rectangulation,2)]):
            print("not disjoint", cl[0])
        if add((r-l) * (t-b) for (l,r,b,t) in rectangulation) != n*n:
            print("not covering")
        rectangulations.append(rectangulation)
        count = count + 1
    if len(Set(rectangulations)) != count:
        print("not injective")
    print("ok")

### DRAWING

def draw_rectangulation(rectangulation):
    r"""
    Draw the rectangulations whose coordinates are given by (left, right, down, up).

    EXAMPLES::
        sage: perm = Permutation([7, 5, 1, 14, 8, 6, 15, 11, 4, 2, 10, 9, 16, 13, 3, 12])
        sage: draw_rectangulation(strong_insertion(perm))
        Launched png viewer for Graphics object consisting of 16 graphics primitives
    """
    res = []
    for (l, r, d, u) in rectangulation:
        res.append(polygon2d([[l,d], [l,u], [r,u], [r,d]], fill=None, axes=False))
    return add(res)

def perm_rectangulation_to_tikz(perm, print_permutation=False):
    r"""
    Return a tikz code for rectangulations.

    EXAMPLES::
        sage: perm = Permutation([7, 5, 1, 14, 8, 6, 15, 11, 4, 2, 10, 9, 16, 13, 3, 12])
        sage: print(to_tikz(strong_insertion(perm)))
        \begin{tikzpicture}[very thick]
            \draw (0,8) rectangle (1,16);
            \draw (1,13) rectangle (2,16);
            \draw (2,13) rectangle (13,16);
            \draw (1,8) rectangle (8,13);
            \draw (0,5) rectangle (5,8);
            \draw (5,5) rectangle (8,8);
            \draw (0,0) rectangle (7,5);
            \draw (7,3) rectangle (8,5);
            \draw (8,10) rectangle (13,13);
            \draw (8,7) rectangle (13,10);
            \draw (8,3) rectangle (13,7);
            \draw (13,12) rectangle (16,16);
            \draw (13,3) rectangle (16,12);
            \draw (7,0) rectangle (11,3);
            \draw (11,0) rectangle (15,3);
            \draw (15,0) rectangle (16,3);
        \end{tikzpicture}
    """
    n = len(perm)
    res = '\\begin{tikzpicture}[very thick, scale=.5]\n'
    for (l, r, d, u) in strong_insertion(perm):
        res = res + '    \\draw (' + str(l) + ',' + str(d) + ') rectangle (' + str(r) + ',' + str(u) + ');\n'
    if print_permutation:
        res = res + '    \\node at (' + str(.5 * n) + ',-.5) {' + ''.join([str(perm(i)) for i in range(1, n+1)]) + '};\n'
    res = res + '\\end{tikzpicture}\n'
    return res

def save_tikz(perm, filename, print_permutation=False):
    r"""
    Save the tikz file of the rectangulation given by the strong insertion of perm.

    EXAMPLES::
        sage: n = 6
        ....: for perm in rectangulations_congruence(n).quotient():
        ....:     save_tikz(perm, "../rectangulotopes/programmation/tikzRectangulations/rectangulations" + str(n) + "/rectangulation_" + ''.join([str(perm(i)) for i in range(1, n+1)]) + ".tikz")
        ....:     save_tikz(perm, "../rectangulotopes/programmation/tikzRectangulations/rectangulations" + str(n) + "/rectangulations" + str(n) + ".tikz")
        ....:
    """
    f = open(filename, "a")
    f.write(to_tikz(perm, print_permutation=print_permutation))
    f.close()
