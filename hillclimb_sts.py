from itertools import combinations
from sts_utils import *
from copy import deepcopy


''''
'''
def is_Pasch(P):
    S = set()
    for t in P:
        S = S.union(t)
    if len(S) != 6:
        return False
    return True


def is_CA(T):
    S = set()
    for t in T:
        S = S.union(t)
    if len(S) != 8:
        return False
    # counting degrees
    deg = {}
    for t in T:
        for i in range(3):
            if t[i] not in deg.keys():
                deg[t[i]] = 1
            else:
                deg[t[i]] += 1
    if sorted(deg.values()) != [3,2,2,2,2,2,1,1]:
        return False
    return True

def element_degrees(P, S = None):
    deg = {}
    for s in S:
        deg[s] = 0
    for B in P:
        for i in range(3):
            deg[B[i]] += 1
    deglist = sorted(deg.values())
    return deglist

def block_degrees(P):
    deglist = []
    for B in P:
        P.remove(B)
        deglist.append(len(filter(lambda T: bool(set(T) & set(B)), P)))
        P.append(B)
    return sorted(deglist)


'''
   Testing if new triple is suitable
'''
def hts_test_triple(P, t):
    # Test no Pasch and no C_A configurations are built with the new triple
    # Get only non-empty intersections
    nonvoid = list(filter(lambda b: bool(set(t) & set(b)), P))
    comb = combinations(range(len(nonvoid)), 3)
    for c in comb:
        SS = [nonvoid[c[0]], nonvoid[c[1]], nonvoid[c[2]], t]
        S = set()
        for B in SS:
            S = S.union(B)
        if len(S) == 6:            # Pasch found
            return False
        if len(S) > 8:
            continue
        [P.remove(B) for B in [nonvoid[c[0]], nonvoid[c[1]], nonvoid[c[2]]]]
        for B in P:
            SS.append(B)
            Sc = set(S)
            Sc = Sc.union(B)
            deg = {}
            if len(Sc) != 8:
                SS.remove(B)
                continue
            else:
                if element_degrees(SS, Sc) == [1,1,2,2,2,2,2,3] and block_degrees(SS) == [2,3,3,4,4]:
                    return False
            SS.remove(B)
        [P.append(B) for B in [nonvoid[c[0]], nonvoid[c[1]], nonvoid[c[2]]]]
    return True


'''
   Initialize live points dictionary
'''
def init_live_points(X, P):
    live = dict([(i, set(X)) for i in X])     # Dictionary of live_points
    for t in P:
        x, y, z = t
        live[x].remove(y), live[x].remove(z)
        live[y].remove(x), live[y].remove(z)
        live[z].remove(x), live[z].remove(y)
    for i in X:
        live[i].remove(i)
        if len(live[i]) == 0:
            del live[i]
    return live

'''
   Updating live points after adding B0 and removing B1
'''
def update_live_points(live, B0, B1):
    if B1 is not None:
        x, y, z = B1
        if x not in live.keys(): live[x] = set()
        if y not in live.keys(): live[y] = set()
        if z not in live.keys(): live[z] = set()
        live[x].add(y), live[x].add(z)
        live[y].add(x), live[y].add(z)
        live[z].add(x), live[z].add(y)
    x, y, z = B0
    live[x].remove(y), live[x].remove(z)
    live[y].remove(x), live[y].remove(z)
    live[z].remove(x), live[z].remove(y)
    if len(live[x]) < 2: del live[x]
    if len(live[y]) < 2: del live[y]
    if len(live[z]) < 2: del live[z]
    return

'''
   Hill-climbing algorithm for Hall triple systems
'''
def hts_hill_climb(k, base = [], maxiters = None):
    v = 3 ** k
    b = (v ** 2 - v) // 6  # Complete number of blocks
    P = deepcopy(base)
    X = list(range(v))              # Set X
    # Initialize live points
    live_points = init_live_points(X, P)
    # Iterate until completing
    if maxiters is None:
        maxiters = 2**32
    iters = 0
    while len(P) < b and iters < maxiters:
        x = np.random.choice(list(live_points.keys()))
        y, z = sorted(np.random.choice(list(live_points[x]), 2, replace=False))
        B0 = (x,y,z)
        B1 = get_block_with_pair(P, (y, z))
        if B1 in base:
            continue
        if B1 is not None:
            P.remove(B1)
        if hts_test_triple(P, B0):
            P.append(tuple(sorted(B0)))
            update_live_points(live_points, B0, B1)
        elif B1 is not None:
            P.append(B1)
        iters += 1
#        if iters % 100 == 0:
#            print(iters)
    return P


#np.random.seed(2)
k = 3
v = 3**k
n = v * (v-1) // 6
B = []
'''T = affine_sts(2).T
B = deepcopy(T)
B.extend([(t[0] + 9, t[1] + 9, t[2] + 9) for t in T])
B.extend([(t[0] + 18, t[1] + 18, t[2] + 18) for t in T])'''
S = []
while len(S) < n:
    S = hts_hill_climb(k, base = B, maxiters=600)
    print(len(S))
S = SteinerTripleSystem(3**k, S)
print(S)
print(cycle_vector(S))