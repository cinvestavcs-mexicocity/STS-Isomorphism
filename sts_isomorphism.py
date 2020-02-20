from Quasigroup import *
from sts_utils import *
from itertools import *
import time

'''
'''
def quasigroup_from_sts(S):
    # Ordered pairs only
    X = S.X[:]
    mult_table = {}
    for t in S.T:
        mult_table[(t[0], t[1])] = t[2]
        mult_table[(t[0], t[2])] = t[1]
        mult_table[(t[1], t[2])] = t[0]
    return SteinerQuasigroup(X, mult_table)


'''
Generate while constructing the whole image
'''
def verify_quasigroup_isomorphism(Q1, Q2, f):
    d = dict(f)
    R = set([e[0] for e in f])
    new_elements = set([f[0][0]])
    while len(d) < len(Q1.X):
        y = new_elements.pop()
        for x in R:
            z = Q1.operate(x, y)
            w = Q2.operate(d[x], d[y])
            if z not in d.keys():
                if w in d.values():
                    return False
                d[z] = w
                new_elements.add(z)
            elif d[z] != w:
                return False
        R.add(y)
    return d

'''
def extend_generators(Q, gens):
    R = set(Q.X[:])
    H = Q.generate(gens)
    R = R.difference(H.X)
    if len(R) == 0:
        return Q
    gens.append(np.random.choice(list(R)))
    return H'''

'''
'''
def find_quasigroup_generators(Q):
    # Random pair to start
    X = set(Q.X[:])
    H = SteinerQuasigroup()
    G = []
    while H.order != Q.order:
        y = np.random.choice(list(X))
        G.append(y)
        H = Q.extend(H, y)
        X = X.difference(H.X)
    return sorted(G)

'''
    Returns a permutation defining the isomorphism
def miller_algorithm(Q1, Q2):
    if type(Q1) == type(SteinerTripleSystem(7)) and type(Q2) == type(SteinerTripleSystem(7)):
        Q1, Q2 = quasigroup_from_sts(Q1), quasigroup_from_sts(Q2)
    elif not (type(Q1) == type(SteinerQuasigroup()) and type(Q2) == type(SteinerQuasigroup())):
        raise ValueError('must be quasigroup or STS')
    return find_quasigroup_isomorphism(Q1, Q2)'''

def sts_miller(S1, S2):
    #start = time.time()
    Q1 = quasigroup_from_sts(S1)
    Q2 = quasigroup_from_sts(S2)
    #print('To quasi: {}'.format(time.time() - start))
    # Check orders
    if Q1.order != Q2.order:
        return False
    # Find a set of generators for Q1
    #start = time.time()
    G1 = find_quasigroup_generators(Q1)
    #print('Gens: {}'.format(time.time() - start))
    # For every set of m elements check if it is a well defined isomorphism
    # Iterate over all possible permutations
    for G2 in permutations(Q2.X, len(G1)):
        iso = list(zip(G1, G2))
        #start = time.time()
        d = verify_quasigroup_isomorphism(Q1, Q2, iso)
        #print('Verify: {}'.format(time.time() - start))
        if d:
            phi = list(map(lambda i: d[i], range(Q1.order)))
            if S1.permute(phi) == S2:
                return phi
    return False
    #return miller_algorithm(Q1, Q2)

'''
    Tries one round of Alg.1 to complete a partial bijection.
    input
        S0, S1: Steiner triple systems of same length
        partial_phi: partially defined bijection represented with a list. Not defined elements have value -1
    returns
        an extended bijection phi
        FALSE if the bijection cannot be extended to an isomorphism between STSs
'''
def try_isomorphism_round(S0, S1, partial_phi):
    n = len(S0.X)
    phi = partial_phi[:]
    for x in range(n):
        for y in range(x+1, n):
            # Verify if bijection is defined on pair (x,y), otherwise try with a new pair
            u, v = phi[x], phi[y]
            if u != -1 and v != -1:
                # Find the triplet containing (x, y)
                t = S0.get_triplet_with_pair(x, y)
                #print('xy found with {}, {} : {}'.format(x, y, z))
                # Get the third element
                z = list(filter(lambda z: z not in [x,y], t))[0]
                # if the third element is undefined, then define it with the triple (u, v, w) in S1
                w = phi[z]
                if w == -1:
                    t = S1.get_triplet_with_pair(u, v)
                    #print('ab found with {}, {}: {}'.format(a, b, t))
                    partial_phi[z] = list(filter(lambda x: x not in [u,v], t))[0]
                else:
                    # If phi was defined on w, verify that definition is an isomorphism
                    if tuple(sorted([u, v, w])) not in S1.T:
                        return False
        #return partial_phi
        #if -1 not in partial_phi:
    return partial_phi

def try_isomorphism(S0, S1, partial_phi):
    phi = try_isomorphism_round(S0, S1, partial_phi)
    while phi != partial_phi:
        partial_phi = phi[:]
        phi = try_isomorphism_round(S0, S1, phi)
    return phi


def get_max_cycle(S, iters = 0):
    pairs = combinations(S.X, 2)
    n = S.order
    maxlen = 0
    maxG = None
    maxpair = None
    if iters == 0: iters = n * (n-1) // 2
    i = 0
    while i < iters:
        pair = next(pairs)
        Gij = S.cycle_graph(pair[0], pair[1])
        if len(Gij[0]) > maxlen:
            maxlen = len(Gij[0])
            maxG = deepcopy(Gij)
            maxpair = pair
        # If the cycle is at least half the size of the maximum length, stop
        if maxlen >= (S.order - 3) // 2:
            break
        i += 1
    return [maxpair, sorted(maxG, key=len, reverse=1)]


def generators_from_cycle(S, xy, Q = None):
    x, y = xy
    Gxy = S.cycle_graph(x, y)
    if Q is None:
        Q = quasigroup_from_sts(S)
    R = [x, y, Gxy[0][0]]
    H = Q.generate(R)
    i = 1
    while H.order != Q.order:
        if i >= len(Gxy):
            print('')
        if Gxy[i][0] in H.X:
            i += 1
        else:
            v = Gxy[i][0]
            R.append(v)
            H = Q.extend(H, v)
    return R



'''
'''
def sts_improved_miller(S1, S2, cycles_iters = 0):
    Q1, Q2 = quasigroup_from_sts(S1), quasigroup_from_sts(S2)

    start = time.time()
    (x, y), G = get_max_cycle(S1, cycles_iters)
    #print('Max cycle: {}'.format(time.time() - start))
    start = time.time()
    gens1 = generators_from_cycle(S1, (x, y), Q1)
    m = len(gens1)
    #print('Cycle gens: {}'.format(time.time() - start))
    cl1 = tuple(sorted(map(len, G)))

    # Iterate over isomorphic cycle graphs in S2
    for pair in combinations(S1.X, 2):
        start = time.time()
        i, j = pair
        Gij = S2.cycle_graph(i, j)
        cl2 = tuple(sorted(map(len, Gij)))
        if cl1 == cl2:
            genset = [C[0] for C in Gij]
            for perm in permutations(genset, len(gens1)-2):
                gens2 = [i, j]
                gens2.extend([idx for idx in perm])
                iso = list(zip(gens1, gens2))
                d = verify_quasigroup_isomorphism(Q1, Q2, iso)
                if d:
                    phi = list(map(lambda x: d[x], range(Q1.order)))
                    if S1.permute(phi) == S2:
                        return phi
                else:
                    gens2[0] = j
                    gens2[1] = i
                    iso = list(zip(gens1, gens2))
                    d = verify_quasigroup_isomorphism(Q1, Q2, iso)
                    if d:
                        phi = list(map(lambda x: d[x], range(Q1.order)))
                        if S1.permute(phi) == S2:
                            return phi
        #print('gens_from_cycle: ', time.time() - start)

        '''V = [C[0] for C in Gij]             # A vertex from each cycle
            perms = combinations(V, m - 2)
            for perm in perms:
                gens2 = [i, j]
                gens2.extend(perm)
                #print('Cycle gens: {}'.format(time.time() - start))
                iso = list(zip(gens1, gens2))
                d = verify_quasigroup_isomorphism(Q1, Q2, iso)
                if d:
                    phi = list(map(lambda i: d[i], range(Q1.order)))
                    if S1.permute(phi) == S2:
                        return phi

                gens2 = [j, i]
                gens2.extend(perm)
                # print('Cycle gens: {}'.format(time.time() - start))
                iso = list(zip(gens1, gens2))
                d = verify_quasigroup_isomorphism(Q1, Q2, iso)
                if d:
                    phi = list(map(lambda i: d[i], range(Q1.order)))
                    if S1.permute(phi) == S2:
                        return phi
            print('{}: {}'.format(pair, time.time() - start))'''

        '''
            # Extend over all possibilities
            gens2 = generators_from_cycle(S2, pair, Q2)
            perms = permutations(gens2)
            for perm in perms:
                # One element per cycle
                # possible iterations of initial gens
                iso = list(zip(gens1, perm))
                #start = time.time()
                phi = verify_quasigroup_isomorphism(Q1, Q2, iso)
                #print('Iso verify: {}'.format(time.time() - start))
                if phi and S1.permute(phi) == S2:
                    return phi'''
    return False


    '''# First part: Find the graph with longest cycle
    n = S1.order
    maxlen = 0
    maxC = None
    cl = None
    Q1, Q2 = quasigroup_from_sts(S1), quasigroup_from_sts(S2)
    for pair in combinations(S1.X, 2):
        Gij = S1.cycle_graph(pair[0], pair[1])
        for C in Gij:
            if len(C) > maxlen:
                maxlen = len(C)
                maxC = C
                cl = tuple(sorted(map(len, Gij)))
            # If the cycle is at least half the size of the maximum length, stop
            if maxlen >= (n - 3) // 2:
                break

    # Start generators from long cycle
    H = Quasigroup()
    gens = C[:3]
    while H.order < Q1.order:
        H = extend_generators(Q1, gens)'''



''' # Second part: Find an isomorphic graph in S2
    for pair in combinations(S1.X, 2):
        Gij = S1.cycle_graph(pair[0], pair[1])
        if cl == tuple(sorted(map(len, Gij))):
            for C in Gij:
                if len(C) == maxlen:'''





# import sts
# # All together
#
# '''
#     Tries one round of Alg.1 to complete a partial bijection.
#     input
#         S0, S1: Steiner triple systems of same length
#         partial_phi: partially defined bijection represented with a list. Not defined elements have value -1
#     returns
#         an extended bijection phi
#         FALSE if the bijection cannot be extended to an isomorphism between STSs
# '''
# def try_isomorphism_round(S0, S1, partial_phi):
#     n = len(S0.X)
#     phi = partial_phi[:]
#     for x in range(n):
#         for y in range(x+1, n):
#             # Verify if bijection is defined on pair (x,y), otherwise try with a new pair
#             u, v = phi[x], phi[y]
#             if u != -1 and v != -1:
#                 # Find the triplet containing (x, y)
#                 t = S0.get_triplet_with_pair(x, y)
#                 #print('xy found with {}, {} : {}'.format(x, y, z))
#                 # Get the third element
#                 z = list(filter(lambda x: x not in [x,y], t))[0]
#                 # if the third element is undefined, then define it with the triple (u, v, w) in S1
#                 if w == -1:
#                     t = S1.get_triplet_with_pair(u, v)
#                     #print('ab found with {}, {}: {}'.format(a, b, t))
#                     partial[z] = list(filter(lambda x: x not in [u,v], t))[0]
#                 else:
#                     # If phi was defined on w, verify that definition is an isomorphism
#                     if sorted([u, v, w]) not in S1.T:
#                         return False
#         if -1 not in partial:
#             return partial
#     return False
#
# '''
#     Tries to complete a bijection until:
#         - phi is completed
#         - phi extension remains the same
#         - phi cannot be extended to a STS isomorphism
#     input
#         S0, S1: Steiner triple systems of same length
#         partial_phi: partially defined bijection represented with a list. Not defined elements have value -1
#     returns
#         an extended bijection phi
#         FALSE if the bijection cannot be extended to an isomorphism between STSs
# '''
# def try_isomorphism(S0, S1, partial_phi):
#     phi = try_isomorphism_round(S0, S1, partial_phi)
#     while phi != partial_phi:
#         partial_phi = phi[:]
#         phi = try_isomorphism_round(S0, S1, phi)
#     return phi
#
# '''
#     Find all possible isomorphisms between S0 an S1 using pivots of the cycle graphs Swap0 and Swap1
#     Swap0 an Swap1 are triples of the form: (Cycles, a, b)
# '''
# def get_isomorphism_from_cycles(S0, S1, Swap0, Swap1):
#     # Verify if the cycle graphs are isomorphic (a simple count must work)
#     Gab, a, b = Swap0
#     Gcd, c, d = Swap1
#     lg = len(Gab)
#     n = len(S0.n)
#     isomorphisms = []
#     lens1, lens2 = list(map(len, Gab)), list(map(len, Gcd))
#     if sorted(lens1) != sorted(lens2):
#         return isomorphisms
#
#     idx1 = np.argsort(lens1)[::-1]
#     idx2 = np.argsort(lens2)[::-1]
#     # We start with longest cycle
#     for i in range(lg):
#         C0 = Gab[idx1[i]]
#         # Pivots of the first cycle
#         P0 = np.array(list(map(lambda t: t[1], C0)))
#         for j in range(lg):
#             C1 = Gcd[j]
#             cl = len(C1)
#             if cl != len(C0):
#                 break
#             # Try with every element in the cycle (forward and backward)
#             for k in range(cl):
#                 phi = [-1] * n
#                 # Pivots of the second cycle
#                 P1 = np.roll(list(map(lambda t: t[1], C1)), k)
#                 for l in range(len(P0)):
#                     phi[P0[l]] = P1[l]
#                 #print(phi)
#                 phi = try_isomorphism(S0, S1, phi)
#                 if phi:
#                     isomorphisms.append(tuple(phi))
#
#                 phi = [-1] * n
#                 # Pivots of the second cycle reversed
#                 P1 = P1[::-1]
#                 for l in range(len(P0)):
#                     phi[P0[l]] = P1[l]
#                 #print(phi)
#                 phi = try_isomorphism(S0, S1, phi)
#                 if phi:
#                     isomorphisms.append(tuple(phi))
#
#     return isomorphisms
#
# '''
#     Find all isomorphisms between Triples S0 and S1 by creating permutations using the cycle graphs
# '''
# def find_isomorphisms(S0, S1):
#     n = len(S0.X)
#     isomorphisms = set()
#     # For every graph
#     for a in range(n-1):
#         for b in range(a+1, n):
#             # Get the cycles of the corresponding pair
#             Gab = S0.cycle_graph(a, b)
#             # We are checking for every graph in the second system
#             for c in range(n-1):
#                 for d in range(c+1, n):
#                     Gcd = S1.cycle_graph(c, d)
#                     isomorphisms.update(set(get_isomorphism_from_cycles(S0, S1, (Gab, a, b), (Gcd, c, d))))
#     return sorted(list(isomorphisms))
#
# '''
#     Find any isomorphism. We use our method
# '''
# def find_isomorphisms(S0, S1):
#     n = len(S0.X)
#     isomorphisms = set()
#     '''
#         Look for a big cycle graph. Longest is of length n-3. Half that size might suffice
#         (which is the min required? Testing with (n-3))'''
#     lmax, l, r = 0, 0, (n-3) // 2
#     for a in range(n-1):
#         for b in range(a+1, n):
#             # Get the cycles of the corresponding pair
#             Gab = S0.cycle_graph(a, b)
#             l = max(map(len, Gab.pivots))
#             if l > lmax:
#                 lmax = l
#             # We are checking for every graph in the second system
#             for c in range(n-1):
#                 for d in range(c+1, n):
#                     Gcd = S1.cycle_graph(c, d)
#                     isomorphisms.update(set(get_isomorphism_from_cycles(S0, S1, (Gab, a, b), (Gcd, c, d))))
#     return sorted(list(isomorphisms))
#
# '''
#     Returns the first isomorphism found
# '''
# def find_an_isomorphism(S0, S1):
#     n = len(S0.X)
#     # For every graph
#     for a in range(n-1):
#         for b in range(a+1, n):
#             # Get the cycles of the corresponding pair
#             Gab = S0.cycle_graph(a, b)
#             # We are checking for every graph in the second system
#             for c in range(n-1):
#                 for d in range(c+1, n):
#                     Gcd = S1.cycle_graph(c, d)
#                     isomorphisms = get_isomorphism_from_cycles(S0, S1, (Gab, a, b), (Gcd, c, d))
#                     if len(isomorphisms) > 0:
#                         return isomorphisms[0]
#     return False
'''
S = SteinerTripleSystem(15)
gens = S.X[:3]
Q = quasigroup_from_sts(S)
H = SteinerQuasigroup()
for y in gens:
    H = Q.extend(H, y)'''