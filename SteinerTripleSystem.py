import numpy as np, random
from CycleGraph import *

def is_sts(X, T):
    v = len(X)
    k = len(T)
    n = v*(v-1)//6
    if k != n:
        raise ValueError('The number of triples is incorrect: ' + str(k) + ' != ' + str(n))
    S = set()
    for t in T:
        S.add(frozenset([t[0], t[1]]))
        S.add(frozenset([t[0], t[2]]))
        S.add(frozenset([t[1], t[2]]))
    if len(S) != v*(v-1)//2:
        #print S
        return False
    return True

''' SteineTripleSystem class.
'''
class SteinerTripleSystem(object):
    ''' Initializes a STS of order n. If a set of triples T is given, a test is run to check if it defines a valid STS.
        If the triples are not given then they are generated
    '''
    def __init__(self, n, T = None):
        self.order = n
        self.X = range(n)
        if (T == None):
            self.__generate__()
        elif (is_sts(self.X, T)):
            self.T = [tuple(sorted(t)) for t in T]
        else:
            raise ValueError('Not a valid STS')

    ''' Generation of a triple system with a given algorithm
    '''
    def __generate__(self):
        n = self.order
        if n % 6 == 3:
            t = (n - 3) // 6
            Z = range(2 * t + 1)
            T = lambda x: x[0] + (2 * t + 1) * x[1]

            sts = [[(i, 0), (i, 1), (i, 2)] for i in Z] + [
                [(i, k), (j, k), (((t + 1) * (i + j)) % (2 * t + 1), (k + 1) % 3)] for k in range(3) for i in Z for j in
                Z if i != j]

        elif n % 6 == 1:
            t = (n - 1) // 6
            N = range(2 * t)
            #T = lambda x, y: x + y * t * 2 if (x, y) != (-1, -1) else n - 1
            T = lambda x: x[0] + x[1] * t * 2 if x != (-1, -1) else n - 1

            L1 = lambda i, j: (i + j) % (int((n - 1) / 3))
            L = lambda i, j: L1(i, j) // 2 if L1(i, j) % 2 == 0 else t + (L1(i, j) - 1) // 2

            sts = [[(i, 0), (i, 1), (i, 2)] for i in range(t)] + \
                  [[(-1, -1), (i, k), (i - t, (k + 1) % 3)] for i in range(t, 2 * t) for k in [0, 1, 2]] + \
                  [[(i, k), (j, k), (L(i, j), (k + 1) % 3)] for k in [0, 1, 2] for i in N for j in N if i < j]
        else:
            raise ValueError("Steiner triple systems only exist for n = 1 mod 6 or n = 3 mod 6")
        #self.T = [T(x) for x in sts]
        self.T = sorted(set(list(map(lambda x: tuple(sorted(set(map(T, x)))), sts))))

    ''' Checks strict equality (same elements)
    '''
    def __equals__(self, S):
        return self.order == S.n and self.T == S.T

    ''' Return a copy of the set of triples
    '''
    def get_triples(self):
        return [t[:] for t in self.T]

    ''' Finds the triple having the pair (a, b)
    '''
    def get_triplet_with_pair(self, a, b):
        return list(filter(lambda t: a in t and b in t, self.T))[0]

    ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    '''''''''''''''''''''''''''' new cycles '''''''''''''''''''''''''''
    ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    '''  Creates the cycle graph with pair (x, y)    
    '''
    # Get cycle graph Gab
    def cycle_graph(self, x, y):
        T1 = self.get_triples()

        # All triples having x or y but not both
        Txy = list(filter(lambda t: (x in t or y in t), T1))
        for txy in Txy:
            if x in txy and y in txy:
                Txy.remove(txy)
                break

        # Find connected components
        Vxy = sorted(list(self.X))                  # Vertex set
        for i in txy: Vxy.remove(i)
        d = dict([i, []] for i in Vxy)

        for t in Txy:
            t = list(t)
            if x in t: t.remove(x)
            else: t.remove(y)
            d[t[0]].append(t[1])
            d[t[1]].append(t[0])

        C = []
        while len(d) > 0:
            i, c, v = 1, [Vxy[0], d[Vxy[0]][1]], d[Vxy[0]][1]
            del d[Vxy[0]]
            while v in d.keys():
                a, b = d[c[i]]
                if a in c: c.append(b)
                else: c.append(a)
                del d[c[i]]
                i += 1
                v = c[-1]
            c = c[:-1]
            for v in c: Vxy.remove(v)
            C.append(c)
        return sorted(C, key=len, reverse=1)

    '''  Creates the cycle switching graph with pair (a, b)    
    '''
    # Get cycle graph Gab
    def switch(self, a, b):
        Gab = []
        T1 = self.get_triples()
        
        # All triples having 'a' or 'b' but not both
        Tab = list(filter(lambda t: (a in t or b in t) and not (a in t and b in t), T1))
        # Continue until no more triples remain to be added
        while Tab != []:
            C = []
            s = a
            ab = a+b
            nxt_t = next(t for t in Tab if (a in t))
            p = next(x for x in nxt_t if (x != a))
            P = []
            # Verify if cycle is complete
            while p not in P:
                P.append(p)
                # Append triple and pivot
                C.append([nxt_t, p])
                # Toggle elements a and b
                s = ab-s
                nxt_t = next(t for t in Tab if (p in t and s in t))
                p = next(x for x in nxt_t if (x != s and x != p))
                Tab.remove(nxt_t) 
            Gab.append(CycleGraphComponent(a, b, P, C))
        return Gab
    
    ''' Get all cycle switchs
    '''
    def all_cycle_switchs(self):
        G = {}
        for i in range(self.order-1):
            for j in range(i+1, self.order):
                G[(i,j)] = self.cycle_switch(i,j)
        return G
    
    '''
        Applies permutation phi of elements to produce an isomorphic STS
        parameters
            phi: Permutation to be applied
        returns
            Isomorphic STS
    ''' 
    def permute(self, phi):
        T1 = self.get_triples()
        perm_sts = map(lambda t: sorted(map(lambda x: phi[x], t)), T1)
        #perm_sts = []
        #for t in T1:
        #    perm_sts.append((phi[t[0]], phi[t[1]], phi[t[2]]))
        perm_sts = sorted(map(tuple, perm_sts))
        return SteinerTripleSystem(self.order, perm_sts)

    '''
        Applies a random permutation of elements to produce an isomorphic STS
        returns
            1. Random permutation phi
            2. Isomorphic STS
    '''
    def random_permute(self):
        phi = np.random.permutation(self.X)
        return phi, self.permute(phi)

    '''
        Applies the cycle switch operation to produce another (possibly non-isomorphic) STS.
        If first argument is provide, then the second is ignored 
        input
            C: Component to perform the cycle switch
            t: triple (a,b,p) for finding Cycle switch Gab an component with pivot p or
        returns
            2. New STS resulting from the Switch cycle operation
    '''
    def cycle_switch(self, Cabp=None, t = None):
        if Cabp is None:
            if t is None:
                raise ValueError('Cycle component or pair and pivot required')
            elif t in self.T:
                raise ValueError('Pivot and pair must not be triple of the system')
            elif np.unique(t).size < 3:
                raise ValueError('Elements in triple must be different')
            else:
                a, b, p = t
                Gab = self.switch(a, b)
                for Cabp in Gab:
                    if p in Cabp.pivots:
                        break
        elif not type(Cabp) == CycleGraphComponent:
            raise ValueError('C must be a CycleSwitchComponent')
        C = Cabp.C[:]
        #print('C: ' + str(C))
        a, b = Cabp.pair[:]
        T1 = self.get_triples()
        T2 = [t1[:] for t1 in T1]
        for t in [c[0] for c in C]:
            T2.remove(t)
            t = list(t)
            toadd = []
            if a in t:
                t.remove(a)
                toadd.append(b)
            elif b in t:
                t.remove(b)
                toadd.append(a)
            t.extend(toadd)
            t.sort()
            T2.append(t)
        T2 = sorted([sorted(t) for t in T2])
        return SteinerTripleSystem(self.order, T2)
    
    '''
        Applies the cycle switch operation to a random pair
        inputs
            noisom: If true, a not isomorphic operation is required the resulting STS is non-isomorphic
                    to this. Otherwise the cycle switch could create an isomorphic STS (or not)
    '''
    def random_cycle_switch(self, notisom = False):
        if self.order < 10 and notisom:
            raise ValueError('A not-isomorphic cycle switch cannot be performed with this size')
        while True:
            # Random cycle graph
            a, b = random.sample(self.X, 2)
            Gab = self.switch(a, b)
            # If result is required to be non-isomorphic
            if notisom and len(Gab) == 1:
                continue
            # Select a random component
            i = random.choice(range(len(Gab)))
            S2 = self.cycle_switch(Gab[i])
            return (a, b, Gab[i].pivots[0]), S2

    def cycle_structure(self):
        return


    def __eq__(self, other):
        return set(map(frozenset, self.T)) == set(map(frozenset, other.T))


'''n = 13
S1 = SteinerTripleSystem(n)
print(S1.cycle_graph(0, 2)[0])
S2 = S1.cycle_switch(t = (0,2,1))
print(S1.T)
print(S2.T)
for i in range(n):
    for j in range(i+1, n):
        G1ij = S1.cycle_graph(i, j)
        G2ij = S2.cycle_graph(i, j)

#        G1ij = np.asarray(G1ij[0])
#        G2ij = np.roll(G2ij[0], -1)
        if not np.array_equal(G1ij, G2ij):
            print('{}: {}'.format((i, j), G1ij))
            print('{}: {}'.format((i, j), G2ij))
'''



