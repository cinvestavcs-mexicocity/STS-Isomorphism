from sts_isomorphism import *
import copy
TEST1 = 0
TEST2 = 0


from SteinerTripleSystem import *
d = dict(zip('0123456789abcde', range(15)))
n = 15
'''S = SteinerTripleSystem(n)
S1 = S.cycle_switch(t=(0,1,7))
for i in range(n-1):
    for j in range(i+1, n):
        GC = S1.cycle_graph(i,j)
        GS = S1.switch(i, j)
        for k in range (len(GS)):
            S2 = S1.cycle_switch(GS[k])
            if S2 == S:
                print(i, j, min(GC[k]))'''


filepath = '/Users/edgar/Dropbox/sts{}list.txt'.format(n)
with open(filepath) as fp:
    line = fp.readline()
    S = {}
    while line:
        k, t0, t1, t2 = line.split(' ')[:4]
        k = int(k[:-1])
        T = [(d[t0[i]], d[t1[i]], d[t2[i]]) for i in range(len(t0))]
        S[k] = SteinerTripleSystem(n, T)
        line = fp.readline()

def cycle_list(S):
    Cij = {}
    n = S.order
    for i in range(n):
        for j in range(i+1, n):
            Gij = S.cycle_graph(i, j)
            key = tuple(sorted(map(len, Gij)))
            if key in Cij.keys():
                Cij[key] += 1
            else:
                Cij[key] = 1
    return Cij


# canonical
def canonical_sts(S):
    n = S.order
    phi = list(range(n))
    T0 = sorted(S.T)[:(n+1) // 2]
    for i in range(1, S.order, 2):
        phi[T0[i // 2][1]] = i
        phi[T0[i // 2][2]] = i+1
    return S.permute(phi)

def transform_cycles(S1, t):
    n = S1.order
    x, y, v = t
    C1 = {}
    for i in range(n):
        for j in range(i + 1, n):
            C1[(i, j)] = [G.pivots for G in S1.switch(i, j)]

    S2 = S1.cycle_switch(t=[x, y, v])
    C2 = {}
    for i in range(15):
        for j in range(i + 1, 15):
            C2[(i, j)] = [G.pivots for G in S2.switch(i, j)]

    i = 0
    O = {}
    keys = sorted(C1.keys())
    for k in keys:
        if C1[k] == C2[k]:
            O[k] = C1[k]
            del C1[k]
            del C2[k]

    for k in sorted(C1.keys()):
        print('{}: {} -- {}'.format(k, C1[k], C2[k]))
    print ('\n***********\n')
    for k in sorted(O.keys()):
        print('{}: {}'.format(k, O[k]))

    print(len(O.keys()) + len(C1.keys()))


for i in sorted(S.keys()):
    print('{}: {}'.format(i, cycle_list(S[i])))

# Check how sw affects cycles
if TEST1:
    #for i in range(15):
    #    for j in range(i+1, 15):
    #        print((i,j), S[3].cycle_graph(i, j))
    transform_cycles(S[1], (0, 1, 3))
    S2 = S[1].cycle_switch(t = (0, 1, 3))
    transform_cycles(S2, (0, 3, 10))


# Find isomorphism that changes only last elements
if TEST2:
    S1 = S[1]
    S2 = S1.cycle_switch(t=[5, 6, 9])
    phi = try_isomorphism(S2, S[2], [0,1,2,3,4,5,6,7,-1,-1,-1,-1,-1,-1,-1])
    print(phi)
    print(sts_isomorphism(S2, S[2]))

# Find commutative switches
TEST3 = False
if TEST3:
    S1 = S[1]
    for a in range(S1.order):
        for b in range(a+1, S1.order):
            tab = sorted(S1.get_triplet_with_pair(a, b))
            #  For each possible switch
            for v in S1.X:
                if len(np.unique((a, b, v))) == 3 and tab != sorted((a, b, v)):
                    # check every other sw
                    for x in range(S1.order):
                        for y in range(x + 1, S1.order):
                            txy = sorted(S1.get_triplet_with_pair(x, y))
                            for w in S1.X:
                                if len(np.unique((x, y, w))) == 3 and txy != sorted((x, y, w)):
                                    try:
                                        S2 = S1.cycle_switch(t = (a, b, v))
                                        S2 = S2.cycle_switch(t = (x, y, w))
                                        S3 = S1.cycle_switch(t = (x, y, w))
                                        S3 = S3.cycle_switch(t = (a, b, v))
                                        if S2 == S3:
                                            print((a, b, v), (x, y, w))
                                    except:
                                        print('Error:')






#cl = [cycle_list(S[1])]
#for i in range(1, len(S)):
    #cl.append(cycle_list(S[i]))
    #print('S_{}: {}'.format(i, cycle_list(S[i])))
    #if (12,) not in cl[i].keys():
    #    print('S_{}: {}'.format(i, cl[i]))
    #for j in range(i-1):
    #   if np.array_equal(cl[i], cl[j]):
    #        print(j, i, cl[j])
'''

print(S[6].cycle_graph(0,1))
for i in range(1, 80):
    print(i)
    for j in range(i-1):
        if sts_isomorphism(S[i], S[j]):
            print('{}: {}'.format(j, S[j]))
            print('{}: {}'.format(i, S[i])) 

S1 = SteinerTripleSystem(15)
iso_class, cycle_lists = [S1], [cycle_list(S1)]
i, idx = 0, 0
while idx < 80:
    print(idx)
    S1 = iso_class[idx]
    for i1 in range(14):
        for i2 in range(i1+1, 15):
            cs = S1.switch(i1, i2)
            if len(cs) > 1:
                for C in cs:
                    S2 = S1.cycle_switch(C)
                    for Si in iso_class:
                        found = False
                        if sts_isomorphism(Si, S2):
                            found = True
                            break
                    if not found:
                        print('****{}****'.format(i+1))
                        iso_class.append(S2)
                        cycle_lists.append(cycle_list(S2))
                        first = True
                        for j in range(len(iso_class)-1):
                            if np.array_equal(cycle_lists[j], cycle_lists[-1]):
                                if first:
                                    print(cycle_lists[-1])
                                    first = False
                                print(cycle_lists[j])
                        i += 1
    idx += 1





'''


