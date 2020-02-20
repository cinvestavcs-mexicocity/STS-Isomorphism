from sts_isomorphism import *
from sts_utils import *
import time

MILLER_TEST = 0
SWITCH_TEST = 0
TRADES_TEST = 0
LOOKING_NONISOMORPHIC = 0

#np.random.seed(0)

if MILLER_TEST:
    S1 = SteinerTripleSystem(15)
    print(sorted(S1.T))
    # Permutation and isomorphic system
    f1, S2 = S1.random_permute()
    # Permutation found by isomorphism test (f1 and f2 are not equal in general)
    f2 = sts_miller(S1, S2)
    # Testing equality of triple systems after applying both isomorphisms
    print(S1 == S2)
    print(S1.permute(f1) == S2)
    print(S1.permute(f2) == S2)

    # Difficult instances for the isomorphism problem
    S1 = pg2k_sts(4)
    f1, S2 = S1.random_permute()
    f2 = sts_miller(S1, S2)
    print(S1 == S2)
    print(S1.permute(f2) == S2)
    print(S1.permute(f2) == S2)

if SWITCH_TEST:
    #np.random.seed(None)
    fail = 0
    miller_mean1 = 0
    improved_mean1 = 0
    miller_mean2 = 0
    improved_mean2 = 0
    # Get a big system with small cycles
    for k in range(4, 12):
        print('************{}**************'.format(k))
        S1 = pg2k_sts(k)

        # Checking isomorphisms
        for i in range(10):
            S2 = S1.random_permute()[1]
            start = time.time()
            f = sts_miller(S1, S2)
            miller_mean1 += time.time() - start
            print('Miller TIME: {} --- {}'.format(time.time() - start, S1.permute(f) == S2))
            if not f:
                print('FAIL Miller')
                f = sts_miller(S1, S2)

            start = time.time()
            f = sts_improved_miller(S1, S2, cycles_iters=20)
            improved_mean1 += time.time() - start
            print('Improved TIME: {} --- {}'.format(time.time() - start, S1.permute(f) == S2))
            if not f:
                print('FAIL Improv')

            # Choose a random triple not contained in the block set
            x, y, z = np.random.choice(S1.X, 3, replace=0)
            while z == operate(S1.T, (x, y)):
                x, y, z = np.random.choice(S1.X, 3, replace=0)
            S3 = S1.cycle_switch(t = (x,y,z))
            #cl = cycle_vector(S2)
            # Check difficulty
            phi, S4 = S3.random_permute()
            '''start = time.time()
            f = sts_miller(S3, S4)
            miller_mean2 += time.time() - start
            if not f or not S3.permute(f) == S4:
                print('FAIL Miller')
                #f = sts_miller(S3, S4)
                ss = S3.permute(f)
                result = ss == S4
            else:
                print('Miller TIME: {} --- {}'.format(time.time() - start, S3.permute(f) == S4))'''
            start = time.time()
            f = sts_improved_miller(S3, S4, cycles_iters=20)
            improved_mean2 += time.time() - start
            if not f or not S3.permute(f) == S4:
                print('FAIL Improv')
                #f = sts_improved_miller(S3, S4, cycles_iters=20)
            else:
                print('Improved TIME: {} --- {}'.format(time.time() - start, S3.permute(f) == S4))
        print('***Mean Miller TIME PG: {} --- {}:'.format(k, miller_mean1 / 10))
        print('***Mean Improv TIME PG {} --- {}:'.format(k, improved_mean1 / 10))
        print('***Mean Miller TIME SW: {} --- {}:'.format(k, miller_mean2 / 10))
        print('***Mean Improv TIME SW: {} --- {}:'.format(k, improved_mean2 / 10))

'''        if S2.permute(f) != S3:
            print(phi, f)
            f = sts_improved_miller(S2, S3)
            fail += 1
            break
    print('Failures: {}'.format(fail))'''


if TRADES_TEST:
    # Get a big system with small cycles
    np.random.seed(0)
    S = affine_sts(4)
    P1 = generate_subsystem(S.T, S.T[:2])
    start = time.time()
    print(sts_improved_miller(S, S.random_permute()[1], 50))
    print('TIME: {}'.format(time.time() - start))
    fail = 0
    for i in range(100):
        # Perform a trade
        P2 = P1.random_permute()[1]
        S2 = trade(S, P1.T, P2.T)
        #if len(cl) < 3:
            # Check difficulty
        phi, S3 = S2.random_permute()
        start = time.time()
        f = sts_improved_miller(S2, S3, cycles_iters=50)
        print(time.time() - start)
        if f and S2.permute(f) == S3:
            print(f)
        else:
            print(phi, f)
            #f = sts_improved_miller(S2, S3)
            fail += 1
            break
    print('Failures: {}'.format(fail))


if LOOKING_NONISOMORPHIC:
    # Get a big system with small cycles
    #np.random.seed(0)
    S = affine_sts(5)
    P1 = generate_subsystem(S.T, S.T[:2])
    f = open('/home/edgar/Documents/STS_TEST.txt', 'w+')
    # First STS
    while True:
        # Perform a trade
        P2 = P1.random_permute()[1]
        S2 = trade(S, P1.T, P2.T)
        cv = cycle_vector(S2)
        if len(cv) < 3:
            break

    f.write("S1: {}\n{}".format(cv, S.T))
    # Looking non-isomorphic
    while True:
        # Perform a trade
        P2 = P1.random_permute()[1]
        S3 = trade(S, P1.T, P2.T)
        while cycle_vector(S3) != cv:
            P2 = P1.random_permute()[1]
            S3 = trade(S, P1.T, P2.T)
        f = sts_improved_miller(S2, S3)
        #print(f)
        if not f:
            print('found')
            f.write("P2: {}\nS2: {}".format(P2.T, S3.T))


TEST_ISO_JOIN = 0
if TEST_ISO_JOIN:
    S = pg2k_sts(6)
    print('{}:\t{}'.format('-,-,-', cycle_vector(S)))
    S = S.cycle_switch(t=(0,1,3))
    print('{}:\t{}'.format((0,1,3), cycle_vector(S)))
    S = S.cycle_switch(t=(10, 61, 63))
    print('{}:\t{}'.format((0, 1, 7), cycle_vector(S)))

'''    S2 = S
    for i in range(1,len(G01)):
        x, y, v = G01[i-1][0], G01[i][1], 1
        S2 = S2.cycle_switch(t = (x, y, v))
        G2_01 = S2.cycle_graph(0, 1)
        print('{}: {}'.format(len(G2_01), G2_01))'''

AFF_NAFF_TEST = 1
if AFF_NAFF_TEST:
    #np.random.seed(1)
    for k in range(7, 12):
        start = time.time()
        S1 = affine_sts(k)
        print('Affine Const {}: {}'.format(k, time.time() - start))
        start = time.time()
        _, S1 = S1.random_permute()
        print('Permute: {}: {}'.format(k, time.time() - start))

        start = time.time()
        T1 = algebraic_hts(k)
        print('Algeraic Const {}: {} --- {}'.format(k, time.time() - start, find_quasigroup_generators(quasigroup_from_sts(T1))))
        _, T1 = T1.random_permute()
        #_,
        miller1 = 0
        miller2 = 0
        imp1 = 0
        imp2 = 0
        for iter in range(10):
            # Affine triple system
            _, S2 = S1.random_permute()
            # Miller
            start = time.time()
            phi = sts_miller(S1, S2)
            miller1 += time.time() - start
            print('Miller - Affine {}: {} --- {}'.format(k, time.time() - start, S1.permute(phi) == S2))
            if not S1.permute(phi) == S2:
                print('FAILLLLL')
                break
            # Improved Miller
            start = time.time()
            phi = sts_improved_miller(S1, S2, 1)
            imp1 += time.time() - start
            print('Improved - Affine {}: {} --- {}'.format(k, time.time() - start, S1.permute(phi) == S2))

            # Non-affine
            _, T2 = T1.random_permute()
            # Miller
            '''start = time.time()
            psi = sts_miller(T1, T2)
            miller2 += time.time() - start
            print('Miller - NonAffine {}: {} --- {}'.format(k, time.time() - start, T1.permute(psi) == T2))'''
            # Improved Miller
            '''start = time.time()
            psi = sts_improved_miller(T1, T2, 1)
            imp2 += time.time() - start
            print('Improved - NonAffine {}: {} --- {}'.format(k, time.time() - start, T1.permute(psi) == T2))'''

        print('***********{}: {} --- {} ||| {} --- {}'.format(k, miller1/10, miller2/10, imp1/10, imp2/10))


'''
    from lookfor import *

    np.random.seed(0)
    base = SteinerTripleSystem(3**2).T
    T = SteinerTripleSystem(3**2).T
    T2, T3 = [], []
    for t in T:
        t = tuple([x + 9 for x in t])
        T2.append(t)
        t = tuple([x + 9 for x in t])
        T3.append(t)
    for t in T2:
        base.append(t)
    for t in T3:
        base.append(t)
    S = hill_climb(3**3, base)

    #S = SteinerTripleSystem(3**3, [(0, 1, 5), (0, 2, 4), (0, 3, 6), (0, 7, 8), (0, 9, 24), (0, 10, 26), (0, 11, 25), (0, 12, 20), (0, 13, 18), (0, 14, 23), (0, 15, 19), (0, 16, 21), (0, 17, 22), (1, 2, 3), (1, 4, 7), (1, 6, 8), (1, 9, 20), (1, 10, 25), (1, 11, 19), (1, 12, 21), (1, 13, 22), (1, 14, 18), (1, 15, 24), (1, 16, 23), (1, 17, 26), (2, 5, 8), (2, 6, 7), (2, 9, 19), (2, 10, 23), (2, 11, 26), (2, 12, 18), (2, 13, 24), (2, 14, 22), (2, 15, 21), (2, 16, 20), (2, 17, 25), (3, 4, 8), (3, 5, 7), (3, 9, 23), (3, 10, 24), (3, 11, 22), (3, 12, 19), (3, 13, 21), (3, 14, 25), (3, 15, 18), (3, 16, 26), (3, 17, 20), (4, 5, 6), (4, 10, 21), (4, 11, 24), (4, 12, 23), (4, 13, 20), (4, 14, 19), (4, 15, 26), (4, 16, 22), (4, 17, 18), (5, 9, 18), (5, 10, 19), (5, 11, 23), (5, 12, 22), (5, 13, 25), (5, 14, 26), (5, 15, 20), (5, 16, 24), (5, 17, 21), (6, 9, 22), (6, 10, 20), (6, 11, 18), (6, 12, 24), (6, 13, 26), (6, 14, 21), (6, 15, 23), (6, 16, 25), (6, 17, 19), (7, 9, 26), (7, 10, 18), (7, 11, 21), (7, 12, 25), (7, 13, 23), (7, 14, 20), (7, 15, 22), (7, 16, 19), (7, 17, 24), (8, 9, 21), (8, 10, 22), (8, 11, 20), (8, 12, 26), (8, 13, 19), (8, 14, 24), (8, 15, 25), (8, 16, 18), (8, 17, 23), (9, 10, 14), (9, 11, 13), (9, 12, 15), (9, 16, 17), (10, 11, 12), (10, 13, 16), (10, 15, 17), (11, 14, 17), (11, 15, 16), (12, 13, 17), (12, 14, 16), (13, 14, 15), (18, 19, 23), (18, 20, 22), (18, 21, 24), (18, 25, 26), (19, 20, 21), (19, 22, 25), (19, 24, 26), (20, 23, 26), (20, 24, 25), (21, 22, 26), (21, 23, 25), (22, 23, 24), (4, 9, 25)])
    print(S.T)
    print(cycle_list(S))

    for i in range(27):
        for j in range(i + 1, 27):
            Gij = S.cycle_graph(i, j)
            print('{}, {}: {}'.format(i, j, sorted(map(len, Gij))))
    #S = SteinerTripleSystem(3**3, [(0, 1, 2), (0, 3, 4), (0, 5, 6), (0, 7, 8), (0, 9, 10), (0, 11, 12), (0, 13, 14), (0, 15, 16), (0, 17, 18), (0, 19, 20), (0, 21, 22), (0, 23, 24), (0, 25, 26), (1, 3, 7), (1, 4, 8), (1, 5, 22), (1, 6, 20), (1, 9, 15), (1, 10, 14), (1, 11, 24), (1, 12, 21), (1, 13, 17), (1, 16, 18), (1, 19, 25), (1, 23, 26), (2, 3, 6), (2, 4, 22), (2, 5, 24), (2, 7, 15), (2, 8, 12), (2, 9, 18), (2, 10, 19), (2, 11, 26), (2, 13, 23), (2, 14, 20), (2, 16, 25), (2, 17, 21), (3, 5, 17), (3, 8, 24), (3, 9, 21), (3, 10, 15), (3, 11, 22), (3, 12, 19), (3, 13, 26), (3, 14, 16), (3, 18, 20), (3, 23, 25), (4, 5, 18), (4, 6, 23), (4, 7, 19), (4, 9, 12), (4, 10, 11), (4, 13, 20), (4, 14, 25), (4, 15, 17), (4, 16, 21), (4, 24, 26), (5, 7, 20), (5, 8, 23), (5, 9, 14), (5, 10, 12), (5, 11, 13), (5, 15, 25), (5, 16, 26), (5, 19, 21), (6, 7, 13), (6, 8, 21), (6, 9, 24), (6, 10, 16), (6, 11, 15), (6, 12, 25), (6, 14, 22), (6, 17, 26), (6, 18, 19), (7, 9, 26), (7, 10, 23), (7, 11, 25), (7, 12, 16), (7, 14, 24), (7, 17, 22), (7, 18, 21), (8, 9, 17), (8, 10, 22), (8, 11, 14), (8, 13, 25), (8, 15, 20), (8, 16, 19), (8, 18, 26), (9, 11, 19), (9, 13, 16), (9, 20, 23), (9, 22, 25), (10, 13, 21), (10, 17, 24), (10, 18, 25), (10, 20, 26), (11, 16, 17), (11, 18, 23), (11, 20, 21), (12, 13, 15), (12, 14, 26), (12, 17, 23), (12, 18, 24), (12, 20, 22), (13, 18, 22), (13, 19, 24), (14, 15, 18), (14, 17, 19), (14, 21, 23), (15, 19, 23), (15, 21, 26), (15, 22, 24), (16, 20, 24), (16, 22, 23), (17, 20, 25), (19, 22, 26), (21, 24, 25)])
    #print(cycle_list(S))

    # np.random.seed(0)
    n = 5
    S = af_sts(3, n)
    T = deepcopy(S.T)
    SS = generated([T[0], T[1]], T)
    # print('SS: {}'.format(SS))
    P = SteinerTripleSystem(3 ** 2, SS)
    # cl = cycle_list(S)
    
    while cycle_list(S2) != cl:
        T = deepcopy(S.T)
        phi, P = P.random_permute()
        # Replace
        for t in SS:
            T.remove(t)
        for t in P.T:
            T.append(t)
        S2 = SteinerTripleSystem(3**n, T)

    # Replace
    for i in range(100):
        phi, P = P.random_permute()
        T = deepcopy(S.T)
        for t in SS:
            T.remove(t)
        for t in P.T:
            T.append(t)
        S2 = SteinerTripleSystem(3 ** n, T)
        print(cycle_list(S2))

    
    for i in range(9, 80):
        print(i)
        cg = S2.cycle_graph(1, i)
        for C in cg:
            if len(C) == 12:
                break
        if len(C) != 12:
            continue
        #print(cg)
        choose = [1, i]
        for j in range(6):
            x = choose[j % 2]
            y = C[j]
            z = C[j+6]
            S3 = S2.cycle_switch(t=(x, y, z))
            cl3 = cycle_list(S3)
            if cl[(6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6)] < cl3[(6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6)] or len(cl3.keys()) <= 3:
                print(x, y, z)
                print(cl3)

    
    # print(cg)
    print(cycle_list(S3))

    cg = S3.cycle_graph(0, 9)
    print(cg)
    S3 = S3.cycle_switch(t=(3, 7, 0))
    S3 = S3.cycle_switch(t=(5, 4, 0))
    cg = S3.cycle_graph(0, 9)
    print(cg)
    print(cycle_list(S3))

    for i in range(3 ** n):
        for j in range(i + 1, 3 ** n):
            Gij = S2.cycle_graph(i, j)
            if tuple(map(len, Gij)) != (6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6):
                print('{}, {}: {}'.format(i, j, sorted(map(len, Gij))))
    # S = SteinerTripleSystem(3**3, [(0, 1, 2), (0, 3, 4), (0, 5, 6), (0, 7, 8), (0, 9, 10), (0, 11, 12), (0, 13, 14), (0, 15, 16), (0, 17, 18), (0, 19, 20), (0, 21, 22), (0, 23, 24), (0, 25, 26), (1, 3, 7), (1, 4, 8), (1, 5, 22), (1, 6, 20), (1, 9, 15), (1, 10, 14), (1, 11, 24), (1, 12, 21), (1, 13, 17), (1, 16, 18), (1, 19, 25), (1, 23, 26), (2, 3, 6), (2, 4, 22), (2, 5, 24), (2, 7, 15), (2, 8, 12), (2, 9, 18), (2, 10, 19), (2, 11, 26), (2, 13, 23), (2, 14, 20), (2, 16, 25), (2, 17, 21), (3, 5, 17), (3, 8, 24), (3, 9, 21), (3, 10, 15), (3, 11, 22), (3, 12, 19), (3, 13, 26), (3, 14, 16), (3, 18, 20), (3, 23, 25), (4, 5, 18), (4, 6, 23), (4, 7, 19), (4, 9, 12), (4, 10, 11), (4, 13, 20), (4, 14, 25), (4, 15, 17), (4, 16, 21), (4, 24, 26), (5, 7, 20), (5, 8, 23), (5, 9, 14), (5, 10, 12), (5, 11, 13), (5, 15, 25), (5, 16, 26), (5, 19, 21), (6, 7, 13), (6, 8, 21), (6, 9, 24), (6, 10, 16), (6, 11, 15), (6, 12, 25), (6, 14, 22), (6, 17, 26), (6, 18, 19), (7, 9, 26), (7, 10, 23), (7, 11, 25), (7, 12, 16), (7, 14, 24), (7, 17, 22), (7, 18, 21), (8, 9, 17), (8, 10, 22), (8, 11, 14), (8, 13, 25), (8, 15, 20), (8, 16, 19), (8, 18, 26), (9, 11, 19), (9, 13, 16), (9, 20, 23), (9, 22, 25), (10, 13, 21), (10, 17, 24), (10, 18, 25), (10, 20, 26), (11, 16, 17), (11, 18, 23), (11, 20, 21), (12, 13, 15), (12, 14, 26), (12, 17, 23), (12, 18, 24), (12, 20, 22), (13, 18, 22), (13, 19, 24), (14, 15, 18), (14, 17, 19), (14, 21, 23), (15, 19, 23), (15, 21, 26), (15, 22, 24), (16, 20, 24), (16, 22, 23), (17, 20, 25), (19, 22, 26), (21, 24, 25)])
    # print(cycle_list(S))
    # print(phi)
    # print(sorted(S2.T))
    # print(sorted(S.T))
    # print(sts_isomorphism(S, S2))'''
