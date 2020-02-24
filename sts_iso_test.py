from sts_isomorphism import *
from sts_utils import *

def zkp_init():
    S0 = algebraic_hts(4)
    t, S1 = S0.random_cycle_switch(notisom=True)
    return S0, S1, t


def zkp_round(S0, S1):
    phi, T0 = S0.random_permute()
    T1 = S1.permute(phi)
    return T0, T1, phi

# Key generation
S0, S1, t = zkp_init()
print('t: ', t)

# challenge step
T0, T1, s = zkp_round(S0, S1)
st = [s[x] for x in t]
print('st: ', st)

# permutation equals cycle switch
print('sigma(S1) == sw(T, sigma(t))? : {}'.format(S1.permute(s) == T0.cycle_switch(t = st)))


print('cycle vector S0: ', cycle_vector(S0))
print('cycle vector S1: ', cycle_vector(S1))
print('cycle vector T0: ', cycle_vector(T0))
print('cycle vector T1: ', cycle_vector(T1))

print('isomorphism: {}'.format(sts_improved_miller(S1, T1, 10)))
print('isomorphism: {}'.format(sts_miller(S1, T1)))


'''
# A STS
S1 = SteinerTripleSystem(15)
S2 = S1.cycle_switch(t = (0,5,6))
print('S1: {}'.format(S1.T))
print('G05: {}'.format(S1.cycle_graph(0, 5)))
print('S2: {}'.format(S2.T))
print('G05: {}'.format(S2.cycle_graph(0, 5)))

print('Are isomorphic?')
print(cycle_vector(S1))
print(cycle_vector(S2))'''