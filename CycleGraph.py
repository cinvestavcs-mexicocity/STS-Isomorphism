''' Cycle Graph
'''
class CycleGraphComponent(object):
    ''' Initializes the cycle graph with the given pivots and element pair (a,b)
    '''

    def __init__(self, a, b, pivots, C):
        self.pair = (a, b)
        self.pivots = pivots[:]
        self.C = C

    def get_pivots(self):
        return self.pivots[:]

    def __len__(self):
        return len(self.pivots)