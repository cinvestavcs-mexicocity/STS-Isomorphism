

class SteinerQuasigroup(object):
    def __init__(self, X_set=None, mult_dict=None):
        if X_set == None:
            X_set = []
            mult_dict = {}
        self.X = X_set
        self.order = len(X_set)
        if self.verify_quasi(mult_dict):
            self.op = mult_dict.copy()

    def verify_quasi(self, mult_table):
        # Trust by now. Only verify if it is latin square
        return True

    def operate(self, x, y):
        if x == y:
            return x
        elif x > y:
            return self.op[(y, x)]
        else:
            return self.op[(x, y)]

    def generate(self, Y):
        H = SteinerQuasigroup()
        for y in Y:
            H = self.extend(H, y)
        return H

    def extend(self, H, y):
        X = set(H.X)
        mtable = H.op
        if y in X:
            return H
        new_elements = set([y])
        while len(new_elements) > 0:
            y = new_elements.pop()
            for x in X:
                z = self.operate(x, y)
                if z not in X:
                    new_elements.add(z)
                    if x > y:
                        mtable[(y, x)] = z
                    elif x < y:
                        mtable[(x, y)] = z
            X.add(y)
        return SteinerQuasigroup(X, mtable)
