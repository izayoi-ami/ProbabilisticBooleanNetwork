from sage.structure.sage_object import SageObject
import operator

class BooleanFunction(SageObject):
    def __init__(self, domain, values):
        self.domain = domain
        sorted(self.domain)
        self.K = len(domain)
        self.N = 2^(self.K)
        if isinstance(values, list):
            self.values = map(GF(2), values) + [GF(2)(0)] * (self.N - len(values))
        else:
            self.set_value_function(values)

    def set_value_function(self, f):
        self.values = map(GF(2),map(f, range(self.N)))

    def set_state_value(self, x, val):
        self.values[self.__to_internal_state(x)] = val

    def evalute(self, k):
        lk = self.translate_global_state_to_local_state(k)
        return self.values[lk]

    def local_evaluate(self, k):
        return self.values[self.__to_internal_state(k)]

    def translate_global_state_to_local_list(self, x):
        xs = self.__to_internal_list(x)
        return [xs[k] if k<len(xs) else 0 for k in self.domain]

    def translate_global_state_to_local_state(self, x):
        return self.__to_internal_state(self.translate_global_state_to_local_list(x))

    def __to_internal_state(self, xs):
        if isinstance(xs, list):
            return ZZ(xs,2)
        return xs

    def __to_internal_list(self, k):
        if isinstance(k, list):
            return k
        return ZZ(k).digits(2, padto=self.dom_size())
    
    def dom_size(self):
        return len(self.domain)

    def combine(self, other, op):
        domain = uniq( self.domain + other.domain )
        sorted(domain)
        xis = [domain.index(x) for x in self.domain]
        yis = [domain.index(y) for y in other.domain]
        values = []
        K = len(domain)
        N = 2^K
        for k in range(N):
            zs = ZZ(k).digits(base=2, padto=K)
            xs = [zs[x] for x in xis]
            ys = [zs[y] for y in yis]
            t = op(self.local_evaluate(xs), other.local_evaluate(ys))
            values.append(t)
        return BooleanFunction(domain, values)

    def __add__(self,other):
        return self.combine(other, lambda x,y: x or y)

    def __mul__(self,other):
        return self.combine(other, lambda x,y: x and y)

    def __repr__(self):
        return "{}-BF Inputs:{}\tValues:{}".format(self.K, self.domain, self.values)

    def __call__(self,xs):
        return self.evalute(xs)

    def __getitem__(self,x):
        return self.local_evaluate(x)

    def __setitem__(self, x, val):
        self.set_state_value(x,val)

    @staticmethod
    def random_element(dom):
        N = 2^len(dom)
        vs = [GF(2).random_element() for _ in range(N)]
        return BooleanFunction(dom, vs)
