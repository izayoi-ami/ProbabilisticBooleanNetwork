from sage.structure.sage_object import SageObject
import operator

class BooleanFunctionSet(SageObject):
    def __init__(self, fs, ps=None):
        self.initialize()
        self.fs = fs
        self.set_probability_vector(ps)

    def initialize(self):
        self._local_dis = None 
        self._global_state_active_p = dict()

    def set_probability_vector(self, ps):
        if ps is None:
            ps = [1] * self.size() 
        self.ps = vector(ps).apply_map(abs).normalized(1)
        self.X = GeneralDiscreteDistribution(self.ps)
        self.initialize()

    def append_function(self, f, p=None):
        if p is None:
            p = 1/self.size()
        self.fs.append(f)
        v = self.ps.list()
        v.append(p)
        self.set_probability_vector(v)

    def random_function(self):
        k = self.X.get_random_element()
        return self[k]

    def global_state_activated_probability(self, x):
        p_true = self._global_state_active_p.get(x)
        if p_true is None:
            p_true = sum([p for f,p in zip(self.fs, self.ps) if f(x)])
            self._global_state_active_p[x] = p_true
        return p_true

    def local_state_activated_probability(self, x):
        if self._local_dis is not None:
            return self._local_dis[x]
        p_true = sum([self.ps[k] for k in range(self.size()) if self.fs[k][self.translate_local_state_to_internal_state(x,k)] ])
        return p_true

    def local_state_distribution(self):
        if self._local_dis is not None:
            return self._local_dis
        K = self.local_state_space()
        self._local_dis = [self.local_state_activated_probability(x) for x in K]
        return self._local_dis

    def translate_local_state_to_internal_state(self,x,k):
        """ translate to the k^th function's internal state
        """
        dom = self.domain()
        fd = [dom.index(var) for var in self.fs[k].domain]
        xs = self.__to_list(x)
        s = self.__to_state([xs[i] for i in fd])
        return s


    def __to_list(self, x):
        return ZZ(x).digits(2, padto=self.dom_size())

    def size(self):
        return len(self.fs)

    def domain(self):
        dom = uniq(sum([f.domain for f in self.fs],[]))
        sorted(dom)
        return dom

    def dom_size(self):
        return len(self.domain())

    def local_state_size(self):
        return 2^self.dom_size()

    def local_state_space(self):
        return range(self.local_state_size())

    def get_function(self, k):
        return self.fs[k]

    def evaluate(self, x, k = None):
        #If k is given, evalutes the k^th function at global state x.
        #Otherwise, evalutes a random function.
        if k is not None:
            return self.fs[k](x)
        return self.random_function()(x)

    def local_evaluate(self, x, k=None):
        return self.fs[k][x]

    def conjunctive_function(self):
        """ Return a BooleanFunction which is the conjunction of all functions in the set
        """
        return prod(self.fs, self.fs[0])

    def disjunctive_function(self):
        """ Return a BooleanFunction which is the disjunction of all functions in the set
        """
        return sum(self.fs, self.fs[0])

    def conjunction(self):
        return BooleanFunctionSet([self.conjunctive_function()])

    def disjunction(self):
        return BooleanFunctionSet([self.disjunctive_function()])

    def inspect(self):
        print "\n".join(["{}\t({})".format(repr(f),p) for p,f in zip(self.ps, self.fs)])

    def full_repr(self):
        return "\n".join(["{}\t({})".format(repr(f),p) for p,f in zip(self.ps, self.fs)])

    def __to_state(self, xs):
        if isinstance(xs, list):
            return ZZ(xs,2)
        return xs

    def __to_list(self, k):
        if isinstance(k, list):
            return k
        return ZZ(k).digits(2, padto=self.dom_size())

    def __getitem__(self, k):
        return self.get_function(k)

    def __call__(self, x):
        return self.evaluate(x)

    def __repr__(self):
        return "Boolean set of {} functions:\t {}".format(self.size(), self.ps)

    @staticmethod
    def random_element(dom, min_variables=1, max_variables=None, min_set_size = 1,  max_set_size=5):
        if max_variables is None: max_variables = len(dom)
        n = Combinations(range(min_set_size, max_set_size+1),1).random_element()
        n = n[0]
        fs = []
        for _ in range(n):
            k = Combinations(range(min_variables, max_variables+1),1).random_element()
            k = k[0]
            fs.append( BooleanFunction.random_element(list(Subsets(dom,k).random_element())) )
        ps = random_vector(RDF,n,0,1).normalized(1)
        return BooleanFunctionSet(fs, ps)

        
