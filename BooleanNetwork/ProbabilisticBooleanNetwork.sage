from sage.structure.sage_object import SageObject
import operator

class ProbabilisticBooleanNetwork(SageObject):
    """ Network (DiGraph) with V vertices, each having one BooleanFunctionSet
        Number of System State is N
    """

    def __init__(self, fs, name={},  state=0, sync=True, fname="", perturbation = "one", param = None):
        """ A list of BooleanFunctionSet
        """
        import tempfile
        self.__initialize()
        self.fs = fs
        self.name = name
        if fname == "":
            _, fname = tempfile.mkstemp(suffix = ".sobj", dir="./data/tmp/")
        self.fname = fname
        self.state = state
        self.set_sync_mode(sync)
        self.set_perturbation(mode = perturbation, param=param)

    def __initialize(self):
        self.name = {}
        self.fs = []
        self.sync = True
        self._network_graph = None

        self._perturbation_mode = None
        self._one_bit_perturbation_matrix = None
        self._full_perturbation_matrix = None

        self.state = 0

        self._transition_graph = None
        self._state_out_neighbourhood = ProbabilisticBooleanNetwork.ldict()
        self._state_probability = ProbabilisticBooleanNetwork.mdict()
        self._transition_dict = None
        self._transition_pmatrix = None
        self._sync_computed = False
        self._sync_scc = None
        self._sync_bscc = None
        self._sync_basin = None
        self._sync_plot_order = None
        self._sync_bscc_plot_order = None

        self._async_transition_graph = None
        self._async_state_out_neighbourhood = ProbabilisticBooleanNetwork.ldict()
        self._async_state_probability = ProbabilisticBooleanNetwork.mdict()
        self._async_transition_dict = None
        self._async_transition_pmatrix = None
        self._async_computed = False
        self._async_scc = None
        self._async_bscc = None
        self._async_basin = None
        self._async_plot_order = None
        self._async_bscc_plot_order = None
        self.__clear_perturbation_state()

    def set_sync_mode(self, is_sync):
        self.sync = is_sync
    
    def set_perturbation(self, param = None, mode="one"):
        """ Set parameters for perturbation
            "one":  One bit perturbation mode
                    Param: [p,q] 
                    p: joint probability for one bit perturbation,
                    q: staying probability
            "full": Full perturbation, V identically independent Bernoulli Distributions
                    Param: p , the parameter for Bernoulli Distributions
        """
        if mode == "one": return self.set_one_perturbation(param)
        if mode == "full": return self.set_full_perturbation(param)
        if mode == "off": return self.set_one_perturbation([0,0])
        raise ValueError("Unknown perturbation mode")

    def V(self):
        return len(self.fs)

    def N(self):
        return 2^self.V()

    def node_in_neighbourhood(self, k):
        return self.fs[k].domain()

    def node_indegree(self, k):
        return self.fs[k].size()

    def state_space(self):
       return range(self.N())

    def nodes(self):
        return range(self.V())
    
    def steady_state(self):
        return self.iterative_steady_state()

    def iterative_steady_state(self, state=None, eps=1e-7, max_iteration=10^7, power = None, sync=None):
        if power is None:
            mode = self.perturbation_mode()
            if mode == "one": power, eps = 1, 1e-10
            if mode == "full": power = 200
        N = self.N()
        a = random_vector(RDF, N,0,1).normalized(1)
        P = self.probability_matrix(sync = sync)
        b = a*P
        for _ in range(max_iteration):
            P = P^power
            a, b = b, b*P
            if norm(a-b) < eps:
                break
        return b

    def sccs(self):
        return self.strongly_connected_components()

    def bsccs(self):
        return self.bottom_strongly_connected_components()
    
    def strongly_connected_components(self):
        if self.sync: return self.sync_strongly_connected_components()
        return self.async_strongly_connected_components()

    def bottom_strongly_connected_components(self):
        if self.sync: return self.sync_bottom_strongly_connected_components()
        return self.async_bottom_strongly_connected_components()

    def basins(self):
        if self.sync: return self.sync_basins()
        return self.async_basins()
    
    def basin_plot(self, f):
        if self.sync: return self.sync_basin_plot(f)
        return self.async_basin_plot(f)

    def basin_data_plot(self, L):
        return self.basin_plot(lambda x:L[x])

    def bscc_plot(self, f):
        if self.sync: return self.sync_bscc_plot(f)
        return self.async_bscc_plot(f)

    def bscc_data_plot(self, L):
        return self.bscc_plot(lambda x:L[x])

    def sampling_basin_plot(self, state=None, length=None, rel=0.1):
        return self.basin_data_plot(self.sampling_probability_vector(state, length, rel))

    def sampling_bscc_plot(self, state=None, length=None, rel=0.1):
        return self.bscc_data_plot(self.sampling_probability_vector(state, length, rel))

    def sampling_plot(self, state=None, length=None, rel=0.1, f=None, pos=None, sep=None):
        pv = map(f, self.sampling_probability_vector(state, length, rel))
        return self.data_plot(pv, pos, sep)

    def probability_matrix(self, sync=None):
        if sync is None:
            sync = self.sync
        if sync: return self.sync_probability_matrix()
        return self.async_probability_matrix()

    def perturbation_mode(self):
        if self._perturbation_mode in ["one", "full", "off"]:
            return self._perturbation_mode
        if self._perturbation_mode is None:
            self.set_perturbation()
        raise ValueError("Unknown perturbation mode")

    def network_graph(self):
        if self._network_graph is None:
            es = {k:self.node_in_neighbourhood(k) for k in self.nodes()}
            D = DiGraph(es)
            self._network_graph = D.reverse()
        return self._network_graph

    def next_state(self, state=None):
        if state is None: state = self.state
        if self.sync: return self.sync_next_state(state)
        return self.async_next_state(state)

    def next(self):
        self.state = self(self.state)
        return self.state

    def plot(self, f, pos=None, sep=None):
        width = 1
        fpx = map(f, map(pos, self.state_space()))
        ymin, ymax = min(fpx), max(fpx)
        ymin = ymin - 0.01*ymin
        ymax = ymax + 0.01*ymax
        Ls=[]
        bc = bar_chart(fpx, width=width)
        if sep is not None:
            Ls=[line([(k*width,ymin),(k*width,ymax)], linestyle="--", color="green") for k in sep if k < self.N()]
            bc = bc + sum(Ls) 
        return bc

    def data_plot(self, L, f=None, pos=None, sep=None):
        pv = map(f, L)
        return self.plot(lambda x:pv[x], pos, sep)

    def set_state(self, state):
        self.state = self.__to_state(state)

    def generate_trajectory(self, state = None, length=None, rel = 0.1):
        """ Generate trajectory of given (relative) length
            state: Initial state (default: system state)
            length: Absolute length (default: None)
            rel: Length relative to state space (default:0.1)
        """
        if length is None: length = ceil(rel * self.N())
        if state is None: state = self.state
        xs = []
        for _ in range(length):
            xs.append(state)
            state = self(state)
        return xs

    def sampling_probability_vector(self, state=None, length= None, rel=0.1):
        from scipy.stats import relfreq
        xs = self.generate_trajectory(state, length, rel)
        N = self.N()
        rel = relfreq(xs, N, (0,N))
        return vector(rel[0])

    def perturbated_next_state(self, state=None):
        ns = None
        if state is None: state = self.state
        perturbation = {
            "one" : self.one_perturbation_next_state,
            "full" : self.full_perturbation_next_state
        }
        ns = perturbation[self._perturbation_mode](self._param, state)
        #print "{}:from {} to {}".format(self._perturbation_mode, state, ns)
        return ns

    def full_perturbation_next_state(self, param, state=None):
        if state is None: state = self.state
        p = param
        X = GeneralDiscreteDistribution([1-p, p])
        xs = [X.get_random_element() for _ in self.nodes()]
        x = int(ZZ(xs,2))
        state = int(state)
        if x == 0: return None
        return state ^^ x

    def one_perturbation_next_state(self, param, state=None):
        if state is None: state = self.state
        p,q = param
        X = GeneralDiscreteDistribution([1-(p+q),p,q])
        r = X.get_random_element()
        if r == 0: return None
        if r == 2: return state
        xs = self.one_bit_neighbourhood(state)
        return self.__random_uniform_element(xs)

    def set_one_perturbation(self, param=None):
        p,q = 0.01, 0.01
        if param is not None: p,q = param
        self._perturbation_mode = "one"
        self.__clear_perturbation_state()
        self._param = [p,q]
        self._one_bit_perturbation_matrix = None

    def set_full_perturbation(self, param=None):
        p = 0.001
        if param is not None: p = param
        self._perturbation_mode = "full"
        self.__clear_perturbation_state()
        self._param = p
        self._full_perturbation_matrix = None

    def perturbation_matrix(self):
        if self._perturbation_mode is None: self.set_one_perturbation()
        if self._perturbation_mode == "one": return self.one_bit_perturbation_matrix()
        if self._perturbation_mode == "full": return self.full_perturbation_matrix()
        raise ValueError("Unknown perturbation mode")

    def full_perturbation_matrix(self):
        if self._full_perturbation_matrix is None:
            p = self._param
            self._full_perturbation_matrix = matrix(RDF, self.N(), lambda x,y: self.__full_perturbation_probability(p,x,y))
        return self._full_perturbation_matrix

    def one_bit_perturbation_matrix(self):
        if self._one_bit_perturbation_matrix is None:
            p,q = self._param
            pn = p / self.V()
            md = ProbabilisticBooleanNetwork.ddict()
            for x in self.state_space():
                for y in self.one_bit_neighbourhood(x):
                    md[(x,y)] = pn
                    md[(x,x)] = q
            self._one_bit_perturbation_matrix = matrix(RDF, self.N(), md)
        return self._one_bit_perturbation_matrix


    def sync_probability_matrix(self):
        if self._sync_probability_matrix is not None:
            return self._sync_probability_matrix
        pm = self.perturbation_matrix()
        tm = self.sync_transition_matrix()
        if self._perturbation_mode == "one":
            p,q = self._param
            self._sync_probability_matrix = pm + (1-p-q)*tm
        if self._perturbation_mode == "full":
            p = self._param
            self._sync_probability_matrix = pm + (1-p)^(self.V()) * tm
        return self._sync_probability_matrix
    
    def sync_transition_graph(self):
        if self._transition_graph is None:
            es = {k:self.sync_state_out_neighbourhood(k) for k in self.state_space()}
            D = DiGraph(es)
            self._transition_graph = D
        return self._transition_graph

    def sync_state_transition_probability(self, x, y):
        if not self._sync_computed:
            self.compute_sync_state()
        p = self._state_probability[x].get(y, 0.0)
        return p

    def sync_state_transition_vector(self, x):
        if not self._sync_computed:
            self.compute_sync_state()
        return self._state_probability[x].copy()

    def sync_state_out_neighbourhood(self, x):
        outs = self._state_out_neighbourhood.get(x) 
        if outs is None:
            xss = []
            yss = []
            for f in self.fs:
                p = f.global_state_activated_probability(x)
                xs , ys = self.__possible_bits(p)
                xss.append(xs)
                yss.append(ys)
            xss = cartesian_product(xss)
            yss = cartesian_product(yss)
            outs = [self.__to_state(list(xs)) for xs in xss]
            ps = [prod(list(ys)) for ys in yss]
            for o,p in zip(outs, ps):
                self._state_probability[x][o] += p
            self._state_out_neighbourhood[x].extend(outs)
        return outs
    
    def sync_transition_matrix(self):
        if not self._sync_computed:
            self.compute_sync_state()
        return self._transition_pmatrix

    def sync_next_transition_state(self, state=None):
        if state is None: state = self.state
        xs = [f(state) for f in self.fs]
        return self.__to_state(xs)

    def sync_next_state(self, state=None):
        if state is None: state = self.state
        ns = self.perturbated_next_state(state)
        if ns is None: ns = self.sync_next_transition_state(state)
        return ns

    def compute_sync_state(self):
        for x in self.state_space():
            self.sync_state_out_neighbourhood(x)
        self._sync_computed = True
        self._transition_dict = ProbabilisticBooleanNetwork.ddict()
        for x in self.state_space():
            d = self.sync_state_transition_vector(x)
            for y,p in d.iteritems():
                self._transition_dict[(x,y)] = p
        self._transition_pmatrix = matrix(RDF, self.N(), self._transition_dict)

    def sync_strongly_connected_components(self):
        if self._sync_scc is None:
            self._sync_scc = self.sync_transition_graph().strongly_connected_components()
        return self._sync_scc

    def sync_bottom_strongly_connected_components(self):
        if self._sync_bscc is not None:
            return self._sync_bscc
        G = self.sync_transition_graph()
        sccs = self.sync_strongly_connected_components()
        bsccs = [scc for scc in sccs if self.__is_bscc(G, scc)]
        self._sync_bscc = bsccs
        return bsccs

    def sync_basins(self):
        if self._sync_basin is not None:
            return self._sync_basin
        self._sync_plot_order = []
        k = 0
        bscc = self.sync_bottom_strongly_connected_components()
        G = self.sync_transition_graph()
        basins = []
        for scc in bscc:
            g = G.breadth_first_search(scc, ignore_direction = True)
            L = list(g)
            basins.append(L)
            self._sync_plot_order.extend(L)
        self._sync_basin = basins
        return basins

    def sync_plot_pos(self, x):
        if self._sync_plot_order is None:
            self.sync_basins()
        return self._sync_plot_order[x]

    def sync_bscc_plot_pos(self, x):
        if self._sync_bscc_plot_order is None:
            bsccs = self.sync_bottom_strongly_connected_components()
            empty = [True] * self.N()
            pos = []
            G = self.sync_transition_graph()
            for bscc in bsccs:
                pos.extend(bscc)
                for v in bscc:
                    empty[v] = False
            for bscc in bsccs:
                g = G.breadth_first_search(bscc, ignore_direction = True)
                L = filter(lambda x:empty[x], list(g))
                pos.extend(L)
                for v in L:
                    empty[v] = False
            self._sync_bscc_plot_order = pos
        return self._sync_bscc_plot_order[x]

    def sync_bscc_plot(self, f):
        from numpy import cumsum
        pos = self.sync_bscc_plot_pos
        bscc = self.sync_bottom_strongly_connected_components()
        sep = [len(b) for b in bscc]
        sep = cumsum(sep)
        return self.plot(f, pos, sep)

    def sync_bscc_data_plot(self, L):
        return self.sync_bscc_plot(lambda x:L[x])

    def sync_basin_plot(self, f):
        """ f: StateSpace -> Real
            Return a plot order by basins

            or a list of N data
        """
        from numpy import cumsum
        width = 1
        basins = self.sync_basins()
        fpx = [f(self.sync_plot_pos(x)) for x in self.state_space()]
        pos = [len(b) for b in basins]
        pos = cumsum(pos)
        ymax = max(fpx)
        ymin = min(fpx)
        ymin = ymin - 0.01*ymin
        ymax = ymax + 0.01*ymax
        Ls=[line([(k*width,ymin),(k*width,ymax)], linestyle="--", color="green") for k in pos if k <= self.N()]
        bc = bar_chart(fpx, width=width)+sum(Ls)
        return bc

    def sync_basin_data_plot(self, L):
        return self.sync_basin_plot(lambda x:L[x])

    def async_probability_matrix(self):
        if self._async_probability_matrix is not None:
            return self._async_probability_matrix
        pm = self.perturbation_matrix()
        tm = self.async_transition_matrix()
        if self._perturbation_mode == "one":
            p,q = self._param
            self._async_probability_matrix = pm + (1-p-q)*tm
        if self._perturbation_mode == "full":
            p = self._param
            self._async_probability_matrix = pm + (1-p)^(self.V()) * tm
        return self._async_probability_matrix
    
    def async_transition_graph(self):
        if self._async_transition_graph is None:
            es = {k:self.async_state_out_neighbourhood(k) for k in self.state_space()}
            D = DiGraph(es)
            self._async_transition_graph = D
        return self._async_transition_graph

    def async_state_out_neighbourhood(self, x):
        outs = self._async_state_out_neighbourhood.get(x)
        if outs is None:
            self.compute_async_state()
        return self._async_state_out_neighbourhood.get(x)

    def async_state_transition_probability(self, x, y):
        if not self._async_computed:
            self.compute_async_state()
        p = self._async_state_probability[x].get(y, 0.0)
        return p

    def async_state_transition_vector(self, x):
        if not self._async_computed:
            self.compute_async_state()
        return self._async_state_probability[x].copy()

    def async_transition_matrix(self):
        if not self._async_computed:
            self.compute_async_state()
        return self._async_transition_pmatrix

    def async_next_transition_state(self, state=None):
        if state is None: state = self.state
        x = self.sync_next_transition_state(state)
        xs = self.split_sync_transitions_to_one_step_async(state, x)
        return self.__random_uniform_element(xs)

    def async_next_state(self, state=None):
        if state is None: state = self.state
        ns = self.perturbated_next_state(state)
        if ns is None: ns = self.async_next_transition_state(state)
        return ns

    def compute_async_state(self):
        for x in self.state_space():
            self.__async_state_partial_computation(x)
        self._async_computed = True
        self._async_transition_dict = ProbabilisticBooleanNetwork.ddict()
        for x in self.state_space():
            d = self.async_state_transition_vector(x)
            for y,p in d.iteritems():
                self._async_transition_dict[(x,y)] = p
        self._async_transition_pmatrix = matrix(RDF, self.N(), self._async_transition_dict)

    def async_strongly_connected_components(self):
        if self._async_scc is None:
            self._async_scc = self.async_transition_graph().strongly_connected_components()
        return self._async_scc

    def async_bottom_strongly_connected_components(self):
        if self._async_bscc is not None:
            return self._async_bscc
        G = self.async_transition_graph()
        sccs = self.async_strongly_connected_components()
        bsccs = [scc for scc in sccs if self.__is_bscc(G, scc)]
        self._async_bscc = bsccs
        return bsccs

    def async_basins(self):
        if self._async_basin is not None:
            return self._async_basin
        self._async_plot_order = []
        k = 0
        bscc = self.async_bottom_strongly_connected_components()
        G = self.async_transition_graph()
        basins = []
        for scc in bscc:
            g = G.breadth_first_search(scc, ignore_direction = True)
            L = list(g)
            basins.append(L)
            self._async_plot_order.extend(L)
        self._async_basin = basins
        return basins

    def async_plot_pos(self, x):
        if self._async_plot_order is None:
            self.async_basins()
        return self._async_plot_order[x]

    def async_bscc_plot_pos(self, x):
        if self._async_bscc_plot_order is None:
            bsccs = self.async_bottom_strongly_connected_components()
            empty = [True] * self.N()
            pos = []
            G = self.async_transition_graph()
            for bscc in bsccs:
                pos.extend(bscc)
                for v in bscc:
                    empty[v] = False
            for bscc in bsccs:
                g = G.breadth_first_search(bscc, ignore_direction = True)
                L = filter(lambda x:empty[x], list(g))
                pos.extend(L)
                for v in L:
                    empty[v] = False
            self._async_bscc_plot_order = pos
        return self._async_bscc_plot_order[x]

    def async_bscc_plot(self, f):
        from numpy import cumsum
        pos = self.async_bscc_plot_pos
        bscc = self.async_bottom_strongly_connected_components()
        sep = [len(b) for b in bscc]
        sep = cumsum(sep)
        return self.plot(f, pos, sep)

    def async_bscc_data_plot(self, L):
        return self.async_bscc_plot(lambda x:L[x])

    def async_basin_plot(self, f):
        """ f: StateSpace -> Real
            Return a plot order by basins

            or a list of N data
        """
        from numpy import cumsum
        width = 1
        basins = self.async_basins()
        fpx = [f(self.async_plot_pos(x)) for x in self.state_space()]
        pos = [len(b) for b in basins]
        pos = cumsum(pos)
        ymax = max(fpx)
        ymin = min(fpx)
        ymin = ymin - 0.01*ymin
        ymax = ymax + 0.01*ymax
        Ls=[line([(k*width,ymin),(k*width,ymax)], linestyle="--", color="green") for k in pos if k <= self.N()]
        bc = bar_chart(fpx, width=width)+sum(Ls)
        return bc

    def async_basin_data_plot(self, L):
        return self.async_basin_plot(lambda x:L[x])

    def export(self):
        self.export_all_graphs()
        self.export_all_bscc()


    def export_graph(self, graph, ext, filename=None):
        import subprocess
        if filename is None:
            filename = self.fname
        root = "{}.{}".format(filename, ext)
        dotfile = "{}.dot".format(root)
        graph.graphviz_to_file_named(dotfile)
        subprocess.Popen(["dot", "-Tps", "-O", dotfile])
        #subprocess.Popen(["dot", "-Tpng", "-O", dotfile])

    def export_sync_state_graph(self, filename=None):
        self.export_graph(self.sync_transition_graph(), "sync", filename)

    def export_async_state_graph(self, filename=None):
        self.export_graph(self.async_transition_graph(), "async", filename)

    def export_network_graph(self, filename=None):
        self.export_graph(self.network_graph(), "net", filename)

    def export_all_graphs(self, filename=None):
        self.export_network_graph(filename)
        self.export_sync_state_graph(filename)
        self.export_async_state_graph(filename)

    def export_bscc(self, ext, filename=None):
        if filename is None:
            filename = self.fname
        root = "{}.{}".format(filename, ext)
        bsccfile = "{}.bscc ".format(root)
        bscc = None
        if ext == "sync":
            bscc = self.sync_bottom_strongly_connected_components()
        elif ext == "async":
            bscc = self.async_bottom_strongly_connected_components()
        else: return None

        with open(bsccfile, "w") as f:
            f.write("\n".join(map(str, bscc)))

    def export_sync_bscc(self, filename=None):
        self.export_bscc("sync", filename)

    def export_async_bscc(self, filename=None):
        self.export_bscc("async", filename)

    def export_all_bscc(self, filename=None):
        self.export_sync_bscc(filename)
        self.export_async_bscc(filename)


    def conjunction(self):
        fs = [f.conjunction() for f in self.fs]
        return ProbabilisticBooleanNetwork(fs, name=self.name, fname="{}.con".format(self.fname))

    def disjunction(self):
        fs = [f.disjunction() for f in self.fs]
        return ProbabilisticBooleanNetwork(fs, name=self.name, fname="{}.dis".format(self.fname))

    def inspect(self):
        for k in self.nodes():
            print "Node {}:".format(k)
            self.fs[k].inspect()

    def check_sum(self):
        for k in self.state_space():
            v = self.sync_state_transition_vector(k)
            print "Sync Node {}:\t Sum: {}\n{}".format(k,sum([v[j] for j in v]), v)
        for k in self.state_space():
            v = self.async_state_transition_vector(k)
            print "Async Node {}:\t Sum: {}\n{}".format(k,sum([v[j] for j in v]), v)

    def mode(self):
        if self.sync: return "Synchronous"
        return "Asynchronous"

    def __repr__(self):
        mode = self.perturbation_mode()
        if mode == "one":
            if sum(self._param) == 0: mode = "off"

        head = "{} {}-Node Network\tPerturbation_mode:{}\tparam:{}".format(self.mode(), self.V(), mode, self._param)

        body ="\n".join(["Node {}:\n{}".format(k, self.fs[k].full_repr()) for k in self.nodes()])
        return "\n".join([head,body])

    def __call__(self, state=None, sync=None):
        if sync is None:
            if state is None: return self.next()
            return self.next_state(state)
        if sync: return self.sync_next_state(state)
        else: return self.async_next_state(state)

    def __async_state_partial_computation(self, x):
        ss = self.sync_state_out_neighbourhood(x)
        for s in ss:
            p = self.sync_state_transition_probability(x,s)
            ns = self.split_sync_transitions_to_one_step_async(x, s)
            L = len(ns)
            self._async_state_out_neighbourhood[x].extend(ns)
            for n in ns:
                self._async_state_probability[x][n] += p/L
        self._async_state_out_neighbourhood[x] = uniq(self._async_state_out_neighbourhood[x])

    def __possible_bits(self, p):
        if p == 0 or p == 1:
            return [p], [1]
        return [0,1], [1-p, p]

    def __to_state(self,xs):
        if isinstance(xs, list):
            return ZZ(xs,2)
        return xs

    def __to_list(self, x):
        if isinstance(x, list):
            return x
        return ZZ(x).digits(2, padto=self.nodes())

    def one_bit_neighbourhood(self, x):
        x = int(x)
        return [ x^^(1<<b) for b in self.nodes()]

    def __full_perturbation_probability(self,p, x,y):
        if x == y: return 0
        k = self.hamming_distance(x,y)
        return p^k * (1-p)^(self.V()-k)

    def __clear_perturbation_state(self):
        self._sync_probability_matrix = None
        self._async_probability_matrix = None
        self._param = None

    @staticmethod
    def __is_bscc(G, scc):
        """ Check if a given SCC of a graph G is a BSCC
        """
        L = len(scc)
        g = G.depth_first_search(scc)
        try:
            for _ in range(L+1):
                g.next()
            return False
        except StopIteration:
            return True
        return False

    @staticmethod
    def __basin_generator(self, G, bscc):
        return G.breadth_first_search(bscc, ignore_direction=True)

    @staticmethod
    def __random_uniform_element(L):
        X = GeneralDiscreteDistribution([1]* len(L))
        return L[X.get_random_element()]
    
    @staticmethod
    def split_sync_transitions_to_one_step_async(x, y):
        if x == y: return [x]
        x = int(x)
        y = int(y)
        t = x.__xor__(y)
        ds = ZZ(t).digits(2)
        bs = range(len(ds))
        ns = [ x^^(1<<b) for d, b in zip(ds, bs) if d == 1]
        return ns    

    @staticmethod
    def ldict():
        from collections import defaultdict
        return defaultdict(list)

    @staticmethod
    def ddict():
        from collections import defaultdict
        return defaultdict(RDF)

    @staticmethod
    def mdict():
        from collections import defaultdict
        return defaultdict(ProbabilisticBooleanNetwork.ddict)

    @staticmethod
    def hamming_distance(x,y):
        x = int(x)
        y = int(y)
        return ZZ(x.__xor__(y)).digits(2).count(1)
    
    @staticmethod
    def load_PBN_file(filename):
        lines = []
        with open(filename, "rU") as f:
            lines = f.readlines()
        M = int(lines[1].strip())
        nf = map(int, lines[3].split())
        N = sum(nf)
        nv = map(int, lines[5].split())
        fvalues = [map(int, l.split()) for l in lines[7:7+N]]
        fdoms = [map(int, l.split()) for l in lines[1+7+N:1+7+2*N]]
        pss = [map(float, l.split()) for l in lines[2+7+2*N:2+7+2*N+M]]
        fs = [BooleanFunction(d,v) for d,v in zip(fdoms,fvalues)]
        i = 0
        bfs = []
        for k in range(M):
            bfs.append(BooleanFunctionSet(fs[i:i+nf[k]], pss[k]))
            i = i+ nf[k]
        return ProbabilisticBooleanNetwork(bfs, fname=filename)

    @staticmethod
    def random_network(nodes, state=None, sync=None, perturbation=None, param=None, min_in_degree=1, max_in_degree=None, min_functions=1, max_functions=5, min_vars=1, max_vars=5):
        V=range(nodes)
        N=2^nodes
        if state is None: state = Integers(N).random_element()
        if sync is None: sync = GF(2).random_element()
        if perturbation is None: 
            perturbation = ProbabilisticBooleanNetwork.__random_uniform_element(["one", "full", "off"])
            if perturbation == "one":
                p,q = [RR.random_element(min=0, max=0.5), RR.random_element(min=0, max=0.1)]
                param = [p,q]
            elif perturbation == "full":
                p = RR.random_element(min=0, max=0.1)
                param = p

        if max_in_degree is None: max_in_degree = nodes
        fs = []
        for _ in V:
            dom = list(Subsets(V, max_in_degree).random_element())
            min_v = min(min_vars, len(dom))
            max_v = min(max(min_v, max_vars), len(dom))

            fs.append(BooleanFunctionSet.random_element(dom, min_v, max_v, min_functions, max_functions))
        
        return ProbabilisticBooleanNetwork(fs, state=state, sync=sync, perturbation = perturbation, param = param)
        

