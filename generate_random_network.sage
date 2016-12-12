import sys
import os
import subprocess

load("import.sage")
if len(sys.argv) < 2:
    print "Expected: Nodes"
    sys.exit()
for fname in sys.argv[1:]:
    print fname
    n = ZZ.random_element(4,8)
    pbn = ProbabilisticBooleanNetwork.random_network(n, max_functions=1)
    pbn.export()
    with open("{}.txt".format(pbn.fname), "w") as f:
        f.write(repr(pbn))
    for m in ["sync", "async"]:
        for pmode, param in zip(["one","full"],[[0.3,0.1],0.03]):
            pbn.set_sync_mode(m == "sync")
            pbn.set_perturbation(mode = pmode, param = param)
            pv = pbn.steady_state()
            save(pbn.bscc_data_plot(pv),"{}.{}.{}.bscc.ps".format(pbn.fname,m,pmode))
            #save(pbn.basin_data_plot(pv),"{}.{}.basin.ps".format(pbn.fname,m))
            save(pbn.sampling_bscc_plot(rel=0.1),"{}.{}.{}.sample.bscc.sparse.ps".format(pbn.fname,m,pmode))
            save(pbn.sampling_bscc_plot(rel=10),"{}.{}.{}.sample.bscc.dense.ps".format(pbn.fname,m,pmode))
