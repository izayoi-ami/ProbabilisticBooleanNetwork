import sys
import os
import subprocess

load("import.sage")
if len(sys.argv) < 2:
    print "Expected: list of PBNs"
    sys.exit()
for fname in sys.argv[1:]:
    print fname
    pbn = ProbabilisticBooleanNetwork.load_PBN_file(fname)
    G = pbn.sync_transition_graph()
    AG = pbn.async_transition_graph()
    dis = pbn.disjunction()
    pbn.export_all_graphs()
    pbn.export_all_bscc()
    dis.export_all_graphs()
    dis.export_all_bscc()
