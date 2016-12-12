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
    pbn.export_all_graphs()
    pbn.export_all_bscc()
