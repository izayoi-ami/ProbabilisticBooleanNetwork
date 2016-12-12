# ProbabilisticBooleanNetwork
This is a sage library for Probabilistic Boolean Network analysis.
It currently accepts ASSA-PBN file.

**Usage**

Clone this repo by the command `git clone`.

**Requirements**

* Sage 7.3( or above) ([Install](http://www.sagemath.org/))


**Use as a sage library**

1. `cd ProbabilisticBooleanNetwork`
1. `sage`
1. `load("import.sage")`

**Generate N random PBN**

1. `cd ProbabilisticBooleanNetwork`
1. `sage generate_random_network.sage {1..N}`
1. The generated networks can be found in the folder `./data/tmp/`

**Analyze a given ASSA-PBN network**

The pbn file is expected in low-level format.

1. `cd ProbabilisticBooleanNetwork`
1. `sage analyze.sage PBN_FILE`

**TODO**

1. Add colors for BSCC nodes
2. Better API call for network generation
