# K-path-Selection
This repository generates an RRG topology file and also selects K-best optimized paths for each source and destination switch pairs and 
generate the routing file using those selected K best paths.

Compilation :
gcc driver.c model_engine_new.c jellyfish_tr.c jellyfish.h helper.c helper.h topology.h topology_var.h -std=c99 -o kpath


Parameters: nTOR , radix , concentration, bwPS, bwSS, K
Example value: 153 24 1 1 1 87
 
The example RRG topology consists of 153 switches, each switch has 24 switch to switch ports and K=87. Other values needs to be 1 for K path selection.
