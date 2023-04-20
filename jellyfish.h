#ifndef JELLYFISH_H
#define JELLYFISH_H

#include "topology.h"


#define JFISH_MAX_NUM_SWITCH 4900
#define JFISH_MAX_NUM_SWITCH_TIMES_K 20000
#define JFISH_MAX_NUM_PATH 400


// Routing schemes for the Jellyfish
// Use block 81 through 100

#define JELLYFISH_SHORTESTPATH_ROUTING    81
#define JELLYFISH_RANDOM1_ROUTING         82
#define JELLYFISH_KPATH_ROUTING           83
#define JELLYFISH_NEWKPATH_ROUTING        84
#define JELLYFISH_ALLPATH_ROUTING         85
#define JELLYFISH_KPATH_LLSKR_MIXED_ROUTING 86
void (* jellyfish_routing_algorithm)(int src, int dst, int *path);
void (* jellyfish_model_routing_algorithm)(int src, int dst, int *len, int *rsrc, int *rdst);
void jellyfish_min_routing(int src, int dst, int *path);
void jellyfish_random1_routing(int src, int dst, int *path);
void jellyfish_model_min_routing(int src, int dst, int *len, int *rsrc, int *rdst);
void jellyfish_model_random1_routing(int src, int dst, int *len, int *rsrc, int *rdst);
int jellyfish_shortestpath_routing(int src, int dst, int *path, int wgraph[MAX_NODE][2*MAX_DEGREE]);
void jellyfish_topology_init(int, int , int, long long int , long long int, int, unsigned int, int, int );
void jellyfish_print_topology();
void jellyfish_dump_topology_binary();
void jellyfish_dump_routing_binary();
void jellyfish_readkpath(int , int , Path []);
void jellyfish_readnewkpath(int , int , Path []);
void jellyfish_print_path(int, int, Path);
void jellyfish_init_path(Path);
void jellyfish_read_routing_to_memory1();
void jellyfish_read_topology();
void jellyfish_calc_avgSP();
int get_r(int swid); //returns network radix of the given switch id   
int get_firstPE(int swid);


/* this routing finds k shortest paths between two nodes */   
int jellyfish_kpath_routing(int src, int dst,             
                              int kpath[JFISH_MAX_NUM_PATH][MAX_PATH_LEN], int k);


/* jellyfish extern variable */

extern int jellyfish_nTOR;   // The Number of TOR Switches
extern int jellyfish_r;      // The number of ports per switch-to-switch link
extern int jellyfish_p;      // The number of servers per TOR switches
extern int jellyfish_x;      // number of left over servers
extern int jellyfish_nPE;    // total servers: nPE = nTOR*p + x;
extern long long int jellyfish_BWPS;
extern long long int jellyfish_BWSS;


#endif
