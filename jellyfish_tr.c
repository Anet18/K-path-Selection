#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "topology.h"


int jellyfish_nTOR;   // The Number of TOR Switches
int jellyfish_r;      // The MAXIMUM number of ports per switch-to-switch link (some switches can have -1 link)
int jellyfish_nPE;    // The number of total processing nodes
//int jellyfish_p;      // The number of servers per TOR switches
int jellyfish_x;
int jellyfish_nPorts;  //total number of ports per switch

int jellyfish_ths = 2;// threshold of path spread value
		      //use KPATH for spreads lower than threshold;

struct Histogram_data{
  int value;
  int frequency;
}; 

long long int jellyfish_BWPS;
long long int jellyfish_BWSS;

#ifndef NEWR_MAX_HOP
#define NEWR_MAX_HOP 3
#endif

#define MAX_JNUM_PATH 1000000
#define MAX_QUEUE_SIZE 1000000

#ifdef _ALL_TO_ALL_PATTERN
#define MAX_PQ_SIZE 50000000  // It could be N*N*K
#else
#define MAX_PQ_SIZE 200
#endif


int lgraph[MAX_NODE][2*MAX_DEGREE]; // temporary local graph
int LARGENUM = 100000;

Path pathcache[JFISH_MAX_NUM_SWITCH][JFISH_MAX_NUM_SWITCH_TIMES_K/10000];

#define MAX_SHORT_PATH_LEN 7
typedef short ShortPath[MAX_SHORT_PATH_LEN-1];

#ifdef _JELLYFISH_LARGE_DATA_STRUCTURE
ShortPath shortpathcache[JFISH_MAX_NUM_SWITCH][JFISH_MAX_NUM_SWITCH*JFISH_MAX_NUM_PATH];
#else
ShortPath shortpathcache[500][15000];
#endif

int jellyfish_K = _REAL_PATH_NUM;

struct TopologyEntry{
  int s;   // Source Node
  int p;   // Port
  int r;   // Remote Node
  int b;   // bandwidth
};

typedef struct TreeNodeType{
  int nid;
  struct TreeNodeType *child[MAX_DEGREE];
  struct TreeNodeType *parent;
  int dist[MAX_DEGREE];
  int pindex;
}Node;

typedef struct EdgeType{
  int s;
  int d;
}Edge;

typedef struct PriorityQueueElementType{
  int snode;
  Node *tnode;
  Edge e;
  int dist;
}PQElement;

typedef struct PriorityQueueType{
  PQElement element[MAX_QUEUE_SIZE + 1];
  int size;
}PriorityQueue;

/*************************************************/
int jellyfish_check_graph();
void jellyfish_print_path(int s, int d, Path path);
void jellyfish_init_path(Path ipath);
void jellyfish_init_graph(int xgraph[MAX_NODE][2*MAX_DEGREE]);
void jellyfish_test_path();
int validateJellyFish();
void jellyfish_build_topology();
int checkStatusForAdjustment(int srcTOR, int tors[], int size);
void printTORSwitch(int t);
void doAdjustment(int srcTOR, int port1, int port2);
int areNeighbors(int srctor, int dsttor);
void printTORS(int tors[], int size);
int getSRCIndex(int tors[], int src, int torsize);
void adjustTORS(int tors[], int i, int *torsize);
int dstFull(int dsttor);
void jellyfish_print_topology_summary();
int jellyfish_compute_K_Shortest_Path(int s, int d, Path kpath[], int K);
void jellyfish_copy_graph(int tgraph[MAX_NODE][2*MAX_DEGREE], int sgraph[MAX_NODE][2*MAX_DEGREE]);
void jellyfish_copy_path(Path dpath, Path spath);
int jellyfish_compute_hop_count(Path tpath);
int jellyfish_find_shortest_path(Path allpath[MAX_JNUM_PATH]);
void jellyfish_disconnect_graph(int lgraph[MAX_NODE][2*MAX_DEGREE], Path kpath[], int dstswitch, int size);
void jellyfish_stitch_path(Path target, Path source, int snode);
void jellyfish_make_path(int s, int d, Path spath);
void jellyfish_read_topology();
void jellyfish_read_routing_to_memory();
void jellyfish_read_routing_to_memory1();
static int getswid(int nid);
static int getlocalid(int nid);
int getfirstPE(int swid);
/************************************************/

int get_firstPE(int swid)
{
  if(swid<totPE|| swid >=totNode) {\
    printf("get_firstPE: TOR id out of range\n");
    exit(0);
  }
  return graph[swid][get_r(swid)+0];
}

static int getswid(int nid)
{
  int p, x;

  if ((nid >= totPE) || (nid <0)) {
    printf("getswid: nid %d out of range.\n", nid);
    exit(0);
  }

  p = jellyfish_nPE / jellyfish_nTOR;
  x = jellyfish_nPE - p * jellyfish_nTOR;

  if (nid < x*(p+1)) return nid / (p+1);
  return ((nid - x*(p+1))/p + x);
}

static int getlocalid(int nid)
{
  int p, x;

  if ((nid >= totPE) || (nid <0)) {
    printf("getswid: nid %d out of range.\n", nid);
    exit(0);
  }

  p = jellyfish_nPE / jellyfish_nTOR;
  x = jellyfish_nPE - p * jellyfish_nTOR;

  if (nid < x*(p+1)) return nid % (p+1);
  return ((nid - x*(p+1)) % p);
}


int get_r(int swid)
{
  int p,x;
  if ((swid >= totNode) || (swid <0)) {
    printf("getswswdegree: invalid node id :%d\n", swid);
    exit(0);
  }
  if(swid<totPE){
    printf("getswswdegree: invalid switch id. :%d is a processing node\n", swid);
    exit(0);
  }
  p = jellyfish_nPE / jellyfish_nTOR;
  x = jellyfish_nPE - p * jellyfish_nTOR;
  if(x==0) return jellyfish_nPorts-p;
  if (swid-totPE < x ) return jellyfish_nPorts - (p+1); //the first 'x' switches have 1 less sw-sw link
  else return jellyfish_nPorts-p;

}

/******************Priority Queue Implementation******************/

void printSubTree(Node *node){
  int i;
  if( node == NULL )return;
  printf("Node = %d ",node->nid);
  if(node->parent != NULL )printf("( Parent = %d )\n",node->parent->nid);
  else printf("( Root Node )\n");
  for(i=0; i<MAX_DEGREE; ++i)
    if(node->child[i] != NULL){
      printf("Child[%d] = %d, ",i,node->child[i]->nid);
      printf("\n");
    }
  for(i=0; i<MAX_DEGREE; ++i)
    if(node->child[i] != NULL)
      printSubTree(node->child[i]);
}

void initTreeNode(Node *node, int id){
  int i;
  node->nid = id;
  for(i=0; i<MAX_DEGREE; ++i){
    node->child[i] = NULL;
    node->dist[i] = -1;
  }
  node->parent = NULL;
  node->pindex = 0;
}

int isQEmpty(PriorityQueue *pq){
  return (pq->size == 0);
}

void printPQElement(PQElement pqe){
  printf("snode = %d e = (%d->%d) dist = %d\n",
  pqe.snode,pqe.e.s,pqe.e.d,pqe.dist);
}

void printPQ(PriorityQueue pq){
  int i;
  for(i=1; i< pq.size+1; ++i){
    printf("element[%d] = ",i);
    printPQElement ( pq.element[i] );
  }
}

void initPQElement( PQElement *pqe ){
  pqe->snode = -1;
  pqe->tnode = NULL;
  pqe->e.s   = -1;
  pqe->e.d   = -1;
  pqe->dist  = -1;  
}

void pushPQ( PriorityQueue *pq, PQElement x ){
  int hole = ++pq->size;
  for(; hole > 1 && x.dist < pq->element[ hole/2 ].dist; hole /= 2)
    pq->element[ hole ] = pq->element[ hole/2 ];
  pq->element[ hole ] = x;
}

void updatePQ( PriorityQueue *pq, int snode, int dist ){
  int i,child;
  PQElement pqe;
  
  initPQElement( &pqe );
  
  for(i=0; i<= MAX_QUEUE_SIZE; ++i)if( pq->element[i].snode == snode )break;
  pq->element[ i ].dist = dist;
  
  if( dist < pq->element[i/2].dist ){
    // Need to percolate up
    for(; i > 1 && dist < pq->element[ i/2 ].dist; i /= 2){
      pqe = pq->element[ i ];
      pq->element[ i ] = pq->element[ i/2 ];
      pq->element[ i/2 ] = pqe;
    }
  }
  else{
    for(; i*2 <= pq->size; i = child){
      child = i*2;
      if( child != pq->size && 
                         pq->element[ child+1].dist < pq->element[child].dist )
      child++;
      if(pq->element[ child ].dist < dist){
        pqe = pq->element[ i ];
        pq->element[ i ] = pq->element[ child];
        pq->element[ child] = pqe;
      }
      else
        break;
    }
  }
}

void percolateDown( PriorityQueue *pq, int hole ){
  int child;
  PQElement tmp = pq->element[ hole ] ;
  for(; hole*2 <= pq->size; hole = child){
    child = hole*2;
    if( child != pq->size && pq->element[ child+1].dist < pq->element[child].dist )
      child++;
    if(pq->element[ child ].dist < tmp.dist)
      pq->element[ hole ] = pq->element[ child];
    else
      break;
  }
  pq->element[hole] = tmp;
}

PQElement popPQ( PriorityQueue *pq ){
  PQElement min = pq->element[1];
  
  pq->element[1] = pq->element[ pq->size-- ];
  percolateDown( pq, 1 );
  return min;
}

void initPQ(PriorityQueue *pq){
  int i;
  for(i=0; i<MAX_QUEUE_SIZE+1; ++i) initPQElement ( &pq->element[i] );
  pq->size = 0;
}
/************************************************/

void jellyfish_print_topology_summary(){
  printf("Number of terminal/processing nodes = %d.\n",totPE);
  printf("Number of switch nodes = %d.\n",totSE);
  printf("Total number of nodes = %d.\n",totNode);
  printf("jellyfish_r: %d, jellyfish_p: %d, jellyfish_x: %d\n", 
	 jellyfish_r, totPE/totSE, jellyfish_x);
  printf("Total number of ports = %d.\n",jellyfish_nPorts);



}

void jellyfish_copy_graph(int tgraph[MAX_NODE][2*MAX_DEGREE], int sgraph[MAX_NODE][2*MAX_DEGREE]){
  int i,j;
  for(i=0; i < MAX_NODE ; ++i)
   for(j=0; j < MAX_DEGREE; ++j)
     tgraph[i][j] = sgraph[i][j];
}

void jellyfish_print_topology(){
  int i, j, k;
  
  jellyfish_print_topology_summary();
  
  for(i=0; i<totNode; i++){
    if(i < totPE)
      printf("Terminal Node %d, ",i);
    else{
      printf("Switch Node %d, ",i);
    }
    printf("Connects to :\n");
    for(j=0; graph[i][j] != -1; j++){
      k = graph[i][j];
      if(k < totPE)
        printf("\tterminal Node %d (bw = %lld)\n",k,bandwidth[i][j]);
      else
        printf("\tswitch Node %d (bw = %lld)\n",k,bandwidth[i][j]);
      printf("\tincoming: %d\n",graph[i][j+MAX_DEGREE]);
    }
  }
  
  fflush(0);
}

void adjustTORS(int tors[], int i, int *torsize){
  for(;i < *torsize; i++)tors[i] = tors[i+1];
  (*torsize)--;
}

int getSRCIndex(int tors[], int src, int torsize){
  int i;
  for(i=0; i<torsize; i++)if(src == tors[i])return i;
  return -1;
}

void printTORS(int tors[], int size){
  int i;
  for(i=0; i< size; i++)
    printf("tors[%d] = %d\n",i,tors[i]);
}

int areNeighbors(int srctor, int dsttor){
  int i;
  
  if(srctor == dsttor) return TRUE;
  
  for(i=0; i<get_r(srctor) ;++i)
    if( graph[srctor][i] == dsttor ) return TRUE;
  
  /*
  for(i=0; i<jellyfish_r;++i)
    if( graph[dsttor][i] == srctor ) return TRUE;
  */
  
  return FALSE;
}

int dstFull(int dsttor){
  int i, k;
  int count;

  count = 0;
  for(i=0; i<totSE; i++)
    for(k=0; k< get_r(totPE+i); k++)
      if(graph[totPE + i][k] == dsttor) ++count;


  if(count >= get_r(dsttor)){
    return TRUE;
  }
  
  //printf("c=%d\n",count);
  return FALSE;
}

void doAdjustment(int srcTOR, int port1, int port2){
  int i;
  int notdone = TRUE;
  int a1TOR, a2TOR;
  
  //  printf("Adjusting for srcTOR[%d].\n",srcTOR);
  
  while( notdone ){
    do{   //randomly select a switch currently connected to srcTOR
      a1TOR = totPE + rand() % jellyfish_nTOR;
    }while(areNeighbors(srcTOR, a1TOR));
    
    for(i=0; i< get_r(a1TOR); ++i){
      a2TOR = graph[a1TOR][i];
      if( ! areNeighbors(srcTOR, a2TOR) ){
        //notdone = FALSE;
        graph[srcTOR][port1] = a1TOR;
        bandwidth[srcTOR][port1] = jellyfish_BWSS;
        graph[a1TOR][i] = srcTOR;
        
        //printf("graph[%d][%d] = %d\n",srcTOR,port1,a1TOR);
        //printf("graph[%d][%d] = %d\n",a1TOR,i,srcTOR);
        break;
      }      
    }
    // Find the port on a2TOR at which a1TOR is connected to it
    for(i=0; i<get_r(a2TOR); ++i){
      if(graph[a2TOR][i] == a1TOR){
        notdone = FALSE;
        graph[srcTOR][port2] = a2TOR;
        bandwidth[srcTOR][port2] = jellyfish_BWSS;
        graph[a2TOR][i] = srcTOR;
        
        //printf("graph[%d][%d] = %d\n",srcTOR,port2,a2TOR);
        //printf("graph[%d][%d] = %d\n",a2TOR,i,srcTOR);
        
        break;
      }
    }
  }
}

void printTORSwitch(int t){
  int i;
  for(i=0; i<get_r(t);i++)
    printf("graph[%d][%d] = %d\n",t,i,graph[t][i]);
}

int checkStatusForAdjustment(int srcTOR, int tors[], int size){
  int i=0;
  int allNeighbor = TRUE;
  while( (i< size) && allNeighbor ){
    allNeighbor = areNeighbors( srcTOR, (totPE + tors[i++]) );
  }
  
  return allNeighbor;
}

void jellyfish_build_topology(){
  int i,j,k,nattempt;
  int srcTOR,dstTOR, src_r, dst_r;
  int torsize = jellyfish_nTOR;
  int tors[torsize+1];
  int srcindex, count, dstindex;
  int doadjust;
  int port1, port2;
  int totEdges = 0;
  
  totSE   = jellyfish_nTOR;
  totPE = jellyfish_nPE;
  totNode = totSE + totPE;
  nprocs  = totPE;
  
  if(totNode > MAX_NODE){
    printf("ERROR:The number of Nodes ( = %d ) exceeds the MAX_NODE = %d.\n",totNode,MAX_NODE);
    exit(0);
  }
  
  // Initialize the graph
  for (i=0; i<MAX_NODE; i++)for (j=0; j<MAX_DEGREE; j++)graph[i][j] = -1;
  
  for (i=0; i<MAX_NODE; i++)
    for (j=0; j<MAX_DEGREE; j++){
      graph_m1[i][j] = 0;
    }

  for(i=0; i<torsize; ++i)tors[i] = i;
  tors[torsize] = -1;
  
  for(i=0; i<totPE; i++){
    srcTOR = totPE + getswid(i); // + i / jellyfish_p;
    src_r = get_r(srcTOR);
    graph[i][0] = srcTOR;
    bandwidth[i][0] = jellyfish_BWPS;
    
    //    graph[srcTOR][jellyfish_r+i%jellyfish_p] = i;
    //    bandwidth[srcTOR][jellyfish_r+i%jellyfish_p] = jellyfish_BWPS;
    graph[srcTOR][src_r+getlocalid(i)] = i;
    bandwidth[srcTOR][src_r+getlocalid(i)] = jellyfish_BWPS;
  }
  
  for(i=0; i<totSE; i++){
    srcTOR = i+totPE;
    src_r = get_r(srcTOR);
    srcindex = getSRCIndex(tors,i,torsize);
    
    //    printf("i = %d srctor = %d srcindex = %d\n",i,srcTOR,srcindex);
    //    printTORS(tors, torsize);
    
    if( srcindex != -1 ){
      for(j=0; j<src_r; j++){
        doadjust = FALSE;
        if( graph[srcTOR][j] == -1 ){
	  //          printf("Working on Switch[%d][%d].\n",srcTOR,j);
          nattempt = 0;
          do{
            dstindex = rand() % torsize;
            dstTOR = totPE + tors[ dstindex ];
            dst_r = get_r(dstTOR);
            ++nattempt;
            
            if(nattempt >= 100){
              if(checkStatusForAdjustment(srcTOR, tors, torsize)){
                doadjust = TRUE;
                break;
              }
            }
          }while( areNeighbors(srcTOR, dstTOR) && (! doadjust));
          
          if(doadjust){
	    //            printf("Adjustment required for srcTOR[%d].\n",srcTOR);
            port1 = -1;
            port2 = -1;
            for(k=0; k<src_r; k++){
              if(graph[srcTOR][k] == -1){
                if(port1 == -1)port1 = k;
                else if(port2 == -1){
                  port2 = k;
                  break;
                }
              }
            }
	    //            printf("Port1 = %d, port2 = %d\n",port1,port2);
            doAdjustment(srcTOR, port1, port2);
            
            continue;
          }
          
	  //          printf("\tdstTOR[%d] ",dstTOR);
        
          for(k=0; k<dst_r; k++){
            if(graph[dstTOR][k] == -1){
	      //              printf("\tPort[%d] on dstTOR[%d][%d] chosen.\n",k,dstTOR,k);
              graph[srcTOR][j] = dstTOR;
              bandwidth[srcTOR][j] = jellyfish_BWSS;
              graph[dstTOR][k] = srcTOR;
              bandwidth[dstTOR][k] = jellyfish_BWSS;
	      totEdges+=2;
              break;
            }
          }
          count = 0;
          for(k=0; k<dst_r; k++)if(graph[dstTOR][k] != -1) ++count;
          if(count >= dst_r){
            adjustTORS(tors, dstindex, &torsize);
	    //            printf("\tdstTOR[%d] is done. tors after removal.\n",dstTOR);
	    //            printTORS(tors, torsize);
          }
        }
      }
      // the srcindex might have been changed. hence find it again
      srcindex = getSRCIndex(tors,i,torsize);
      adjustTORS(tors, srcindex, &torsize);
      //      printf("\tsrcTOR[%d] is done. tors after removal.\n",srcTOR);
      //      printTORS(tors, torsize);
    }
  }
  printf("total edge count: %d\n",totEdges);
}

void jellyfish_build_topology_unidirectional(){
  int i, j, k, pp, nattempt;
  int srcTOR,dstTOR, intTOR, src_r, dst_r, int_r;
  int torsize = jellyfish_nTOR;
  int tors[torsize+1];
  int srcindex, done, dstindex, intindex;
  int totEdges = 0;

  totSE   = jellyfish_nTOR;
  //  totPE   = jellyfish_nTOR * jellyfish_p;
  totPE = jellyfish_nPE;
  totNode = totSE + totPE;
  nprocs  = totPE;
  
  if(totNode > MAX_NODE){
    printf("ERROR:The number of Nodes ( = %d ) exceeds the MAX_NODE = %d.\n",totNode,MAX_NODE);
    exit(0);
  }
  
  // Initialize the graph
  for (i=0; i<MAX_NODE; i++)for (j=0; j<2*MAX_DEGREE; j++)graph[i][j] = -1;
  
  for(i=0; i<torsize; ++i)tors[i] = i;
  tors[torsize] = -1;
  
  for(i=0; i<totPE; i++){
    srcTOR = totPE + getswid(i); // + i / jellyfish_p;
    src_r = get_r(srcTOR);
    graph[i][0] = srcTOR;
    bandwidth[i][0] = jellyfish_BWPS;
    graph[i][MAX_DEGREE] = srcTOR;

    //    graph[srcTOR][jellyfish_r+i%jellyfish_p] = i;
    //    bandwidth[srcTOR][jellyfish_r+i%jellyfish_p] = jellyfish_BWPS;
    graph[srcTOR][src_r+getlocalid(i)] = i;
    bandwidth[srcTOR][src_r+getlocalid(i)] = jellyfish_BWPS;
    graph[srcTOR][src_r+getlocalid(i)+MAX_DEGREE] = i;
  }
  
  for(i=0; i<totSE; i++){
    srcTOR = i+totPE;
    src_r = get_r(srcTOR);
    srcindex = getSRCIndex(tors,i,torsize);
    
    //    printf("i = %d srctor = %d srcindex = %d\n",i,srcTOR,srcindex);
    //    printTORS(tors, torsize);
    
    if( srcindex != -1 ){
      for(j=0; j<src_r; j++){
        //doadjust = FALSE;
        if( graph[srcTOR][j] == -1 ){
	  //          printf("Working on Switch[%d][%d].\n",srcTOR,j);
          nattempt = 0;
          do{
            dstindex = rand() % torsize;
            dstTOR = totPE + tors[ dstindex ];
            ++nattempt;
            
            /*if(nattempt >= 100){
              if(checkStatusForAdjustment(srcTOR, tors, torsize)){
                doadjust = TRUE;
                break;
              }
	      }*/
          }while( (areNeighbors(srcTOR, dstTOR) && nattempt < 2 * totSE - 1) || dstFull(dstTOR));
          

	  if(areNeighbors(srcTOR, dstTOR) && nattempt > 2 * totSE - 1) {
	    done = 0;
	    do{
	      intindex = rand() % torsize;
	      intTOR = totPE + tors[ intindex ];
              int_r = get_r(intTOR);
	      for(k=0; k<int_r; k++){
		if(graph[intTOR][k] != -1 && !areNeighbors(srcTOR, graph[intTOR][k]) && !areNeighbors(intTOR, dstTOR)) {
		  done = 1;
		  graph[srcTOR][j] = graph[intTOR][k];
		  bandwidth[srcTOR][j] = jellyfish_BWSS;
		  for(pp=MAX_DEGREE;pp<2*MAX_DEGREE;pp++)
		    if(graph[graph[intTOR][k]][pp] == -1){
		      graph[graph[intTOR][k]][pp] = srcTOR;
		      break;
		    }
		  graph[intTOR][k] = dstTOR;
		  bandwidth[intTOR][k] = jellyfish_BWSS;
		  for(pp=MAX_DEGREE;pp<2*MAX_DEGREE;pp++)
		    if(graph[dstTOR][pp] == -1){
		      graph[dstTOR][pp] = intTOR;
		      break;
		    }
		  
		  totEdges++;
		  break;
		}
	      }
	    }while(!done);
	    
	  }
	  else{
	    graph[srcTOR][j] = dstTOR;
	    bandwidth[srcTOR][j] = jellyfish_BWSS;
	    for(pp=MAX_DEGREE;pp<2*MAX_DEGREE;pp++)
	      if(graph[dstTOR][pp] == -1){
		graph[dstTOR][pp] = srcTOR;
		break;
	      }
	    totEdges++;
	  }
        
          
        }
      }
    }
  }
  printf("total edge count: %d\n",totEdges);
  //jellyfish_print_topology();
}

void jellyfish_read_topology(){
  char topFile[100];
  FILE *fp = NULL;
  int i,j;
  unsigned int stat;
  
  struct TopologyEntry te;
  
  sprintf(topFile, "jellyfish_topology_%d_%d_%d.bin",jellyfish_nTOR, jellyfish_nPE, jellyfish_r);
  
  printf("Reading the topology from the file %s\n",topFile);
  
  if( (fp = fopen(topFile, "rb")) == NULL ){
    printf("Could not open file %s for reading.\n",topFile);
    exit(0);
  }
  
  totSE   = jellyfish_nTOR;
  totPE   = jellyfish_nPE;
  totNode = totSE + totPE;
  nprocs  = totPE;
  
  if(totNode > MAX_NODE){
    printf("ERROR:The number of Nodes ( = %d ) exceeds the MAX_NODE = %d.\n",totNode,MAX_NODE);
    exit(0);
  }
  
  // Initialize the graph
  for (i=0; i<totNode; i++)for (j=0; j<MAX_DEGREE; j++)graph[i][j] = -1;
  
  while(( stat = fread(&te, sizeof(te),1,fp))){
    //printf("s = %d p = %d r = %d b = %d\n",te.s,te.p,te.r,te.b);
    if( (te.s < 0) || (te.p < 0) || (te.r < 0) || (te.b < 0)){
      printf("Error while reading the topology file.\n");
      exit(0);
    }
    graph[te.s][te.p] = te.r;
    bandwidth[te.s][te.p] = te.b;
  }
  
  if(fp != NULL)fclose(fp);
  
  printf("Reading the topology from the file %s done\n",topFile);
}

int validateJellyFish(){

  int p;

  p = jellyfish_nPE / jellyfish_nTOR;
  if (jellyfish_nPE % jellyfish_nTOR != 0) p++;

  /*if( (jellyfish_nTOR * jellyfish_r)%2 == 1 ){
    printf("The product of nTOR and r should be even %d %d.\n", jellyfish_nTOR, jellyfish_r);
    return 1;
  }
  else */
  if((jellyfish_nTOR-1) < jellyfish_r){
    printf("(nTOR-1) should be >= r.\n");
    return 1;
  }
  else if((jellyfish_r + p) > MAX_DEGREE){
    printf("(jellyfish_r + jellyfish_p) exceeds MAX_DEGREE.\n");
    return 1;
  }
  return 0;
}

void jellyfish_test_path(){
  int i,j;
  Path path;
  
  for(i=0; i<totPE; ++i){
    for(j=0; j<totPE; ++j){
      if( i != j ){
        if( jellyfish_shortestpath_routing( i,j,path, graph ) ){
          jellyfish_print_path(i, j, path);
        }
        else{
          printf("Shortest path does not exists between S[%d]=>D[%d].\n",i,j);
        }
      }
    }
  }
}

void jellyfish_param_init( int n, int totports, int totservers, long long int bwPS, long long int bwSS, int rc, unsigned int seed, int dump){
  //  int i,j;
  
  //for(i=0; i<totNode; i++)for (j=0; graph[i][j] != -1; j++)load_graph[i][j] = 0;
  
  srand(seed);
  
  jellyfish_nTOR = n;
  //jellyfish_p = totservers/n;  
  //keep these variables for backward compatibility
  jellyfish_x = totservers - (n*( totservers/n)); 
  //keep these variables for backward compatibility
  jellyfish_nPE = totservers;
  jellyfish_r    = totports - (totservers/n) - ((jellyfish_x)?1:0);
  jellyfish_nPorts = totports;
  jellyfish_BWPS = bwPS;
  jellyfish_BWSS = bwSS;
  
  totSE   = jellyfish_nTOR;
  totPE   = jellyfish_nPE;
  totNode = totSE + totPE;
  nprocs  = totPE;
  
  topology = JELLYFISH;
  

  if(validateJellyFish())exit(0);
  
  if( dump ){
    //printf("building jellyfish topology\n\n");
    //jellyfish_build_topology();
  }
  else{
    jellyfish_read_topology();
    if (rc ==JELLYFISH_KPATH_ROUTING) 
      jellyfish_read_routing_to_memory1();
    else if ( rc ==JELLYFISH_NEWKPATH_ROUTING) 
      jellyfish_read_routing_to_memory1();
  }
  /*
  if ( jellyfish_check_graph() ){
    printf("jellyfish:topology connectivity test **FAILED**.\n");
    exit(0);
  }
  else{
    printf("jellyfish:topology connectivity test passed.\n");
    //printf("jellyfish:topology connectivity test **FAILED**.\n");
  }
  */
  routing  = rc;
  if( routing == JELLYFISH_SHORTESTPATH_ROUTING ){
    routing_algorithm = jellyfish_routing_algorithm = jellyfish_min_routing;
    model_routing_algorithm = jellyfish_model_routing_algorithm = jellyfish_model_min_routing;
  }
  else if( routing == JELLYFISH_RANDOM1_ROUTING ){
    routing_algorithm = jellyfish_routing_algorithm = jellyfish_random1_routing;
    model_routing_algorithm = jellyfish_model_routing_algorithm = jellyfish_model_random1_routing;
  }
  else if( routing == JELLYFISH_KPATH_ROUTING ){
    
  } else if( routing == JELLYFISH_NEWKPATH_ROUTING ){
    
  }
  else{
    printf("Unsupported routing scheme for Jellyfish.\n");
  }
  
  printf("Simulating JELLYFISH with parameters:\n");
  printf("N(number of TOR) = %d\n",jellyfish_nTOR);
  printf("nPorts(Total Number of switch ports ) = %d\n",jellyfish_nPorts);
  printf("r(Number of switch ports connecting other switches) = %d\n",jellyfish_r);
  printf("NPE(Number of servers) = %d\n",jellyfish_nPE);
  printf("BWPS(Bandwidth processor to switch) = %lld\n",jellyfish_BWPS);
  printf("BWPS(Bandwidth switch to switch) = %lld\n",jellyfish_BWSS);
  
  //jellyfish_print_topology();
  
  jellyfish_print_topology_summary();
}

//void jellyfish_topology_init( int n, int r, int p, int x, long long int bwPS, long long int bwSS, int rc, unsigned int seed, int dump, int directional){
void jellyfish_topology_init( int n, int totports, int totservers, long long int bwPS, long long int bwSS, int rc, unsigned int seed, int dump, int directional){
  int i,j, k;
  
  
  srand(seed);
  

  //jellyfish_nTOR = n;
  //jellyfish_p = p;
  //jellyfish_x = x;
  //jellyfish_nPE = n*p + x;
  //jellyfish_r    = r;

  jellyfish_nTOR = n;
  //jellyfish_p = totservers/n;  //keep these variables for backward compatibility
  jellyfish_x = totservers - (n*( totservers/n)); //keep these variables for backward compatibility
  jellyfish_nPE = totservers;
  jellyfish_r    = totports - ( totservers/n) - ((jellyfish_x)?1:0);
  jellyfish_nPorts = totports;

  jellyfish_BWPS = bwPS;
  jellyfish_BWSS = bwSS;

  totSE   = jellyfish_nTOR;
  totPE = jellyfish_nPE;
  totNode = jellyfish_nTOR + totPE;
  nprocs = totPE;

  for (i=0; i<totSE; i++)
    for (j=0; j<totSE*JFISH_MAX_NUM_PATH; j++)
      for (k=0; k<MAX_SHORT_PATH_LEN; k++) shortpathcache[i][j][k] = -1;

  //  for(i=0; i<totNode; i++)for (j=0; graph[i][j] != -1; j++)load_graph[i][j] = 0;

  topology = JELLYFISH;
  
  if(validateJellyFish()) exit(0);
  
  if( dump ){
    if(directional == 1)
      jellyfish_build_topology_unidirectional();
    else
      jellyfish_build_topology();
    printf("Jellyfish topology build successful\n\n");
    //jellyfish_print_topology();
  }
  else{
    jellyfish_read_topology();
    if (rc ==JELLYFISH_KPATH_ROUTING) 
      jellyfish_read_routing_to_memory1();
    else if ( rc ==JELLYFISH_NEWKPATH_ROUTING) 
      jellyfish_read_routing_to_memory1();
    else if (rc == JELLYFISH_SHORTESTPATH_ROUTING)
      jellyfish_read_routing_to_memory1();
  }
  
  if ( jellyfish_check_graph() ){
    printf("jellyfish:topology connectivity test **FAILED**.\n");
    exit(0);
  }
  else{
    printf("jellyfish:topology connectivity test passed.\n");
    //printf("jellyfish:topology connectivity test **FAILED**.\n");
  }
  
  routing  = rc;
  if( routing == JELLYFISH_SHORTESTPATH_ROUTING ){
    routing_algorithm = jellyfish_routing_algorithm = jellyfish_min_routing;
    model_routing_algorithm = jellyfish_model_routing_algorithm = jellyfish_model_min_routing;
    routingType = SINGLEPATH;
  }
  else if( routing == JELLYFISH_RANDOM1_ROUTING ){
    routing_algorithm = jellyfish_routing_algorithm = jellyfish_random1_routing;
    model_routing_algorithm = jellyfish_model_routing_algorithm = jellyfish_model_random1_routing;
    routingType = SINGLEPATH;
  }
  else if( routing == JELLYFISH_KPATH_ROUTING ){
    routingType = MULTIPATH;
  } else if( routing == JELLYFISH_NEWKPATH_ROUTING ){
    routingType = MULTIPATH;
  } else if( routing == JELLYFISH_KPATH_LLSKR_MIXED_ROUTING ){
    routingType = MULTIPATH;
  }
  else{
    printf("Unsupported routing scheme for Jellyfish.\n");
  }
  
  if(directional == 0) printf("Simulating JELLYFISH with parameters:\n");
  else printf("Simulating UNIJELLYFISH with parameters:\n");
  printf("N(number of TOR) = %d\n",jellyfish_nTOR);
  printf("r(Number of switch ports connecting other switches) = %d\n",jellyfish_r);
  printf("nPorts(Total Number of switch ports ) = %d\n",jellyfish_nPorts);
  printf("nPE(Number of servers) = %d\n",jellyfish_nPE);
  printf("BWPS(Bandwidth processor to switch) = %lld\n",jellyfish_BWPS);
  printf("BWPS(Bandwidth switch to switch) = %lld\n",jellyfish_BWSS);
  
  //jellyfish_print_topology();
  
  jellyfish_print_topology_summary();
}

typedef struct FTableValueType{
  int fromNode;
  int cost;
}FTablevalue;

typedef struct QueueType{
  int element[MAX_NODE];
  int head;
  int tail;
  int size;
}Queue;

void reset(Queue *q){
  q->head = 0;
  q->tail = 0;
  q->size = 0;
}

void printQ(Queue *q){
  for(;q->head < q->tail; q->head++) 
    printf("Queue[%d] = %d\n",q->head, q->element[q->head]);
}

int push(Queue *q, int x){
  if( q->size < MAX_QUEUE_SIZE ){
    q->element[q->tail++] = x;
    q->size++;
    //printf("Pushed value = %x, now head = %d, tail = %d,size = %d\n",
    //x,q->head,q->tail,q->size);
  }
  else{
    printf("Queue has reached capacity, capacity = %d currentsize = %d.\n",MAX_QUEUE_SIZE, q->size);
    //printQ( q );
    
    exit(0);
  }
  return 0;
}

int pop(Queue *q){
  if( q->size > 0 ){
    q->element[ q->head++ ] = -1;
    q->size--;
    
    //printf("Popped now head = %d, tail = %d,size = %d\n",q->head,q->tail,q->size);
  }
  else{
    printf("Queue is already empty.\n");
    exit(0);
  }
  
  return 0;
}

int getSize(Queue *q){
  return q->size;
}

int isEmpty(Queue *q){
  return (q->size == 0);
}

int getFront(Queue *q){
  return q->element[q->head];
}

void jellyfish_init_graph(int xgraph[MAX_NODE][2*MAX_DEGREE]){
  int i,j;
  for(i=0; i<MAX_NODE; ++i)
    for(j=0; j<MAX_DEGREE; ++j)
      xgraph[i][j] = -1;
}

void jellyfish_init_path(Path ipath){
  int i;
  for(i=0; i< MAX_PATH_LEN; ++i) ipath[i] = -1;
}

void jellyfish_print_path(int s, int d, Path path){
  int i;
  
  printf("Path between Src[%d] ==> Dst[%d]\t",s,d);
  for(i=0; (path[i] != -1) && (i<MAX_PATH_LEN); ++i) printf("%d,",path[i]);
  printf("  hops: %d \n",i-3);
}

// This routine checks the connectivity of the graph
int jellyfish_check_graph(){
  int i,j;
  int p, x;
  int pp;

  p = jellyfish_nPE / jellyfish_nTOR;
  x = jellyfish_nPE - jellyfish_nTOR * p;  

  for(i=0; i< totNode; ++i){
    if(i < totPE){
      // Check the connectivity of the terminals
      
      // terminal should always be connected to a switch
      if( graph[i][0] == -1 ) {
        printf("terminal Node[%d] not connected to switch.\n",i);
        return -1;         
      }
      
      // terminal should not be connected to terminal
      if( graph[i][0] < totPE ) {
        printf("terminal Node[%d] connected to terminal Node[%d].\n",i,graph[i][0]);
        return -1;
      }
    }
    else{
      // Check the connectivity of the switches
      
      // Switch[0] should always be connected to a switch
      if( graph[i][0] == -1 ) {
        printf("switch Node[%d] not connected to switch.\n",i);
        return -1;
      }
      
      // Switch[ 0 - r] should always be connected to a switch
      for(j=0; j < get_r(i); ++j)
        if( graph[i][j] == -1 ) {
          printf("Switch[%d]Port[%d] not connected to switch.\n",i,j);
          return -1;
        }

      // Switch should have an incoming link from its neighbor
      for(j=0; j < get_r(i); ++j)
        if( areNeighbors(graph[i][j],i) == FALSE ) {
          printf("Switch[%d]Port[%d] has outgoing link to switch %d, but has no incoming link from it.\n",i,j,graph[i][j]);
          return -1;
        }
      
      // Switch[ 0 - r] should always be connected to a switch not terminal
      for(j=0; j < get_r(i); ++j)
        if( graph[i][j] < totPE ) {
          printf("Switch[%d]Port[%d] should not be connected to terminal.\n",i,j);
          return -1;
        }

      if (i < totPE + x) pp = p+1;
      else pp = p;
        
      // Switch[ r - (r+p)] should always be connected to a terminal
      for(j=0; j < pp; ++j)
        if( graph[i][get_r(i) + j] == -1 ) {
          printf("Switch[%d]Port[%d] not connected to terminal.\n",i,j);
          return -1;
        }
        
      // Switch[ r - (r+p)] should always be connected to a terminal not switch
      for(j=0; j < pp; ++j)
        if( graph[i][get_r(i) + j] >= totPE ) {
          printf("Switch[%d]Port[%d] should not be connected to switch.\n",i,j);
          return -1;
        }
    }
  }
  return 0;
}

int jellyfish_shortestpath_routing(int srcswitch, int dstswitch, Path path, int wgraph[MAX_NODE][2*MAX_DEGREE]){
  Queue q;
  int swid,rid;
  int i, 
      j, 
      ts,   // this source
      nhop=0;
  int unvisited[totSE];
  FTablevalue ftable[totSE];
  int cost, fromNode;
  Path tpath;
  int pathexist = TRUE;
  
  //printf("Computing shortest path for (%d => %d).\n",srcswitch,dstswitch);
  // These switch ids are the ids among the switch network (excluding the 
  // processing nodes)
  if((srcswitch < 0) || (srcswitch >= totSE) || (dstswitch < 0) || (dstswitch >= totSE)){
    printf("Either srcswitch=%d or dstswitch=%d or both are not valid switch nodes.\n",srcswitch,dstswitch);
    exit(0);
  }
  
  // Compute the relative id of the switch among the network of switches
  //srcswitch = s / jellyfish_p;
  //dstswitch = d / jellyfish_p;
  
  //printf("SRC-Switch = %d , DST-Switch = %d\n",srcswitch,dstswitch);
  
  jellyfish_init_path( path );
  jellyfish_init_path( tpath );
  
  for(i=0; i< totSE; ++i)unvisited[i] = 1;
  
  for(i=0; i< totSE; ++i){
    if(i != srcswitch) {
      ftable[i].fromNode = -1;
      ftable[i].cost = LARGENUM;
    }
    else {
      ftable[srcswitch].fromNode = srcswitch;
      ftable[srcswitch].cost = 0;
    }
  }
  
  if(srcswitch == dstswitch){
    path[0] = srcswitch + totPE;
    
    return TRUE;
  }  
  
  reset(&q);
  push(&q, srcswitch);
  
  // While the queue is not empty
  while( ! isEmpty(&q) ){
    // get the front element
    ts = getFront(&q);
    
    // pop out the front element
    pop( &q );
    
    swid = ts + totPE;
    
    //printf("Exploring at ts = %d(Switch[%d]):univisted[%d] = %d.\n",
    //ts,swid,ts, unvisited[ts]);
    
    // if the node is node yet visited
    if( unvisited[ts] ){
      // mark it as visited
      unvisited[ts] = 0;
      
      // find the cost in number of hops
      cost = ftable[ ts ].cost;
      fromNode = ftable[ ts ].fromNode; fromNode = fromNode;
      
      //printf("\tCost to Switch[%d] = %d from Node = %d.\n",
      //(ts+totPE),cost,fromNode);
      
      // for each node connected to it
      for(j=0; j < get_r(swid); ++j){
        if(wgraph[swid][j] != -1){
          if( wgraph[swid][j] >= totPE ){
            rid = wgraph[swid][j] - totPE;
            //printf("Considering rid = %d Switch[%d]\n",rid,wgraph[swid][j]);
            if( (cost+1) < ftable[ rid ].cost){
              // This path should be accepted
              ftable[ rid ].cost = cost + 1;
              ftable[ rid ].fromNode = ts;
            }
        
            if( unvisited[ rid ] ){
              push(&q, rid);
            }
          }
        }
      }
    }
  }
  
  pathexist = !unvisited[dstswitch];
  
  if( pathexist ){
    ts = dstswitch;
    //i = 1;
    i = 0;
    //tpath[0] = d;
    while( ts != srcswitch ){
      tpath[i++] = ts+totPE;
      ts = ftable[ts].fromNode;
      ++nhop;
    }
    tpath[i++] = srcswitch + totPE;
    //for(;i>=0;i--)printf("tpath[%d] = %d, ",i,tpath[i]);
    --i;    
    for(j=0; i>=0; i--)path[j++] = tpath[i];
    //for(j=0; j<(i+1);++j) printf("path[%d] = %d, ",j,tpath[j]);
  }
  return pathexist;
}

void jellyfish_model_min_routing(int src, int dst, int *len, int *rsrc, int *rdst){
  int path[MAX_PATH_LEN];
  int i;

  jellyfish_min_routing(src, dst, path);
  
  for (i=0; path[i+1] != -1; i++) {
    rsrc[i] = path[i];
    rdst[i] = path[i+1];
  }
  *len = i;
}

void jellyfish_model_random1_routing(int src, int dst, int *len, int *rsrc, int *rdst){
  int path[MAX_PATH_LEN];
  int i;

  jellyfish_random1_routing(src, dst, path);
  
  for (i=0; path[i+1] != -1; i++) {
    rsrc[i] = path[i];
    rdst[i] = path[i+1];
  }
  *len = i;
}

int jellyfish_compute_hop_count(Path tpath){
  int i;
  //printPath(tpath);
  for(i=0; tpath[i] != -1; ++i);
  //printf("returning %d\n",i);
  return i;
}

int jellyfish_find_shortest_path(Path allpath[MAX_JNUM_PATH]){
  int i,mi,mh = LARGENUM,hop;
  
  for(i=0; i<MAX_JNUM_PATH; ++i)if( allpath[i][0] != -1 ) break;
  
  for(; i < MAX_JNUM_PATH; ++i){
    //printf("Allpath[%d]=",i);
    //jellyfish_print_path(-1,-1,allpath[i]);
    if( allpath[i][0] != -1 ){
      hop = jellyfish_compute_hop_count(allpath[i]);
      if( hop < mh ){
        mh = hop;
        mi = i;
      }
    }
  }
  
  return mi;
}

void jellyfish_disconnect_graph(int lgraph[MAX_NODE][2*MAX_DEGREE], Path kpath[], int dstswitch, int size){
  int i, j, k;
  int s,d;
  //for(k=0; k< size;++k)jellyfish_print_path( -1, -1, kpath[k] );
  
  for(k=0; k< size;++k){
    for(i=0; kpath[k][i+1] != -1; ++i){
      s = kpath[k][i];
      d = kpath[k][i+1];
      
      //if(i >= 5)exit(0);
      for(j=0; j<MAX_PATH_LEN; ++j){
        if( lgraph[s][j] == d ){
          //printf("b4:lgraph[%d][%d]=%d\n",s,j,lgraph[s][j]);
          lgraph[s][j] = -1;
          //printf("AF;lgraph[%d][%d]=%d\n",s,j,lgraph[s][j]);
        }
        if( lgraph[d][j] == s ){
          //printf("b4:lgraph[%d][%d]=%d\n",d,j,lgraph[d][j]);
          lgraph[d][j] = -1;
          //printf("AF:lgraph[%d][%d]=%d\n",d,j,lgraph[d][j]);
        }
      }
    }
  }
}

void jellyfish_stitch_path(Path target, Path source, int snode){
  int i,j=1;
  for(i=0; i<MAX_DEGREE; ++i)if(target[i] == snode)break;
  for(++i;source[j] != -1;++i)target[i] = source[j++];
  for(;i<MAX_DEGREE;++i)target[i] = -1;
}

void jellyfish_copy_path(Path dpath, Path spath){
  int i;
  for(i=0; i<MAX_PATH_LEN;++i)dpath[i] = -1;
  for(i=0; spath[i] != -1;++i)dpath[i] = spath[i];
}

void jellyfish_make_path(int s, int d, Path spath){
  Path tpath;
  int i;
  
  jellyfish_init_path( tpath );
  tpath[0] = s;
  
  for(i=0; spath[i] != -1; ++i)tpath[i+1] = spath[i];
  tpath[i+1] = d;
  
  jellyfish_copy_path(spath, tpath );
}

int jellyfish_compute_K_Shortest_Path(int srcswitch, int dstswitch, Path kpath[], int K){
  Path allpath[MAX_JNUM_PATH];       // Shortest paths under consideration
  Path spath,                       // Shortest path
       tpath, 
       npath;                       // New path
  int kpathindex = 0,               // The current index to the K-paths
      allpathindex = 0,             // The urrent index to allpaths
      allpathsize = 0;              // The current size of the allpaths
  
  //int srcswitch, dstswitch;
  
  int i, j, pathexist;
  int intnode[MAX_PATH_LEN];
  
  // Initialize the paths 
  jellyfish_init_path( spath );
  jellyfish_init_path( tpath );
  jellyfish_init_path( npath );
  
  for(i=0; i<MAX_JNUM_PATH; ++i) jellyfish_init_path( allpath[i] );
  for(i=0; i<K; ++i)            jellyfish_init_path( kpath[i] );
  
  //srcswitch = s / jellyfish_p;
  //dstswitch = d / jellyfish_p;
  
  // Compute the first shortest path
  pathexist = jellyfish_shortestpath_routing(srcswitch,dstswitch,spath,graph);
  
  if( pathexist ){
    // The path exists between the s and d nodes.
  
    // Copy the first shortest path to the shortest paths under consideration
    jellyfish_copy_path(allpath[allpathindex++], spath);
    
    // Increase the size
    ++allpathsize;

    //printf("The first shortest path is :\n");
    //jellyfish_print_path(srcswitch,dstswitch, spath);
    
    // Loop until the K-paths have been found or there are paths under
    // consideration
    while((kpathindex < K) && (allpathsize > 0)){
    
      // Find the index of the shortest path in the allpath
      i = jellyfish_find_shortest_path(allpath);
      
      //printf("Allpath-size for k=%d is %d\n",kpathindex,allpathindex);
    
      //jellyfish_print_path_array(allpath, MAX_JNUM_PATH, srcswitch+totPE, dstswitch+totPE);
      //printf("for k=%d allpath index = %d is the shortest\n",kpathindex,i);
      //jellyfish_print_path(srcswitch, dstswitch, allpath[i] );
      
      // Initialize the path
      jellyfish_init_path( spath );
      
      // Copy the shortest path from allpath to spath for analysis
      jellyfish_copy_path( spath , allpath[ i ]);
      
      // Since this is the shortest amomg the paths under consideration,
      // add it to the K-paths
      jellyfish_copy_path( kpath[kpathindex++], spath);
      
      //printf("KPath[%d] constructed.\n",(kpathindex-1));
      //jellyfish_print_path(srcswitch+totPE, dstswitch+totPE, spath );
      
      // Remove from the allpaths
      //jellyfish_init_path( allpath[i] );
      
      // Need to make them contiguous
      for(; i< allpathsize; ++i)jellyfish_copy_path(allpath[i], allpath[i+1]);
      
      //Decrement the size
      --allpathsize;
      allpathindex = allpathsize;
      
      //printPathArray(kpath);
      //printf("The %d-th path in KPath entered.\n",kpathindex);
      //for(i=0; i<kpathindex; ++i) jellyfish_print_path(kpath[i]);
    
      //jellyfish_print_path(srcswitch, dstswitch, spath);  
      for(i=0; i<MAX_PATH_LEN; ++i) intnode[i] = -1;
      for(i=0; spath[i] != (dstswitch+totPE); ++i)  intnode[i] = spath[i];
    
      //for(j=0; intnode[j] != -1;++j)printf("IntNode[%d] = %d\n",j,intnode[j]);
      //printf("i = %d\n",i);
      // Analyze each link along the path between switches.
      
      for(j=i-1; j >= 0;j--){
        
        jellyfish_init_graph( lgraph );
        
        jellyfish_copy_graph( lgraph, graph);
        
        jellyfish_disconnect_graph(lgraph, kpath, dstswitch+totPE, kpathindex);
        //printGraph(lgraph);
        //jellyfish_print_graph(lgraph);
        //printf("Analysisng between srcsw = %d dstsw = %d\n",intnode[j],dstswitch+totPE);
        
        // Compuet the path if a path exists
        pathexist = jellyfish_shortestpath_routing(intnode[j]-totPE,dstswitch,tpath, lgraph);
        if( pathexist ){
          jellyfish_copy_path( npath, spath);
          jellyfish_stitch_path(npath, tpath, intnode[j]);
      
          //printf("Path computed for (%d=>%d):",intnode[j],dstswitch+totPE);
          //jellyfish_print_path( srcswitch+totPE, dstswitch+totPE, tpath );
          //printf("Stitched path:");
          //jellyfish_print_path( npath );
        
          jellyfish_copy_path(allpath[allpathindex++], npath);
          ++allpathsize;
        }
        else{
          //printf("*****No more path exists between the nodes %d and %d *****.\n",intnode[j],dstswitch+totPE);
        }
      }
    }
    
    //printf("Printing the K(%d)-paths:\n",K);
    //for(i=0; i<kpathindex; ++i) jellyfish_print_path(s, d, kpath[i]);
  }
  else{
    printf("*****No path exists between the nodes %d and %d *****.\n",srcswitch,dstswitch);
    exit(0);
  }
  
  return kpathindex;
}

void jellyfish_dump_topology_binary(){
  char topFile[100];
  FILE *fp = NULL;
  int i,j;
  
  struct TopologyEntry te;
  
  sprintf(topFile, "jellyfish_topology_%d_%d_%d.bin",jellyfish_nTOR, jellyfish_nPE, jellyfish_r);
  printf("Dumping the topology into the file %s\n",topFile);
  
  if( (fp = fopen(topFile, "wb")) == NULL ){
    printf("Could not open file %s for writing.\n",topFile);
    exit(0);
  }
  
  for(i=0; i<totNode; i++){
    for(j=0; graph[i][j]!=-1; j++){
      te.s = i;
      te.p = j;
      te.r = graph[i][j];
      te.b = bandwidth[i][j];
      fwrite(&te,sizeof(struct TopologyEntry),1,fp);
    }
  }
  
  if(fp != NULL)fclose(fp);
  
  printf("Dumping the topology into the file %s done.\n",topFile);
}

typedef struct PsuedoTreeQueueType{
  Node *element[MAX_PQ_SIZE];
  int head;
  int tail;
  int size;
}PsuedoTreeQueue;

void resetPTQ(PsuedoTreeQueue *q){
  q->head = 0;
  q->tail = 0;
  q->size = 0;
}

void printPTQ(PsuedoTreeQueue *q){
  for(;q->head < q->tail; q->head++) 
    printf("Queue[%d] = %d\n",q->head, q->element[q->head]->nid);
}

int pushPTQ(PsuedoTreeQueue *q, Node *x){
  if( q->size < MAX_PQ_SIZE ){
    q->element[q->tail++] = x;
    q->size++;
    //printf("Pushed value = %x, now head = %d, tail = %d,size = %d\n",
    //x,q->head,q->tail,q->size);
  }
  else{
    printf("Priority Queue has reached capacity. Capcity = %d current size = %d\n",MAX_PQ_SIZE,q->size);
    //printQ( q );
    
    exit(0);
  }
  return 0;
}

int popPTQ(PsuedoTreeQueue *q){
  q->head++;
  if( q->size > 0 ){
    q->size--;
    
    //printf("Popped now head = %d, tail = %d,size = %d\n",q->head,q->tail,q->size);
  }
  else{
    printf("Queue is already empty.\n");
    exit(0);
  }
  
  return 0;
}

int getSizePTQ(PsuedoTreeQueue *q){
  return q->size;
}

int isEmptyPTQ(PsuedoTreeQueue *q){
  return (q->size == 0);
}

Node* getFrontPTQ(PsuedoTreeQueue *q){
  return q->element[q->head];
}

//////////////////////////////////////////////////////

int getPathFromTree(Node *root, Node *node, Path path){
  int i=0;
  jellyfish_init_path( path );
  //printf("root node = %d node = %d equality = %d\n",root->nid,node->nid,(node == root));
  
  while( node != root ){
    path[i++] = node->nid;
    node = node->parent;
  }
  path[i++] = root->nid;
  return i;
}


void jellyfish_dump_paths(Node *root, Path paths[totSE * jellyfish_K]){
  int i=0,j,n;
  Node *tnode;
  static PsuedoTreeQueue q;
  Path path,tpath;
  resetPTQ( &q );
  pushPTQ( &q, root ) ;
  int index[totSE];
  
  for(i=0; i<totSE; ++i)index[i] = 0;
  for(i=0;i< totSE * jellyfish_K;++i) jellyfish_init_path( paths[i] );
  
  while( ! isEmptyPTQ(&q) ){
    tnode = getFrontPTQ(&q);
    //printf("tnode = %d\n",tnode->nid);
    n = getPathFromTree( root, tnode, tpath );
    //printf("Hop = %d ",n);
    //jellyfish_print_path( root->nid, -1, tpath );
    jellyfish_init_path( path );
    for(i=0;i<n; ++i) path[i] = tpath[n-i-1];
    //jellyfish_print_path( root->nid, -1, path );
    j = tpath[0]-totPE;
    i = j*jellyfish_K + index[j];
    index[j]++;
    
    //jellyfish_init_path( paths[i] );
    jellyfish_copy_path(paths[i], path);
    
    // pop out the front element
    popPTQ( &q );
  
    // for each node connected to it
    for(j=0; j < tnode->pindex; ++j){
      pushPTQ( &q, tnode->child[j] ) ;
      //printf("Pushed %d\n",tnode->child[j]->nid);
    }
  }
  //for(i=0;i<totSE*K;++i)
  //  jellyfish_print_path(root->nid,-1,paths[i]);
}

void cleanTree(Node *node){
  int i;
  for(i=0;i<MAX_DEGREE;++i){
    if( node->child[i] ){

      cleanTree( node->child[i] );
      
      free( node->child[i] );
      node->child[i] = NULL;
    }
  }
}

//int darray[MAX_NODE][300];
int darray[10][300];

void jellyfish_genkshortestpath(int s, Path paths[totSE * jellyfish_K]){
  Node *root = NULL;
  Node *temp = NULL;
  Node *node;
  static int count[MAX_NODE];  
  static   PriorityQueue pq;  
  int i,j,k,l,v;
  int Naa, Na, Nb, Nc, Nx, dis;
  PQElement pqe;
  PQElement tpqe;
  int found;  
  
  initPQElement( &pqe );
  
  for(i=0; i<MAX_NODE; ++i)for(j=0; j< jellyfish_K; ++j)darray[i][j] = LARGENUM;
      
  for(i=0; i<MAX_NODE; ++i)count[i] = 0;
  
  initPQ( &pq );
  
  // Deal with the root of the pseudo tree
  if( (root = (Node *)malloc(sizeof( Node ))) == NULL ){
    printf("Could not allocate memory for temp.\n");
    exit(0);
  }
  
  initTreeNode( root, s );
  
  pqe.snode = s;
  pqe.tnode = root;
  for(j=0; j < get_r(s); ++j){
    Nx = graph[s][j];
    pqe.e.s = s;
    pqe.e.d = Nx;
    pqe.dist = 1;
    pushPQ( &pq, pqe );
    count[Nx] = 1;
  }  
  
  while ( ! isQEmpty( &pq ) ){
    pqe = popPQ( &pq );
    Naa = pqe.snode; Naa = Naa;
    Na  = pqe.e.s; Na = Na;
    Nb  = pqe.e.d;
    dis = pqe.dist;
    
    // Concatenate
    if( (temp = (Node *)malloc(sizeof( Node ))) == NULL ){
      printf("Could not allocate memory for temp.\n");
      exit(0);
    }    
    initTreeNode( temp, Nb );
    temp->parent = pqe.tnode;
    pqe.tnode->child[ pqe.tnode->pindex++ ] = temp;
    
    // For all edges
    for(j=0; j < get_r(Nb); ++j){
      Nc = graph[Nb][j];
      
      found = FALSE;

      for(node = temp; node != NULL; node = node->parent){
        if(node->nid == Nc){
          found = TRUE;
          break;
        }
      }
      //printSubTree( root );
      
      k = jellyfish_K-1;
      if( (! found) && ( ( dis + 1 ) < darray[Nc][k]) ){
        if( count[Nc] < jellyfish_K ){
          ++count[Nc];
          
          initPQElement( &tpqe );
          tpqe.snode = Nb;
          tpqe.tnode = temp;
          tpqe.e.s = Nb;
          tpqe.e.d = Nc;
          tpqe.dist = dis+1;
          
          pushPQ( &pq, tpqe);
        }
        else{
          updatePQ( &pq, Nb, darray[Nc][k]);
        }
        darray[Nc][k] = dis + 1;
        v = darray[Nc][k];

        for(i=0;i<jellyfish_K;++i){          
          if(darray[Nc][i] > v){
            for(l=jellyfish_K-1;l>=i;l-- )darray[Nc][l] = darray[Nc][l-1];
            darray[Nc][i] = v;
            break;
          }
        }
      }
    }
  }
  
  jellyfish_dump_paths(root, paths );
  //for(i=0;i<totSE*K;++i)
  //  jellyfish_print_path(root->nid,i/K+totPE,paths[i]);
  // Do the cleanup
  cleanTree( root );
  free( root );
  root = NULL;
}


void jellyfish_dump_routing_ascii(){
  int i,j,k;
  char routefile[100];
  FILE *fp = NULL;  
  
  sprintf(routefile, "jellyfish_routing_%d_%d_%d.txt",
	  jellyfish_nTOR, jellyfish_nPE, jellyfish_r);
  
  if( (fp = fopen(routefile,"w")) == NULL ){
    printf("Could not open routefile %s\n for reading.",routefile);
    exit(0);
  }
  
  for(i=0;i<totSE;++i)
    for(j=0;j<totSE*jellyfish_K;++j){
      fprintf(fp, "(%d=>%d) = ", i+totPE,j/jellyfish_K + totPE );
      for(k=0;k<jellyfish_K;++k){
        if( pathcache[i][j][k] != -1)
          fprintf(fp, "%d ",pathcache[i][j][k]);
        //jellyfish_print_path( i+totPE, j/jellyfish_K + totPE, pathcache[i][j] );
      }
      fprintf(fp,"\n");
    }
  if(fp != NULL)fclose(fp);
}



void jellyfish_read_routing_to_memory(){
  int i=0,j=0;
  char routefile[100];
  FILE *fp = NULL;
  unsigned int stat;
  unsigned int nbytes;
  
  sprintf(routefile, "jellyfish_routing_%d_%d_%d.bin",
	  jellyfish_nTOR, jellyfish_nPE, jellyfish_r);
  
  printf("Reading the routing info from %s to memory.\n",routefile);
                                     
  if( (fp = fopen(routefile,"rb")) == NULL ){
    printf("Could not open routefile %s\n for reading.",routefile);
    exit(0);
  }
  nbytes = jellyfish_K * sizeof(Path);
  while( (stat = fread(pathcache[i][j*jellyfish_K], nbytes, 1, fp)) == 1 ){
    j++;
    
    if(j == totSE ){
      i++;
      j = 0;
    }
  }
  
  if(fp != NULL)fclose(fp);
  
  //for(i=0;i<totSE;++i)
  //  for(j=0;j<totSE*jellyfish_K;++j)
  //    jellyfish_print_path( i+totPE, j/jellyfish_K + totPE, pathcache[i][j] );
  
  jellyfish_dump_routing_ascii();
  
  printf("Reading the routing info from %s to memory done.\n",routefile);
}



void print_shortpathcache(int i){
  int j,k,l;
  for(j=0; j < totSE; j++){
    printf("Src[%d]==>Dst[%d]\n",i+totPE,j+totPE);
    for(k=0; k<JFISH_MAX_NUM_PATH; k++){
      for(l=0; 
      (l < MAX_PATH_LEN && shortpathcache[i][j*JFISH_MAX_NUM_PATH+k][l] != -1); 
      l++){
        printf("%d ",(int)shortpathcache[i][j*JFISH_MAX_NUM_PATH+k][l]);
      }
      printf("\n");
    }
  }
}

void jellyfish_read_routing_to_memory1(){
  int i=0,j=0;
  char routefile[100];
  FILE *fp = NULL;
  unsigned int stat;
  unsigned int nbytes;

  sprintf(routefile, "jellyfish_routing_%d_%d_%d.bin",
	  jellyfish_nTOR, totPE, jellyfish_r);
  
  printf("Reading the routing info from %s to memory.\n",routefile);
                                     
  if( (fp = fopen(routefile,"rb")) == NULL ){
    printf("Could not open routefile %s\n for reading.",routefile);
    exit(0);
  }
  nbytes = JFISH_MAX_NUM_PATH * sizeof(ShortPath);
  while( (stat = fread(shortpathcache[i][j*JFISH_MAX_NUM_PATH], nbytes, 1, fp)) == 1 ){
    //    print_shortpathcache(i);
    j++;
    
    if(j == totSE ){
      i++;
      j = 0;
    }
  }
  
  
  if (i != jellyfish_nTOR) {
    printf("Read routes from file failed.\n");
    exit(0);
  }
  
  
  if(fp != NULL)fclose(fp);
  
  //  jellyfish_dump_routing_ascii();
  
  printf("Reading the routing info from %s to memory done.\n",routefile);
}


//void jellyfish_readkpath(int src, int dst, Path kpath[ jellyfish_K ]){
void jellyfish_readkpath(int src, int dst, Path kpath[ ]){
  unsigned int nbytes;
  int srcswitch = getswid(src);// / jellyfish_p;
  int dstswitch = getswid(dst);// / jellyfish_p;
  int i,j;
  
  for(i=0;i<jellyfish_K;++i)jellyfish_init_path( kpath[i] );
  if( src == dst ){
    kpath[0][0] = src;
    return;
  }
  if(srcswitch == dstswitch){
    kpath[0][0] = src;
    kpath[0][1] = srcswitch+totPE;
    kpath[0][2] = dst;
    
    return;
  }
  
  nbytes = jellyfish_K * sizeof(Path);
  memcpy(kpath, pathcache[srcswitch][dstswitch*jellyfish_K], nbytes );
  
  for(i=0;i<jellyfish_K;++i){
    for(j=0; j<MAX_PATH_LEN;++j)if(kpath[i][j] == -1)break;
      
    kpath[i][j+1] = dst;
    for(; j>0; --j)kpath[i][j] = kpath[i][j-1];
    kpath[i][0] = src;
  }
  
  //  for(i=0; i<jellyfish_K;++i) jellyfish_print_path( src, dst, kpath[i] );
}


//void jellyfish_readnewkpath(int src, int dst, Path kpath[ jellyfish_K ]){
void jellyfish_readnewkpath(int src, int dst, Path kpath[ ]){
  int srcswitch = getswid(src); // / jellyfish_p;
  int dstswitch = getswid(dst); // / jellyfish_p;
  int i,j;
  int ld = getlocalid(dst); //dst % jellyfish_p;
  int npath;
  int nshortpath=0;

  //printf("src: %d, srcswitch: %d", src, srcswitch);
  //printf("dst: %d, dstswitch: %d", dst, dstswitch);
  //  printf("ssw = %d dsw = %d, jellyfish_K %d, %d,  %d, ld =%d\n", srcswitch, dstswitch,
  //	 jellyfish_K, NEWR_MAX_HOP, _REAL_PATH_NUM, ld);

  for(i=0;i<jellyfish_K;++i)jellyfish_init_path( kpath[i] );

  if( src == dst ){
    kpath[0][0] = src;
    return;
  }
  if(srcswitch == dstswitch){
    kpath[0][0] = src;
    kpath[0][1] = srcswitch+totPE;
    kpath[0][2] = dst;

    return;
  }

/*  if( topology != JELLYFLY ){
    if ((NEWR_MAX_HOP != 2) &&  (NEWR_MAX_HOP != 3)) {
      printf("new jellyfish routing, NEWR_MAX_HOP should be 2 (small) or 3 (large).\n");
      exit(0);
    }
  }*/

  npath = 0;
  for (i=0; i<JFISH_MAX_NUM_PATH; i++) {
    if (shortpathcache[srcswitch][dstswitch*JFISH_MAX_NUM_PATH+i][NEWR_MAX_HOP] == -1) nshortpath++;
    if (shortpathcache[srcswitch][dstswitch*JFISH_MAX_NUM_PATH+i][NEWR_MAX_HOP+1] != -1)
      break;
    if (shortpathcache[srcswitch][dstswitch*JFISH_MAX_NUM_PATH+i][0] == -1)
      break;
  }
  npath = i;

  //  printf("npath = %d\n", npath);

  /*
  if ((npath < 1) || (npath > 100)) {
    printf("Strange, npath = %d, src = %d, dst = %d, JFISH_M_PATH = %d\n",
	   npath, srcswitch, dstswitch, JFISH_MAX_NUM_PATH);
    npath = 0;
    for (srcswitch = 0; srcswitch < jellyfish_nTOR; srcswitch++)
      for (dstswitch = 0; dstswitch < jellyfish_nTOR; dstswitch++) {
	for (i=0; i<JFISH_MAX_NUM_PATH; i++) {
	  if (shortpathcache[srcswitch][dstswitch*JFISH_MAX_NUM_PATH+i][NEWR_MAX_HOP+1] != -1)
	    break;
	}
	npath = i;
	printf("(%d->%d) npath = %d\n", srcswitch, dstswitch, npath);
      }
    exit(0);
  }
  */

  static int SD_counter=0;
  //jellyfish_ths = jellyfish_K;
  //printf("max hop %d, number of paths: %d, number of short paths %d\n", NEWR_MAX_HOP, npath,nshortpath);
  if (npath < jellyfish_K ){
    //printf("number of paths between %d and %d is: %d, updated to %d, update count %d\n",src,dst,npath,jellyfish_K, SD_counter++);
    npath = jellyfish_K; // no short paths, pick first 8 paths
    jellyfish_kpath_routing(src, dst, kpath, jellyfish_K);

  }
  else if(nshortpath>=jellyfish_K){
    //if there are sufficient short paths available, then use KPATH
    jellyfish_kpath_routing(src, dst, kpath, jellyfish_K);

    
  }else{

    /*while (((shortpathcache[srcswitch][dstswitch*JFISH_MAX_NUM_PATH+npath-1][0] 
      	   == -1)) && (npath >0)) npath --;
    if (npath == 0) {
      printf("jellyfish_read_newkpath: no path from %d to %d\n", srcswitch,
	   dstswitch);
      exit(0);
    }*/


    for (i=0; i<jellyfish_K; i++) {
      kpath[i][0] = src;

      if (routing == JELLYFISH_NEWKPATH_ROUTING) {
        for (j=0; shortpathcache[srcswitch][dstswitch*JFISH_MAX_NUM_PATH+(ld*jellyfish_K+i) % npath][j] != -1; j++) {
          kpath[i][j+1] =   shortpathcache[srcswitch][dstswitch*JFISH_MAX_NUM_PATH+(ld*jellyfish_K+i) % npath][j];
        }
        kpath[i][j+1] = dst;
      } else if (routing == JELLYFISH_KPATH_ROUTING) {
        for (j=0; shortpathcache[srcswitch][dstswitch*JFISH_MAX_NUM_PATH+i][j] != -1; j++) {
	  kpath[i][j+1] =   shortpathcache[srcswitch][dstswitch*JFISH_MAX_NUM_PATH+i][j];
        }
        kpath[i][j+1] = dst;
      }
      else if (routing == JELLYFISH_KPATH_LLSKR_MIXED_ROUTING) {
        if(i<nshortpath){
          for (j=0; shortpathcache[srcswitch][dstswitch*JFISH_MAX_NUM_PATH+i][j] != -1; j++) {
	    kpath[i][j+1] =   shortpathcache[srcswitch][dstswitch*JFISH_MAX_NUM_PATH+i][j];
          }
          kpath[i][j+1] = dst;
        }else{
          for (j=0; shortpathcache[srcswitch][dstswitch*JFISH_MAX_NUM_PATH+(ld*jellyfish_K+i) % (npath-nshortpath)][j] != -1; j++) {
            kpath[i][j+1] =   shortpathcache[srcswitch][dstswitch*JFISH_MAX_NUM_PATH+(ld*jellyfish_K+i) % (npath-nshortpath)][j];
          }
          kpath[i][j+1] = dst;        
        }
      }//endif

    }//i loop
  }

}

/* min_routing */
/* this assume the switch path is alreadys in variable path */

void jellyfish_min_routing(int src, int dst, Path path){
  //Path kpath[jellyfish_K];
  int i;
  int srcswitch = getswid(src); // / jellyfish_p;
  int dstswitch = getswid(dst); // / jellyfish_p;
  //  unsigned int nbytes;
  
  jellyfish_init_path(path);
  if( src == dst ){
    path[0] = src;
    path[1] = -1;
    return;
  }
  if(srcswitch == dstswitch){
    path[0] = src;
    path[1] = srcswitch+totPE;
    path[2] = dst;
    path[3] = -1;
    return;
  }

  //  nbytes = sizeof(Path);
  /*
  memcpy(path, pathcache[srcswitch][dstswitch*jellyfish_K], nbytes );
  */

  for(i=0;i<MAX_PATH_LEN;++i)if(path[i] == -1)break;
  
  path[i+1] = dst;
  for(;i>0;--i)path[i] = path[i-1];
  path[0] = src;
  
  //jellyfish_print_path(src, dst,path);
}

// Computes the random 1 path routing.
// It reads/computes the K paths and randomly chooses one path from them
void jellyfish_random2_routing(int src, int dst, Path path){
  int i,j;
  int srcswitch = getswid(src); // / jellyfish_p;
  int dstswitch = getswid(dst); // / jellyfish_p;
  
  jellyfish_init_path(path);
  if( src == dst ){
    path[0] = src;
    return;
  }
  if(srcswitch == dstswitch){
    path[0] = src;
    path[1] = srcswitch+totPE;
    path[2] = dst;
    return;
  }
  
  Path kpath[ jellyfish_K ];
  
  jellyfish_readkpath(src, dst, kpath);
  
  //for(i=0; i< jellyfish_K; ++i) jellyfish_print_path(src, dst, kpath[i]);

  for(i=0; i< jellyfish_K; ++i)if(kpath[i][0] == -1)break;
  
  // There are i paths for this SD pair
  j = rand() % i;
  memcpy(path,kpath[j],sizeof(Path));
  
  //printf("Chosen path is the %d-th path.\n",j);
  //jellyfish_print_path(src, dst, path);
}


// Computes the random 1 path routing.
// It reads/computes the K paths and randomly chooses one path from them
void jellyfish_random1_routing(int src, int dst, Path path){
  int i,j;
  int ii;

  int srcswitch = getswid(src); // / jellyfish_p;
  int dstswitch = getswid(dst); // / jellyfish_p;

  jellyfish_init_path(path);
  if( src == dst ){
    path[0] = src;
    path[1] = -1;
    return;
  }
  if(srcswitch == dstswitch){
    path[0] = src;
    path[1] = srcswitch+totPE;
    path[2] = dst;
    path[3] = -1;
    return;
  }

  path[0] = src;
  path[1] = srcswitch+totPE;
  path[2] = graph[path[1]][rand()%get_r(path[1])];
  
  Path kpath[ jellyfish_K ];
  
  //if (path[2] - totPE < jellyfish_x) 
    //ii = (path[2]-totPE)*(jellyfish_p+1);
  //else 
    //ii = jellyfish_x * (jellyfish_p+1)+
    //  (path[2] -totPE-jellyfish_x)*jellyfish_p;
  ii = get_firstPE(path[2]);

  jellyfish_readkpath(ii, dst, kpath);

  for(i=0; i< jellyfish_K; ++i)if(kpath[i][0] == -1)break;
  
  // There are i paths for this SD pair
  //  j = rand() % i;
  j = 0; // always shortest path
  //for(i=0; i< jellyfish_K; ++i) jellyfish_print_path(src, dst, kpath[i]);

  for (i=2; kpath[j][i] != -1; i++) path[i+1] = kpath[j][i];    
  path[i+1] = -1;

  //  jellyfish_print_path(src, dst, path);
}

void jellyfish_readkpath1(int src, int dst, Path kpath[jellyfish_K])
{
  Path path;
  int i;
  int srcswitch = getswid(src); // / jellyfish_p;
  int dstswitch = getswid(dst); // / jellyfish_p;

  //  printf("inside readkpath %d -> %d\n", src, dst);

  jellyfish_init_path(path);
  for (i=0; i<jellyfish_K; i++) jellyfish_init_path(kpath[i]);

  if( src == dst ){  // one node
    kpath[0][0] = src;
    return;
  }

  if(srcswitch == dstswitch){ // in one switch
    kpath[0][0] = src;
    kpath[0][1] = srcswitch+totPE;
    kpath[0][2] = dst;
    return;
  }

  for (i=0; i<jellyfish_K; i++) {
    jellyfish_random1_routing(src, dst, kpath[i]);
  }
  //  for(i=0; i<jellyfish_K;++i) jellyfish_print_path( src, dst, kpath[i] );
}

#define MAX_ONE_SWSOURCE_PATH (39*39*39*39+1)
short oneswsourcepath[MAX_ONE_SWSOURCE_PATH][6];
int  swpcount = 0;

#define isPE(x) (((x)>=0) && ((x) < totPE))
#define isSW(x) (((x)>=totPE) && ((x) < totPE+jellyfish_nTOR)) 

void jelly_compute_oneswsourcepath(int srcsw, int hop)
{
  int i;
  int j, k, l;

  int s, e;
  int nsw;
  if (!isSW(srcsw)){
    printf("%d is not a switch.\n", srcsw);
    exit(0);
  }

  if (hop > 5) {
    printf("MAX_hop count for path is limited to 5 at this time.\n");
    exit(0);
  }

  for (i=0; i<MAX_ONE_SWSOURCE_PATH; i++) 
    for (j=0; j<6; j++) oneswsourcepath[i][j] = -1;
  
  oneswsourcepath[0][0] = srcsw;
  oneswsourcepath[0][1] = -1;
  swpcount = 1;

  s = -1; e = 0;
  for (i=1; i<= hop; i++) {
    s = e;  
    e = swpcount;
    for (k=0; k<jellyfish_nPorts; k++) {
      for (j=s; j<e; j++) {
	nsw = graph[oneswsourcepath[j][i-1]][k];
        if (!isSW(nsw)) continue;
        for (l=0; l<i; l++)
	  if (nsw == oneswsourcepath[j][l]) break;
        if (l != i) continue; // find a new path
        for (l=0; l<i; l++)
	  oneswsourcepath[swpcount][l] = oneswsourcepath[j][l];
	oneswsourcepath[swpcount][l] = nsw;
	swpcount++;
        if (swpcount >= MAX_ONE_SWSOURCE_PATH) {
	  printf("Run out of oneswsourcepath memory. %d\n", swpcount);
	  exit(0);
	}
      }
    }
    //    printf("s = %d, e = %d, swpcount = %d\n", s, e, swpcount);
  }


  //  printf("src = %d, swpcount = %d\n", srcsw, swpcount);

  
  /*jellyfish_print_topology();
  for (i=0; i<swpcount+2; i++) {
    for (j=0; j<6; j++) {
      printf("%d ", oneswsourcepath[i][j]);
    }
    printf("\n");
  }
  exit(0);
  */
}

void packswpaths(int srcsw, int K)
{
  int count[MAX_NODE];
  int i, j, k;
  int lsw;

  if (!isSW(srcsw)) {
    printf("%d is not a switch.\n", srcsw);
    exit(0);
  }

  for (i=0; i<MAX_NODE; i++) count[i] = 0;
  for (i=1; i<swpcount; i++) {
    if (oneswsourcepath[i][0] != srcsw) {
      printf("Path is wrong %d %d\n", oneswsourcepath[i][0], srcsw);
      exit(0);
    }
    for (j=0; (oneswsourcepath[i][j+1] != -1) && (j < 6); j++);
    if (j==6) {
      printf("Shortpath too long.\n");
      exit(0);
    }
    lsw = oneswsourcepath[i][j];
    if (count[lsw] >= K) continue;
    for (j=0; oneswsourcepath[i][j] != -1; j++) {
      shortpathcache[srcsw-totPE][(lsw - totPE)*K + count[lsw]][j] = 
	oneswsourcepath[i][j];
    }
    if (shortpathcache[srcsw-totPE][(lsw - totPE)*K + count[lsw]][j-1] != lsw) {
      printf("Something is wrong. %d %d\n",
	     shortpathcache[srcsw-totPE][(lsw - totPE)*K + count[lsw]][j-1], lsw
	     );

      for (j=0; j<6; j++) {
	printf("%d\n",
	       shortpathcache[srcsw-totPE][(lsw - totPE)*K + count[lsw]][j]);
      }
    }

    for (k=j; k<MAX_SHORT_PATH_LEN -1; k++)
      shortpathcache[srcsw-totPE][(lsw - totPE)*K + count[lsw]][k] = -1;

    count[lsw]++;
  }

  for (i=totPE; i<totPE+jellyfish_nTOR; i++) {
    //  printf("count[%d] = %d\n", i, count[i]);
    if (0 && (count[i] < 10) && (i != srcsw)) {
      printf("Strange path distribution, consider changing topology.\n");
      exit(0);
    }
  }

  /*  
  for (i=0; i<jellyfish_nTOR; i++) {
    printf("\nsrc = %d, dst = %d, k = %d\n", srcsw, i+totPE, K);
    for (j=0; j<K; j++) {
      for (k=0; k<MAX_SHORT_PATH_LEN-1; k++) {
	printf("%d ", shortpathcache[srcsw-totPE][i*K+j][k]);
      }
      printf("\n");
    }
    }*/
  
}



//Path allpath[MAX_ONE_SWSOURCE_PATH];
        
void jellyfish_dump_routing_binary(){
  int i;
  //Path kpath[K];
  char routeFile[100];
  FILE *fp = NULL;
  unsigned int stat;
  unsigned int nbytes;
  
  sprintf(routeFile, "jellyfish_routing_%d_%d_%d.bin",jellyfish_nTOR, jellyfish_nPE, 
	  jellyfish_r);
  
  printf("Dumping the routing into the file %s\n",routeFile);
  
  if( (fp = fopen(routeFile, "wb")) == NULL ){
    printf("Could not open file %s for writing.\n",routeFile);
    exit(0);
  }
  /*
    if (routing == JELLYFISH_KPATH_ROUTING) {
    for(i=0; i<totSE; ++i){
      jellyfish_genkshortestpath(i+totPE, allpath );
      nbytes = totSE * jellyfish_K * sizeof(Path);
      if ( (stat = fwrite( allpath, nbytes, 1, fp )) == 0){
	printf("Could not write to file %s.\n",routeFile);
	exit(0);
      }
      else{      
	printf("Switch %d done, Wrote %d chunk of %d bytes to file \n", 
	       i,stat,nbytes);
      }    
    }    
    //printf("Write done for srcswitch = %d\n",i);
    } else if (routing == JELLYFISH_NEWKPATH_ROUTING)
  */ 
  {
    for(i=0; i<totSE; ++i){
      /* memset(&shortpathcache[i],-1,sizeof(shortpathcache[i])), */

      jelly_compute_oneswsourcepath(i+totPE, 4);
      packswpaths(i+totPE, JFISH_MAX_NUM_PATH);
      nbytes = totSE * JFISH_MAX_NUM_PATH * sizeof(ShortPath);
      //      print_shortpathcache(i);      

      if ((stat = fwrite(shortpathcache[i], nbytes, 1, fp )) != 1){
	printf("Could not write to file %s.\n",routeFile);
	exit(0);
      }
      else{
	printf("Switch %d done, Wrote %d chunk of %d bytes to file \n", 
	       i,stat,nbytes);
	       //print_shortpathcache(i);
      }    
    }
  }    
    //printf("Write done for srcswitch = %d\n",i);
  
  printf("Dumping the routing into the file %s done.\n",routeFile);
  
  if(fp != NULL)fclose(fp);
}


int pathIndex2;
qType q2;

/*---------------------------------------------*/


/* this function traverses the network recursively using DFS to find all paths between src and dest */
int jfind_all_paths(int srcNode, int srcSW, int dstNode,int dstSW,  int allpath[JFISH_MAX_NUM_PATH][MAX_PATH_LEN], int allhop[MAX_PATH_LEN], int maxhop)
{

  //printf("jfind_all_paths: current params: srcnode %d srcSW %d dstNode %d dstSW %d maxhop %d pathIndex2 %d\n", srcNode, srcSW, dstNode, dstSW, maxhop, pathIndex2);
  int i, j;
  stackElementType curNode;
  int curNode_r;

  curNode.node = srcSW;
  curNode.path[0] = srcSW;
  curNode.pid = 1;
  curNode_r = get_r(curNode.node);
  
  QPush(&q2, curNode);
  while(QIsEmpty(&q2) == 0 && pathIndex2 < MAX_ALLPATH){
    curNode = QPop(&q2);
    for(i=0;i<curNode_r;i++){
      if(graph[curNode.node][i] == dstSW){
	for(j=0;j<curNode.pid;j++){
	  allpath[pathIndex2][j+1] = curNode.path[j];
	}
	allpath[pathIndex2][0] = srcNode;
	allpath[pathIndex2][j + 1] = dstSW;
	allpath[pathIndex2][j + 2] = dstNode;
	allpath[pathIndex2][j + 3] = -1;
	pathIndex2++;
	allhop[curNode.pid]++;
	//printf("%d,%d, %d\n",pathIndex, curNode.pid, allpath[pathIndex-1][j + 3]);
      }
      else if(curNode.pid < maxhop){
	stackElementType tmpNode;
	tmpNode.node = graph[curNode.node][i];
	for(j=0;j<curNode.pid;j++){
	  if(tmpNode.node != curNode.path[j])
	    tmpNode.path[j] = curNode.path[j];
	  else break;
	}
	if(j == curNode.pid){
	  tmpNode.path[j] = graph[curNode.node][i];
	  tmpNode.pid = j+1;
	  QPush(&q2, tmpNode);
	}
      }
    }
  }

  return pathIndex2;
  
}

/* this function finds k shortest paths between src and dest */
int jfind_k_paths(int srcNode, int srcSW, int dstNode,int dstSW,  int kpath[JFISH_MAX_NUM_PATH][MAX_PATH_LEN])
{
  int i, j;
  stackElementType curNode;
  int curNode_r;

  curNode.node = srcSW;
  curNode.path[0] = srcSW;
  curNode.pid = 1;
  curNode_r = get_r(curNode.node);

  QPush(&q2, curNode);
  while(QIsEmpty(&q2) == 0 && pathIndex2 < jellyfish_K){
    curNode = QPop(&q2);
    for(i=0;i<curNode_r;i++){
      if(graph[curNode.node][i] == dstSW){
	for(j=0;j<curNode.pid;j++){
	  kpath[pathIndex2][j+1] = curNode.path[j];
	}
	kpath[pathIndex2][0] = srcNode;
	kpath[pathIndex2][j + 1] = dstSW;
	kpath[pathIndex2][j + 2] = dstNode;
	kpath[pathIndex2][j + 3] = -1;
	pathIndex2++;
      }
      else{
	stackElementType tmpNode;
	tmpNode.node = graph[curNode.node][i];
	for(j=0;j<curNode.pid;j++){
	  if(tmpNode.node != curNode.path[j])
	    tmpNode.path[j] = curNode.path[j];
	  else break;
	}
	if(j == curNode.pid){
	  tmpNode.path[j] = graph[curNode.node][i];
	  tmpNode.pid = j+1;
	  QPush(&q2, tmpNode);
	}
      }
    }
  }

  return pathIndex2;
  
}


/* this routing finds all paths between two nodes */
int jellyfish_allpath_routing(int src, int dst, 
			      int allpath[JFISH_MAX_NUM_PATH][MAX_PATH_LEN], int allhop[MAX_PATH_LEN], int maxhop)
{
  /*  int i, j;
  int ii, jj, jjj;

  int srcg, dstg;
  int srcsw1, dstsw1;
  int srcsw2, dstsw2;
  int k, kk;*/
  //printf("jellyfish_allpath_routing_2511: src %d dst %d maxhop %d\n", src,dst,maxhop);
  int srcsw, dstsw;
  //int visited[totSE];

  if ((src < 0) || (src >= totPE)) {
    printf("src %d is not a processing node.\n", src);
    exit(0);
  }

  if ((dst < 0) || (dst >= totPE)) {
    printf("dst %d is not a processing node.\n", dst);
    exit(0);
  }

  //srcsw = src / jellyfish_p + totPE;
  //dstsw = dst / jellyfish_p + totPE;
  srcsw = getswid(src)+totPE;
  dstsw = getswid(dst)+totPE;


  if (src == dst) { /* same node*/
    allpath[0][0] = src;
    allpath[0][1] = -1;
    return 0;
  } else if (srcsw == dstsw) {/* same switch */
    allpath[0][0] = src;
    allpath[0][1] = srcsw;
    allpath[0][2] = dst;
    allpath[0][3] = -1;
    return 1;
  } else {
    QInit(&q2);
    pathIndex2 = 0;
    return jfind_all_paths(src, srcsw, dst, dstsw, allpath, allhop, maxhop);
  }
}

/* this routing finds k shortest paths between two nodes */
int jellyfish_kpath_routing(int src, int dst, 
			      int kpath[JFISH_MAX_NUM_PATH][MAX_PATH_LEN], int k)
{
  int srcsw, dstsw;

  if ((src < 0) || (src >= totPE)) {
    printf("src %d is not a processing node.\n", src);
    exit(0);
  }

  if ((dst < 0) || (dst >= totPE)) {
    printf("dst %d is not a processing node.\n", dst);
    exit(0);
  }

  srcsw = getswid(src) + totPE;
  dstsw = getswid(dst) + totPE;


  if (src == dst) { /* same node*/
    kpath[0][0] = src;
    kpath[0][1] = -1;
    return 0;
  } else if (srcsw == dstsw) {/* same switch */
    kpath[0][0] = src;
    kpath[0][1] = srcsw;
    kpath[0][2] = dst;
    kpath[0][3] = -1;
    return 1;
  } else {
    QInit(&q2);
    pathIndex2 = 0;
    return jfind_k_paths(src, srcsw, dst, dstsw, kpath);
  }
}

int jcompareHistogram_data (const void *a, const void *b)
{
  const struct Histogram_data *ia = (struct Histogram_data *)a; // casting pointer
  const struct Histogram_data *ib = (struct Histogram_data *)b; // casting pointer
  if ( ia->value <  ib->value ) return -1;
  else if ( ia->value == ib->value ) return 0;
  else return 1;
}


void jdump_histogram_stat(char *filename){
  Path kpath[JFISH_MAX_NUM_PATH];
  const int HISTOGRAM_SIZE=10000;
  struct Histogram_data link_histogram[HISTOGRAM_SIZE];
  //struct Histogram_data path_histogram[HISTOGRAM_SIZE];
  int max_link_usage;
  int *path_histogram;
  int link_histogram_size=0;
  int path_cost=0, max_path_cost=0;
  int num_path;
  FILE *fd, *ofd1, *ofd2;
  char buf[1000], *ch;
  int i, j, k, ii;
  long long int tot_len;
  int next_s, next_d, prev_s;
  long long int count = 0;
  int totpath=0, tothop=0;
  long long int flid = 0;

  for(i=0;i<HISTOGRAM_SIZE;i++){
    link_histogram[i].value=0;
    link_histogram[i].frequency=0;
  }
  if ((fd = fopen(filename, "r")) == NULL) {
    printf("file %s does not exist.\n", filename);
    exit(0);
  }
  if ((ofd1 = fopen("link_stat", "w")) == NULL) {
    printf("file creation failed\n");
    exit(0);
  }
  if ((ofd2 = fopen("path_stat", "w")) == NULL) {
    printf("file creation failed\n");
    exit(0);
  }
  if( routingType == SINGLEPATH && _REAL_PATH_NUM > 1){
    printf("Single-path routing requested but multi-path routing found.\n");
    printf("Please recompile with the correct number of paths.\n");
    //writeExceptionToModelResult( ofd1 );
    exit(0);
  }
  //_____________________________________________________________________PASS 1
  tot_len = 0;
  prev_s = -1;
  ch = fgets(buf, 1000, fd);
  i = sscanf(buf, "%d %d", &next_s, &next_d);
  
  while (next_s > -1) {
    if (next_s == next_d || graph[next_s][0] == graph[next_d][0]) {
      ch = fgets(buf, 1000, fd);
      i = sscanf(buf, "%d %d", &next_s, &next_d);
      continue;
    }
    if (next_s != next_d) {
      count++;
      if (count % 1000000 == 0)
        printf("count = %lld\n", count);
    }
    if( next_s >= MAX_NPROCS || next_d >= MAX_NPROCS ){
      printf("Source ID(%d) and/or Destination ID(%d) > MAX_NPROCS(%d)\n", 
	     next_s,next_d,MAX_NPROCS);
      //writeExceptionToModelResult( ofd1 );
      exit(0);
    }
    num_path = jellyfish_kpath_routing(next_s, next_d, kpath, 60);
    //num_path = jellyfish_K;

    for (i=0, ii = 0; ii < num_path; i++, ii++) {   // for each of the K paths
      if (kpath[ii][0] == -1) ii = 0;
      if (kpath[ii][0] != next_s) {
	//printf("%d,%d,%d,%d,%d\n",kpath[0][0],kpath[0][1],kpath[0][2],next_s, next_d);
        printf("Path computation gone wrong. The first hop should be\
        the source node %d but is %d\n",next_s,kpath[ii][0]);
        exit(0);
      }
      for (j=0; kpath[ii][j+1] != -1; j++) {    //for each hop of the current path
        for (k=0;(k < MAX_DEGREE) && (graph[kpath[ii][j]][k] != kpath[ii][j+1]); k++);
        if (k >= MAX_DEGREE) {
          printf("Something is wrong XX3, ii = %d j = %d,%d,%d.\n", ii, j, kpath[ii][j], kpath[ii][j+1]);
          exit(0);
        }
        graph_m[kpath[ii][j]][k] = 1;
        graph_m1[kpath[ii][j]][k] ++;
      }
      for (j=0; kpath[ii][j+1] != -1; j++, tot_len++);
    }
    prev_s = next_s;
    ch = fgets(buf, 1000, fd);
    i = sscanf(buf, "%d %d", &next_s, &next_d);
    flid++;
  }
  printf("PASS 1 COMPLETE\n");
  
  for (i=totPE; i<totNode; i++) {               // for every switch
    for (j=0; graph[i][j] != -1; j++) {         //exhaustive search on all links of sel$
      for(k=0;k<link_histogram_size;k++){
        if(graph[i][j]<totPE) break;
        else if(link_histogram[k].value==graph_m1[i][j]){       //count freq of switch-$
          link_histogram[k].frequency++;
          break;
        }
      }
      if(k==link_histogram_size){       //new "link usage" value...added at the end
        link_histogram[k].value = graph_m1[i][j];
        //path_histogram[k].va+lue = graph_m1[i][j];            //required in pass 2
        link_histogram[k].frequency = 1;
        link_histogram_size++;
      }
    }
  }
  int totfreq=0;
  double avgload, totload=0;
  qsort(link_histogram, link_histogram_size, sizeof (struct Histogram_data), jcompareHistogram_data);
  for(i=0;i<link_histogram_size;i++){
    //fprintf(ofd1, "%d %d\n",link_histogram[i].value, link_histogram[i].frequency);
    if(i == link_histogram_size - 1)
      fprintf(ofd1, "%d ",link_histogram[i].value);
    totload += (link_histogram[i].value * link_histogram[i].frequency);
    totfreq += link_histogram[i].frequency; 
  }
  avgload = totload / totfreq;
  fprintf(ofd1, "%3.2f %d\n",avgload, jellyfish_K);
  if( ofd1 )fclose(ofd1);
  /*
  max_link_usage = link_histogram[link_histogram_size-1].value;
  path_histogram = (int *)malloc(max_link_usage*sizeof(int));
  memset (path_histogram, 0, sizeof(int) * max_link_usage);
  fseek(fd, 0, SEEK_SET);
  prev_s = -1;
  ch = fgets(buf, 1000, fd);
  i = sscanf(buf, "%d %d", &next_s, &next_d);
  flid = 0;
  while (next_s > -1)
    {
      if (next_s == next_d || graph[next_s][0] == graph[next_d][0]) {
	ch = fgets(buf, 1000, fd);
	i = sscanf(buf, "%d %d", &next_s, &next_d);
	continue;
      }
      
      jellyfish_kpath_routing(next_s, next_d, kpath, 60);
      num_path = jellyfish_K;      
      for (i=0, ii=0; ii < num_path; i++, ii++)
	{
	  
	  max_path_cost=0;
	  if (kpath[ii][0] == -1) printf("%d,%d,%d\n",next_s,next_d,ii);//ii = 0;
	  for (j=1; kpath[ii][j+1] != -1; j++)      //skip j=0 as it  indicates the edge between PE to SE
	    {
	      for (k=0; (k<MAX_DEGREE) && graph[kpath[ii][j]][k] != kpath[ii][j+1]; k++);
	      if (j >= MAX_DEGREE) {
		printf("Path computation gone wrong for (%d,%d).",next_s,next_d);
		exit(0);
	      }
	      if(graph[kpath[ii][j]][k]>=totPE )
		path_cost = graph_m1[kpath[ii][j]][k];
	      else
		path_cost=0;
	      if(path_cost>max_path_cost && path_cost<=max_link_usage) max_path_cost = path_cost;
	    }
	  //path length calculation
	  tothop+=(j-2);//need further testing on this!
	  totpath++;
	  path_histogram[max_path_cost]++;
	}
      prev_s = next_s;
      ch = fgets(buf, 1000, fd);
      i = sscanf(buf, "%d %d", &next_s, &next_d);
      flid++;
    }
  printf("Items scanned on pass 2: %d\n", (int)flid);
  
  for(i=1;i<=max_link_usage;i++)
    if(path_histogram[i]>0)
      fprintf(ofd2, "%d %d\n",i, path_histogram[i]);
  fprintf(ofd2, "Average path length in this instance: %lf\n", (tothop*1.0)/totpath );
  */
  if( ofd2 )fclose(ofd2);
  if( fd )fclose(fd);
}

int jellyfish_calc_diameter(){
  int i, j, shPathLen, diameter;
  Path kpath[JFISH_MAX_NUM_PATH];

  diameter = 0;

  for(i=totPE;i<totNode;i++){
    for(j=totPE;j<totNode;j++){
      if(i !=j ){
        //pass first PE of the selected switch/router 
        //because the routing function expects PEs as source and destination
        if(j%100==0) printf("(%d,%d)",i-totPE, j-totPE);
	jellyfish_kpath_routing(graph[i][get_r(i)+0], graph[j][get_r(j)+0], kpath, 1);
	shPathLen = jellyfish_compute_hop_count(kpath[0]);
	if(shPathLen > diameter) 
	  {
	    diameter = shPathLen;
	    jellyfish_print_path(i, j, kpath[0]);
	  }
      }
    }
  }
  return diameter-3;
}

void jellyfish_calc_avgSP(){
  int i, j, m, n;
  int totCount = 0;
  Path kpath[JFISH_MAX_NUM_PATH];
  double avg[JFISH_MAX_NUM_PATH];

  for(i=0;i<JFISH_MAX_NUM_PATH;i++){
    avg[i] = 0;
  }

  for(i=totPE;i<totNode;i++){
    for(j=totPE;j<totNode;j++){
      //if((i*totSE + j) % 100000 == 0)
	//printf("%d nodes processed ...\n",i*totSE+j);
      if(i !=j ){
	jellyfish_kpath_routing(get_firstPE(i), get_firstPE(j), kpath, 60);
	for(m=0;m<jellyfish_K;m++){
	  int hopCount = jellyfish_compute_hop_count(kpath[m]);
	  for(n=m;n<jellyfish_K;n++){
	    avg[n] += hopCount - 3;
	  }
	}
	totCount++;
      }
    }
  }
  for(i=0;i<jellyfish_K;i++)
    //printf("%d %3.2f\n", i+1, (avg[i]/totCount));
    printf("%d %3.2f\n", i+1, (avg[i]/(totCount*(i+1))));
}
