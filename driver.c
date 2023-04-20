#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include "topology_var.h"

struct timeval t;
int patharray[NPJFISH_MAX_NUM_PATH][MAX_PATH_LEN];
int allpath[MAX_ALLPATH][MAX_PATH_LEN];

void K_path_with_minimum_linkloads(int K,int p,int r);
void gen_topofile(int r);
void gen_allpath_routefile(char* outputfile, int K);

int main()
{
	//total 7 argument
    // argument will be ./program_name, total_num_switches, number_of_inter_switch_ports, concentration, bwPS, bwSS, K
    int i;
    int nTOR; //number of TOR switches
    int nPE;  //total number of servers
    int totPorts;    // total # of links per TOR
    int p;    // # of server per TOR
    int x;    // total # of extra additional servers distributed among ntor switches
    int K;
    long long int bwPS, bwSS;
    int t1,t2;
    int rc;
    unsigned int seed= 1201;
    int dump_to_file = 1, directional=0;
	int r; //number of inter switch ports
	printf("Enter values : nTOR , radix , concentration, bwPS, bwSS, K\n");
    scanf("%d %d %d %d %d %d",&nTOR,&r,&p,&t1,&t2,&K);
    //nTOR = atoi( argv[1] ); //num of switches
	//r = atoi( argv[2] );
    totPorts    = r+p; //radix of a switch
    //p = atoi(argv[3]); //concentration
    nPE = p * nTOR; //total numbber of compute nodes
    bwPS = t1*((long long int)1000000);
    bwSS = t2*((long long int)1000000);
    //K =  atoi( argv[6]);
    rc = JELLYFISH_KPATH_LLSKR_MIXED_ROUTING;


    jellyfish_topology_init( nTOR, totPorts, nPE, bwPS, bwSS, rc, seed, dump_to_file, directional);
    jellyfish_dump_routing_binary();
    jellyfish_dump_routing_ascii();
    printf("Jellyfish topo building complete\n");
    char a2a_allroutes_filename[100]="";
    sprintf(a2a_allroutes_filename, "rrg_%d_%d_K%d.allroutes",nTOR,r,K);
gen_allpath_routefile(a2a_allroutes_filename, K);
    printf("Calculating K minimum linkload paths\n");
	//K_path_with_minimum_linkloads(K,p,r);
    gen_topofile(r);
    printf("Allroute file generates successfully\n");
	return 0;

  }

void gen_allpath_routefile(char* outputfile, int K)
{
  int i,j,k,jj;
  int path_spread, path_spread_vlb;
  int random_vlb_index;
  int path_len;
  int allhop[MAX_PATH_LEN];
  int onepath[MAX_PATH_LEN];

  FILE *ofd;

  if ((ofd = fopen(outputfile, "w")) == NULL) {
    printf("gen_allpath_routefile: cant open %s  for writing.\n", outputfile);
    exit(0);
  }
  for(i=0;i<totSE;i++){
    for(j=0;j<totSE;j++){
      if(i==j) continue;
        //path_spread = jellyfish_K;
        //jellyfish_readnewkpath(get_firstPE(i+totPE), get_firstPE(j+totPE), patharray);
        path_spread = jellyfish_allpath_routing(get_firstPE(i+totPE), get_firstPE(j+totPE), patharray, allhop, 3);
        int path_spread_4H = jellyfish_allpath_routing(get_firstPE(i+totPE), get_firstPE(j+totPE), allpath, allhop, 4);
        int random_4H_index;
        while(path_spread <K)
        {
          do{
		  random_4H_index = rand() % path_spread_4H;
            for(path_len=0;allpath[random_4H_index][path_len] != -1;path_len++);
          }while(path_len<7); //4 sw-sw hops

          for(int x=0;x<=path_len;x++)  //copy including the -1
            patharray[path_spread][x] = allpath[random_4H_index][x];
          path_spread++;
        }
        path_spread = K;
	  

      printf("%d paths between switches %d and %d\n", path_spread,i,j);
      fprintf(ofd, "%d %d %d\n", i, j,path_spread);

      for(k=0;k<path_spread;k++){
      fprintf(ofd, "\t");
        for(jj=1;patharray[k][jj+1]!=-1;jj++){
          fprintf(ofd, " %d", patharray[k][jj]-totPE);
        }
      fprintf(ofd, "\n");

      }
    }
 }
  if(ofd) fclose(ofd);

}

