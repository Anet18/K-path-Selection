#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "topology.h"
#include "topology_var.h"

// The _STOREDPATH is specified at the makefile
#ifndef _STOREDPATH
#define _STOREDPATH 0
#endif

//this is used to select the minimum number of hops that we consider as minimal paths
#define MIN_HOP 20

#ifdef _LARGE_PATTERN
#define MAX_NUM_PATH 100000000
#else
#define MAX_NUM_PATH 5000000
#endif

#if _STOREDPATH == 1

//int storedpath[MAX_NUM_PATH][MAX_PATH_LEN];

#endif


//double graph_linkload[totSE][totSE];


void K_path_with_minimum_linkloads(int K,int p,int r)
{
	FILE *fd,*ofd;

	double graph_linkload[totSE][totSE];

	//array to hold the selected k_paths
	short k_paths[totSE][totSE][K][MAX_PATH_LEN];

	//this array is used to keep track of already selected paths
	int flag [totSE][totSE];
	for(int i=0;i<totSE;i++)
    {
        for(int j=0;j<totSE;j++)
        {
            flag[i][j]=0;
        }
    }

	short used_path[totSE][totSE][K+1];

	//initialize graph_linkload to 0
	for(int i=0;i<totSE;i++)
	{
		for(int j=0;j<totSE;j++)
		{
			graph_linkload[i][j]=0.0;
		}
	}

	//calculate k paths
	for(int k=0;k<K;k++)
	{
		for(int s=0;s<totSE;s++)
		{
			for(int d=0;d<totSE;d++)
			{
				if(s!=d) // remove source to source path
				{
				    //printf("%d %d\n",s,d);
					int allpath[JFISH_MAX_NUM_PATH][MAX_PATH_LEN];
					int allhop[MAX_PATH_LEN];
					int maxhop = 4;
					int x = jellyfish_allpath_routing(s*p,d*p,allpath,allhop,maxhop);
					
					//printf("x = %d\n",x);
                    /*
					for(int i=0;i<x;i++)
                    {
                        for(int j=0;j<MAX_PATH_LEN;j++)
                        {
                            printf("%d ",allpath[i][j]);
                        }
                        printf("\n");
                    }
					*/
					
					int hop_count[x];
					double max_linkload[x];

					//calculate hop counts for all paths
					for(int i=0;i<x;i++)
					{
						int hop = 0;
						for(int j=0;j<MAX_PATH_LEN;j++)
						{
							if(allpath[i][j]== -1)
                                break;
                            else
								hop++;
						}
						hop_count[i]=hop-3;
						//printf("%d ",hop_count[i]);
					}
                    //printf("\n");
					//set already selected paths hop count to a large value so that it appears at the last
					if(flag[s][d] == 1)
					{
						int t=0;
						while(used_path[s][d][t]!= -1)
						{
							hop_count[used_path[s][d][t]] = 100000;
							t++;
						}
					}

					//calculate max_linkload for all paths
					for(int i=0;i<x;i++)
					{
						double mx = -99999.0;
						int j=1;
							while(allpath[i][j+2] != -1)
							{
							    //printf("%lf",graph_linkload[allpath[i][j]-nPE][allpath[i][j+1]-nPE]);
								if(graph_linkload[allpath[i][j]-totPE][allpath[i][j+1]-totPE]>mx)
								{
									mx = graph_linkload[allpath[i][j]-totPE][allpath[i][j+1]-totPE];
								}
								j++;
							}
						max_linkload[i]= mx;
						//printf("%lf ",max_linkload[i]);
					}
					//printf("\n");

					//calculating cost for all paths
					double cost_array[x];
					for(int i=0;i<x;i++)
					{
						cost_array[i]= hop_count[i]*max_linkload[i];
						//printf("%lf ",cost_array[i]);
					}
                    //printf("\n");
					
					//list all the short hop paths that are not already selected
					int shortpathlist[x];
					for(int i=0;i<x;i++)
					{
						shortpathlist[i]=0;
					}
					int short_path_count = 0;
					for(int i=0;i<x;i++)
					{
						if(hop_count[i]<=MIN_HOP){
							shortpathlist[i]=1;
                                                        if(s==0 && d==1)
                                                        printf("index = %d ",i);
							short_path_count=short_path_count+1;
						}
     					}
                                        if(s==0 && d==1){
                                        printf("\n");
                                        printf("S=%d D=%d C = %d\n",s,d,short_path_count);
                                        }
					//we need to select all the short paths first if they are not already selected then we will select from long paths
					int short_paths_index_list[short_path_count];
					if(short_path_count!=0) //there exists at least 1 unselected short path
					{
						int s1=0;
						for(int i=0;i<x;i++)
						{
							if(shortpathlist[i]==1)
							{
								short_paths_index_list[s1]=i;
                                                                if(s==0 && d==1)
                                                                printf("short path index= %d ",short_paths_index_list[s1]);
								s1++;
							}
						}
						if(s==0 && d==1)
                                                printf("\n");
						//find minimum cost shortest path
						double min = 999999.0;
					    int min_path_index;
					    for(int i=0;i<short_path_count;i++)
					    {
						    if(cost_array[short_paths_index_list[i]]<min)
						    {
							    min = cost_array[short_paths_index_list[i]];
							    min_path_index = short_paths_index_list[i];
						    }
					    }
                                            if(s==0 && d==1)
                                            printf("min path = %d\n",min_path_index);
                                             
						
						//update the linkload of min_path
					    int h=1;
					    while(allpath[min_path_index][h+2] != -1)
					    {
						     graph_linkload[allpath[min_path_index][h]-totPE][allpath[min_path_index][h+1]-totPE] +=(1.0/K);
				                  		     h++;
					    }
						
						//select min_path in kpaths
					    int l=0;
					    //copy whole path including src and dst node and -1.
					    while(allpath[min_path_index][l]!= -1)
					    {
						    k_paths[s][d][k][l]= allpath[min_path_index][l];
						    l++;
					    }
					    k_paths[s][d][k][l] = -1;
						
						//delete this min_path from the pathlist
					    used_path[s][d][k]= min_path_index;
					    used_path[s][d][k+1] = -1;
					    flag[s][d]= 1;	
					}
					else //all shortest paths are selected 
					{
					
					//............................................
					//find the minimum cost path
					double min = 999999.0;
					int min_path_index;
					for(int i=0;i<x;i++)
					{
						if(cost_array[i]<min)
						{
							min = cost_array[i];
							min_path_index = i;
						}
					}
                    //printf("min path %d \n",min_path_index);
				    //update the linkload of min_path
					int h=1;
					while(allpath[min_path_index][h+2] != -1)
					{
						graph_linkload[allpath[min_path_index][h]-totPE][allpath[min_path_index][h+1]-totPE] += 1.0/K;
						h++;
					}

					//select min_path in kpaths
					int l=0;
					//copy whole path including src and dst node and -1.
					while(allpath[min_path_index][l]!= -1)
					{
						k_paths[s][d][k][l]= allpath[min_path_index][l];
						l++;
					}
					k_paths[s][d][k][l] = -1;

					//delete this min_path from the pathlist
					used_path[s][d][k]= min_path_index;
					used_path[s][d][k+1] = -1;
					flag[s][d]= 1;
					
					}//end of else
				}
				//printf("%d %d\n",s,d);
			}
		}
	}

	//print the kpaths for each source and destination switch in the output file
	char outfilename[50];
	sprintf(outfilename, "rrg_%d_%d_K%d.allroutes",totSE,r,K);

	if ((ofd = fopen(outfilename, "w")) == NULL)
	{
		printf("canot open file %s for writing.\n", outfilename);
		exit(0);
   }

	for(int i=0;i<totSE;i++)
	{
		for(int j=0;j<totSE;j++)
		{
			if(i!=j)
			{
				fprintf(ofd, "%d %d %d \n",i,j,K);
				for(int k=0;k<K;k++)
				{
					fprintf(ofd,"\t");
					int h=1;
					while(k_paths[i][j][k][h+1]!= -1)
					{
						fprintf(ofd,"%d ",k_paths[i][j][k][h]-totPE);
						h++;
					}
					fprintf(ofd,"\n");
				}
			}
		}
	}
}

void gen_topofile(int r)
{
    int i,j;
	FILE *ofd;
    char outfilename[50];
	sprintf(outfilename, "rrg_%d_%d.topology",totSE,r);

	if ((ofd = fopen(outfilename, "w")) == NULL)
	{
		printf("canot open file %s for writing.\n", outfilename);
		exit(0);
   }

  for(i=totPE;i<totNode;i++){
    for(j=0;j<MAX_DEGREE;j++){
      if(graph[i][j]>=totPE) fprintf(ofd, "%d ", graph[i][j]-totPE);
    }
    fprintf(ofd, "\n");
  }
  
  if(ofd) fclose(ofd);
  
}





