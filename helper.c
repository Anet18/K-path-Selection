#include "topology.h"

#define MAX_QUEUE_SIZE 1000000

/*-----------Queue Implementation------------------*/

void QInit(qType *myq)
{
  myq->head = myq->tail = -1;
}

int QIsEmpty(qType *myq)
{
  return (myq->head == myq->tail);
}

int QIsFull(qType *myq)
{
  return (myq->tail - MAX_QUEUE_SIZE == myq->head);
}

void QPush(qType *myq, stackElementType element)
{
  if (QIsFull(myq)) {
    fprintf(stderr, "Can't push element on queue: queue is full.%d,%d\n",myq->tail, myq->head);
    exit(1);  /* Exit, returning error code. */
  }
  
  myq->tail++;
  myq->contents[myq->tail % MAX_QUEUE_SIZE] = element;
  /* Put information in array; update top. */

  
}

stackElementType QPop(qType *myq)
{
  if (QIsEmpty(myq)) {
    fprintf(stderr, "Can't pop element from queue: queue is empty.\n");
    exit(1);  /* Exit, returning error code. */
  }

  myq->head++;
  return myq->contents[myq->head % MAX_QUEUE_SIZE];
}

/*------------------------------------------------------*/
/*-------Bisection Bandwidth---------------------------*/

int search(int i, int *arr, int size)
{
  int l = 0; 
  int r = size;
  while (l < r) {
    if (arr[(l + r) / 2] == i) return 1;
    else if (arr[(l+r)/2] > i) r = (l+r)/2;
    else l = (l+r)/2 + 1;
  }
  return 0;
}

void gen_half(int *arr, int N, int p)
{
  int flag[MAX_NODE], flag1[MAX_NODE];
  int set[MAX_NODE];
  int i, j, k;

  for (i=0; i<N; i++) flag[i] = i;
  for (i=0; i<N/2; i++) {
    k = random() % (N-i);
    set[i] = flag[k];
    for (j=k; j<N-i-1; j++) flag[j] = flag[j+1];
  }
  for (i=0; i<N; i++) flag1[i] = 0;
  for (i=0; i<N/2; i++) flag1[set[i]] = 1;
  i = 0;
  for (j=0; j<N; j++) {
    if (flag1[j] == 1) {arr[i] = j+N*p; i++;}
  }
}

