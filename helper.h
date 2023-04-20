#ifndef HELPER_H
#define HELPER_H

#include "topology.h"

#define MAX_QUEUE_SIZE 1000000

typedef struct {
  int node;
  int path[MAX_PATH_LEN];
  int pid;
} stackElementType;

typedef struct {
  stackElementType contents[MAX_QUEUE_SIZE];
  int head;
  int tail;
} qType;

void QInit(qType *myq);
int QIsEmpty(qType *myq);
int QIsFull(qType *myq);
void QPush(qType *myq, stackElementType element);
stackElementType QPop(qType *myq);

#endif
