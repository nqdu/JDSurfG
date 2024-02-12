#pragma once 
#include <omp.h>
void allocate_tasks(int ntasks,int nprocs,int myrank, int &start,int &end);