#ifndef JDSURFG_SHARED_PARALLEL_H_
#define JDSURFG_SHARED_PARALLEL_H_

#include <omp.h>
void allocate_tasks(int ntasks,int nprocs,int myrank, int &start,int &end);

#endif // end JDSURFG_SHARED_PARALLEL_H_
