#ifndef CONVEC_KERNELS_H
#define CONVEC_KERNELS_H

#include <cuda.h>
#include <cuda_runtime.h>

// CPU version
void convec_cpu(int nelem, int nnode, int ngaus, int npoints, int *connec,
                float *N, float *dN, float *w, float *u, float *R);

// Basic GPU version
__global__ void convec_gpuBasic(int nelem, int nnode, int ngaus, int npoints, int *connec,
                                float *N, float *dN, float *w, float *u, float *R);

// Shared memory GPU version 1
__global__ void convec_gpuShared1(int nelem, int nnode, int ngaus, int npoints, int *connec,
                                  float *N, float *dN, float *w, float *u, float *R);

// Shared memory GPU version 2
__global__ void convec_gpuShared2(int nelem, int nnode, int ngaus, int npoints, int *connec,
    float *N, float *dN, float *w, float *u, float *R);


#endif // CONVEC_KERNELS_H