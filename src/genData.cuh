#ifndef GENDATA_H
#define GENDATA_H

#include <cuda.h>
#include <cuda_runtime.h>
#include "defConstants.cuh"

int quadratureData(int ngaus, float *xgp, float *wgp);
int elementData(int ngaus, int nnode, float *xgp, float *N, float *dN);

#endif // GENDATA_H