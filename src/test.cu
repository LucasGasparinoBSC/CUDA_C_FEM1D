#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "defConstants.cuh"
#include "genData.cuh"
#include "convecKernels.cuh"

//__constant__ float xgp_cte[MAX_NGAUS];
__constant__ float wgp_cte[MAX_NGAUS];
//__constant__ float N_cte[MAX_NGAUS*MAX_NNODE];
//__constant__ float dN_cte[MAX_NGAUS*MAX_NNODE];

__global__ void convec_gpuConst(int nelem, int nnode, int ngaus, int npoints, int *connec,
                                float *N, float *dN, float *w, float *u, float *R)
{
    // Create shared memory
    __shared__ float u_shared[8]; // Max 3 nodes per element
    __shared__ float v_shared[8]; // Max 3 Gauss points per element

    // Set ielem and inode according to block and thread indices
    int ielem = blockIdx.x;
    int inode = threadIdx.x;
    int igaus = threadIdx.y;

    // Ensure R is zero
    R[connec[ielem*nnode + inode]] = 0.0f;

    // Fill shared memory
    v_shared[igaus] = 0.0f;
    u_shared[inode] = u[connec[ielem*nnode + inode]];
    __syncthreads();

    // Compute dN*u_shared at each Gauss point
    //for (int jnode = 0; jnode < nnode; jnode++)
    //{
    //	v_shared[igaus] += dN[igaus*nnode + jnode]*u_shared[jnode];
    //}
    //__syncthreads();
    atomicAdd(&v_shared[igaus], dN[igaus*nnode + inode]*u_shared[inode]);

    // Atomically update R
    atomicAdd(&R[connec[ielem*nnode + inode]], wgp_cte[igaus]*N[igaus*nnode + inode]*v_shared[igaus]);
    __syncthreads();
}

int main(void)
{
    // Var for error checking
    int ierr;

    // Set mesh details
    int nelem = 2;
    int nnode = 3;
    int ngaus = 3;
    int npoints = 5;

    // Print mesh info
    printf("*----------*\n");
    printf("nelem = %d\n", nelem);
    printf("nnode = %d\n", nnode);
    printf("ngaus = %d\n", ngaus);
    printf("npoints = %d\n", npoints);
    printf("*----------*\n");

    // Create connectivity table
    int *connec = (int *)malloc(nelem*nnode*sizeof(int));
    connec[0] = 0;
    connec[1] = 1;
    connec[2] = 3;
    connec[3] = 1;
    connec[4] = 2;
    connec[5] = 4;

    // Print connec in a table format
    printf("connec = \n");
    for (int ielem = 0; ielem < nelem; ielem++)
    {
        for (int inode = 0; inode < nnode; inode++)
        {
            printf("%d ", connec[ielem*nnode + inode]);
        }
        printf("\n");
    }
    printf("*----------*\n");

    // Set quadrature points
    float *xgp = (float *)malloc(ngaus*sizeof(float));
    float *wgp = (float *)malloc(ngaus*sizeof(float));
    ierr = quadratureData(ngaus,xgp,wgp);
    if (ierr != 0)
    {
        printf("Error in quadratureData\n");
        return EXIT_FAILURE;
    }
    //wgp[0] = 1.0f;
    //wgp[1] = 1.0f;
    //wgp[2] = 1.0f;

    // Set N and dN
    float *N = (float *)malloc(nnode*ngaus*sizeof(float));
    float *dN = (float *)malloc(nnode*ngaus*sizeof(float));
    for (int igaus = 0; igaus < ngaus; igaus++)
    {
        for (int inode = 0; inode < nnode; inode++)
        {
            N[igaus*nnode + inode] = 1.0f;
            dN[igaus*nnode + inode] = 0.5f;
        }
    }

    // Set initial condition u
    float *u = (float *)malloc(npoints*sizeof(float));
    u[0] = 1.0f;
    u[1] = 2.0f;
    u[2] = 1.0f;
    u[3] = 1.5f;
    u[4] = 1.5f;

    // Print u
    printf("u = \n");
    for (int ipoint = 0; ipoint < npoints; ipoint++)
    {
        printf("%d %f\n", ipoint, u[ipoint]);
    }
    printf("*----------*\n");

    // Call the CPU version of convec
    float *R_cpu = (float *)malloc(npoints*sizeof(float));
    convec_cpu(nelem,nnode,ngaus,npoints,connec,
                N,dN,wgp,u,R_cpu);

    // Print R_cpu
    printf("R_cpu = \n");
    for (int ipoint = 0; ipoint < npoints; ipoint++)
    {
        printf("%d %f\n", ipoint, R_cpu[ipoint]);
    }
    printf("*----------*\n");

    // Create GPU arrays for connec, N, dN, wgp, u, and R_gpu
    int *connec_gpu;
    float *N_gpu, *dN_gpu, *wgp_gpu, *u_gpu, *R_gpu;

    cudaMalloc((void **)&connec_gpu, nelem*nnode*sizeof(int));
    cudaMalloc((void **)&N_gpu, nnode*ngaus*sizeof(float));
    cudaMalloc((void **)&dN_gpu, nnode*ngaus*sizeof(float));
    cudaMalloc((void **)&wgp_gpu, ngaus*sizeof(float));
    cudaMalloc((void **)&u_gpu, npoints*sizeof(float));
    cudaMalloc((void **)&R_gpu, npoints*sizeof(float));

    cudaMemcpy(connec_gpu, connec, nelem*nnode*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(N_gpu, N, nnode*ngaus*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dN_gpu, dN, nnode*ngaus*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(wgp_gpu, wgp, ngaus*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(u_gpu, u, npoints*sizeof(float), cudaMemcpyHostToDevice);

    // Call the basic GPU version of convec
    convec_gpuBasic<<<nelem,nnode>>>(nelem,nnode,ngaus,npoints,connec_gpu,
                                     N_gpu,dN_gpu,wgp_gpu,u_gpu,R_gpu);
    
    // Copy data from GPU to CPU
    float *R_gpuBasic = (float *)malloc(npoints*sizeof(float));
    cudaMemcpy(R_gpuBasic, R_gpu, npoints*sizeof(float), cudaMemcpyDeviceToHost);

    // Print R_gpuBasic
    printf("R_gpuBasic = \n");
    for (int ipoint = 0; ipoint < npoints; ipoint++)
    {
        printf("%d %f\n", ipoint, R_gpuBasic[ipoint]);
    }
    printf("*----------*\n");

    // Call the 1st shared memory GPU version of convec
    convec_gpuShared1<<<nelem,nnode>>>(nelem,nnode,ngaus,npoints,connec_gpu,
                                       N_gpu,dN_gpu,wgp_gpu,u_gpu,R_gpu);

    // Copy data from GPU to CPU
    float *R_gpuShared1 = (float *)malloc(npoints*sizeof(float));
    cudaMemcpy(R_gpuShared1, R_gpu, npoints*sizeof(float), cudaMemcpyDeviceToHost);

    // Print R_gpuShared1
    printf("R_gpuShared1 = \n");
    for (int ipoint = 0; ipoint < npoints; ipoint++)
    {
        printf("%d %f\n", ipoint, R_gpuShared1[ipoint]);
    }
    printf("*----------*\n");

    // Call the 2nd shared memory GPU version of convec
    dim3 block(nnode,ngaus,1);
    dim3 grid(nelem,1,1);
    convec_gpuShared2<<<grid,block>>>(nelem,nnode,ngaus,npoints,connec_gpu,
                                      N_gpu,dN_gpu,wgp_gpu,u_gpu,R_gpu);

    // Copy data from GPU to CPU
    float *R_gpuShared2 = (float *)malloc(npoints*sizeof(float));
    cudaMemcpy(R_gpuShared2, R_gpu, npoints*sizeof(float), cudaMemcpyDeviceToHost);

    // Print R_gpuShared2
    printf("R_gpuShared2 = \n");
    for (int ipoint = 0; ipoint < npoints; ipoint++)
    {
        printf("%d %f\n", ipoint, R_gpuShared2[ipoint]);
    }
    printf("*----------*\n");

    // Fill the constant memory wgp_cte
    cudaMemcpyToSymbol(wgp_cte, wgp, ngaus*sizeof(float), 0, cudaMemcpyHostToDevice);
    
    // Call the constant memory GPU version of convec
    convec_gpuConst<<<grid,block>>>(nelem,nnode,ngaus,npoints,connec_gpu,
                                     N_gpu,dN_gpu,wgp_gpu,u_gpu,R_gpu);

    // Copy data from GPU to CPU
    float *R_gpuConst = (float *)malloc(npoints*sizeof(float));
    cudaMemcpy(R_gpuConst, R_gpu, npoints*sizeof(float), cudaMemcpyDeviceToHost);

    // Print R_gpuConst
    printf("R_gpuConst = \n");
    for (int ipoint = 0; ipoint < npoints; ipoint++)
    {
        printf("%d %f\n", ipoint, R_gpuConst[ipoint]);
    }
    printf("*----------*\n");

    return 0;
}