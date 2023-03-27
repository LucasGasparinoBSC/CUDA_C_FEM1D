#include "convecKernels.cuh"
#include <string.h>
#include <stdio.h>

void convec_cpu(int nelem, int nnode, int ngaus, int npoints, int *connec,
				float *N, float *dN, float *w, float *u, float *R)
{
	// Create variables
	float v;
	// Ensure R is zero using memset
	memset(R, 0, npoints*sizeof(float));
	// Loop over elements
	for (int ielem = 0; ielem < nelem; ielem++)
	{
		printf("ielem  |  igaus  |  jnode  |  v\n");
		// Loop over quadrature points
		for (int igaus = 0; igaus < ngaus; igaus++)
		{
			// Compute v = dot(dN,u)
			v = 0.0f;
			for (int jnode = 0; jnode < nnode; jnode++)
			{
				v += dN[igaus*nnode + jnode]*u[connec[ielem*nnode + jnode]];
				printf("%d  |  %d  |  %d  |  %f\n", ielem, igaus, jnode, v);
			}
			// Print v
			// Compute quadrature and store on global assembly array
			for (int inode = 0; inode < nnode; inode++)
			{
				R[connec[ielem*nnode + inode]] += w[igaus]*N[igaus*nnode + inode]*v;
			}
		}
	}
	printf("*----------*\n");
}

__global__ void convec_gpuBasic(int nelem, int nnode, int ngaus, int npoints, int *connec,
								float *N, float *dN, float *w, float *u, float *R)
{
	// Creatte variables
	float v;

	// Set ielem and inode according to block and thread indices
	int ielem = blockIdx.x;
	int inode = threadIdx.x;

	// Ensure R is zero
	R[connec[ielem*nnode + inode]] = 0.0f;

	// Loop over Gauss points
	for (int igaus = 0; igaus < ngaus; igaus++)
	{
		// Compute v = dot(dN,u)
		v = 0.0f;
		for (int jnode = 0; jnode < nnode; jnode++)
		{
			v += dN[igaus*nnode + jnode]*u[connec[ielem*nnode + jnode]];
		}
		// Atomically update R
		atomicAdd(&R[connec[ielem*nnode + inode]], w[igaus]*N[igaus*nnode + inode]*v);
	}
}

__global__ void convec_gpuShared1(int nelem, int nnode, int ngaus, int npoints, int *connec,
									float *N, float *dN, float *w, float *u, float *R)
{
	// Create variables
	float v;

	// Create shared memory for u
	// Static allocation of 8 to account for max nodes per element
	__shared__ float u_shared[8];

	// Set ielem and inode according to block and thread indices
	int ielem = blockIdx.x;
	int inode = threadIdx.x;

	// Ensure R is zero
	R[connec[ielem*nnode + inode]] = 0.0f;

	// Fill shared memory
	if (inode > nnode) {
		u_shared[inode] = 0.0f;
	}else{
		u_shared[inode] = u[connec[ielem*nnode + inode]];
	}
	__syncthreads();

	// Loop over Gauss points
	for (int igaus = 0; igaus < ngaus; igaus++)
	{
		// Compute v = dot(dN,u)
		v = 0.0f;
		for (int jnode = 0; jnode < nnode; jnode++)
		{
			v += dN[igaus*nnode + jnode]*u_shared[jnode];
		}
		// Atomically update R
		atomicAdd(&R[connec[ielem*nnode + inode]], w[igaus]*N[igaus*nnode + inode]*v);
	}
}

__global__ void convec_gpuShared2(int nelem, int nnode, int ngaus, int npoints, int *connec,
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
	atomicAdd(&v_shared[igaus], dN[igaus*nnode + inode]*u_shared[inode]);
	//__syncthreads();

	// Atomically update R
	atomicAdd(&R[connec[ielem*nnode + inode]], w[igaus]*N[igaus*nnode + inode]*v_shared[igaus]);
}