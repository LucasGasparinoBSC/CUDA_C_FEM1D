#include "genData.cuh"
#include <stdio.h>

int quadratureData(int ngaus, float *xgp, float *wgp)
{
    int ierr = 0;
    // Only compute if ngaus < MAX_NGAUS
    if (ngaus > MAX_NGAUS) {
        ierr = 10;
        printf("Error %d: ngaus > MAX_NGAUS\n", ierr);
    } else {
        // Set xgp and wgp to the correct values (for now just some bs values)
        for (int igaus = 0; igaus < ngaus; igaus++)
        {
            xgp[igaus] = 0.5f;
            wgp[igaus] = 1.0f;
        }
    }
    // Return error code
    return ierr;
}

int elementData(int ngaus, int nnode, float *xgp, float *N, float *dN)
{
    int ierr = 0;
    // Check that ngaus and nnode are valid
    if (ngaus > MAX_NGAUS) {
        ierr = 10;
        printf("Error %d: ngaus > MAX_NGAUS\n", ierr);
        return 10;
    }

    if (nnode > MAX_NNODE) {
        ierr = 11;
        printf("Error %d: nnode > MAX_NNODE\n", ierr);
        return 11;
    }

    // Compute N and dN (for now just some dummy test vars)
    for (int igaus = 0; igaus < ngaus; igaus++)
    {
        for (int inode = 0; inode < nnode; inode++)
        {
            N[igaus*nnode + inode] = 1.0f;
            dN[igaus*nnode + inode] = 0.5f;
        }
    }
    return ierr;
}