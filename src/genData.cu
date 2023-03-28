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
        // Set xgp and wgp to 0 usingg memset
        memset(xgp, 0, sizeof(float)*ngaus);
        memset(wgp, 0, sizeof(float)*ngaus);
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
    return 0;
}