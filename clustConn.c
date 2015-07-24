#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "mex.h"
#include "matrix.h"
#include "math.h"

/*stamatiad.st@gmail.com*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    /*Declerations:*/
    double *allCombs, *idxCombs, *connMat, *out ;
    int compsLen, idxLen, diffComp, clustLen, _i, _j, i,j, counter, matSide;
	
    allCombs = mxGetPr(prhs[0]);
    clustLen = (int)mxGetM(prhs[0]);
    compsLen = (int)mxGetN(prhs[0]);
    idxCombs = mxGetPr(prhs[1]);
    idxLen = (int)mxGetN(prhs[1]);
    connMat = mxGetPr(prhs[2]);
    matSide = (int)mxGetN(prhs[2]);
    
	plhs[0] = mxCreateDoubleMatrix(compsLen,1,mxREAL);
    out = mxGetPr(plhs[0]);

    for(i=0;i<compsLen;i++){
        counter = 0;
        for(j=0;j<idxLen;j++){
            _i = (int)allCombs[i*clustLen+(int)idxCombs[j*2]];
            _j = (int)allCombs[i*clustLen+(int)idxCombs[j*2+1]];
            counter += (int)connMat[ _i * matSide + _j];
        }
        out[i] = counter;
    }

    return;
}
