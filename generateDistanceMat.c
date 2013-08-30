/*THIS FILE COMPILES AS A MATLAB MEX DYNAMIC LIBRARY*/
/*Function to generate intersomatic distance matrix*/
/*Input : 3xN array of 3d points*/
/*Output: NxN array of Distances, MinVal, MaxVal*/

/*stamatiad.st@gmail.com*/

#include <stdio.h>
#include <stdlib.h>

#include "mex.h"
#include "matrix.h"
#include "math.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    /*Declerations:*/
    double *inputPts, *outDistMat ;
    int i,j, NoPts;
    
    
    /*Input: the 3d points' coordinates*/
    inputPts = mxGetPr(prhs[0]); /*concatenated column-wise (Matlab does this..)*/
    NoPts = mxGetN(prhs[0]);
    
    plhs[0] = mxCreateDoubleMatrix(NoPts,NoPts,mxREAL);
    outDistMat = mxGetPr(plhs[0]);

    /*find for all pairs the corresponding distance:*/
    for(i=0;i<NoPts;i++)
        for(j=0;j<NoPts;j++)
            if(i!=j){
        outDistMat[i*NoPts+j] = sqrt((inputPts[i*3]-inputPts[j*3])*(inputPts[i*3]-inputPts[j*3])
        + (inputPts[i*3+1]-inputPts[j*3+1])*(inputPts[i*3+1]-inputPts[j*3+1])
        + (inputPts[i*3+2]-inputPts[j*3+2])*(inputPts[i*3+2]-inputPts[j*3+2]));
            }
    
    return;
}
