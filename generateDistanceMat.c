#include <stdio.h>
#include <stdlib.h>

#include "mex.h"
#include "matrix.h"
#include "math.h"

/*Function to generate intersoma distance matrix*/
/*Input : 3xN array of 3d points*/
/*Output: NxN array of Distances, MinVal, MaxVal*/
/*stefanou@imbb.forth.gr*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    /*Declerations:*/
    double *inputPts, *outDistMat ;
    /*unsigned int *Min, *Max ;*/
    int i,j, NoPts;
    
    
    /*Input: the 3d points' coordinates*/
    inputPts = mxGetPr(prhs[0]); /*concatenated column-wise (Matlab does this..)*/
    NoPts = mxGetN(prhs[0]);
    
    plhs[0] = mxCreateDoubleMatrix(NoPts,NoPts,mxREAL);
    outDistMat = mxGetPr(plhs[0]);
    /*plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    Min = mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
    Max = mxGetPr(plhs[2]);*/


    /*point(i,x) = inputPts[i*3+0], point(i,y) = inputPts[i*3+1],point(i,z) = inputPts[i*3+2]*/
    /*find for all pairs the corresponding distance:*/
    for(i=0;i<NoPts;i++)
        for(j=0;j<NoPts;j++)
            if(i!=j){
        outDistMat[i*NoPts+j] = sqrt((inputPts[i*3]-inputPts[j*3])*(inputPts[i*3]-inputPts[j*3])
        + (inputPts[i*3+1]-inputPts[j*3+1])*(inputPts[i*3+1]-inputPts[j*3+1])
        + (inputPts[i*3+2]-inputPts[j*3+2])*(inputPts[i*3+2]-inputPts[j*3+2]));
            }
    
    /*Search for min/max intersomatic distance inside the network:
    *Min = 10000;
    *Max = 0;
    for(i = 0;i<NoPts*NoPts;i++){
        if(outDistMat[i]> *Max)
            *Max = outDistMat[i];
        if(outDistMat[i]< *Min)
            *Min = outDistMat[i];
    }*/
    
    return;
}
