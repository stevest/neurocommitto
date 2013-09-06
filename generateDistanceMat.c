/*THIS FILE COMPILES AS A MATLAB MEX DYNAMIC LIBRARY*/
/*Function to generate intersomatic distance matrix*/
/*Input :1: 3xN array of 3d points*/
/*      :2: Wrap dimensions flag */
/*      :3: Size of cubic lattice (aproximate, because lattice is not cubic)
/*      :   used to wrap dimensionality*/
/*Output: NxN array of Distances*/

/*stamatiad.st@gmail.com*/

#include <stdio.h>
#include <stdlib.h>

#include "mex.h"
#include "matrix.h"
#include "math.h"

typedef struct Point{
    double x;
    double y;
    double z;
}Point;

double euclidean(Point a, Point b);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    /*Declerations:*/
    double *inputPts, *outDistMat, dist, distTemp, side ;
    int i,j,x,y,z, NoPts,wrapDims=0,start;
    
    
    /*Input: the 3d points' coordinates*/
    inputPts = mxGetPr(prhs[0]); /*concatenated column-wise (Matlab does this..)*/
    NoPts = mxGetN(prhs[0]);
    
    /* Argument checks */
    if(nrhs == 1){
        wrapDims = 0;}
    if(nrhs == 2){
        if(*mxGetPr(prhs[1]) > 0){
            mexErrMsgTxt("Three input arguments required: \n    Please input the side of the cubic lattice in um!");}}
    if(nrhs > 2){
        if(*mxGetPr(prhs[1]) > 0){
            wrapDims = 1;
            side = *mxGetPr(prhs[2]);} /*square side in um)*/
    }
    
    plhs[0] = mxCreateDoubleMatrix(NoPts,NoPts,mxREAL);
    outDistMat = mxGetPr(plhs[0]);
    
    /*find for all pairs the corresponding distance:*/
    start=1;
    for(i=0;i<NoPts;i++){
        for(j=start;j<NoPts;j++){
            if(i!=j){
                Point Pi;
                Point Pj;
                Pi.x = inputPts[i*3];
                Pi.y = inputPts[i*3+1];
                Pi.z = inputPts[i*3+2];
                Pj.x = inputPts[j*3];
                Pj.y = inputPts[j*3+1];
                Pj.z = inputPts[j*3+2];
                if(wrapDims){
                    dist = sqrt(3) * side; /*set to maximum value ie diagonal*/
                    /*for al permutations*/
                    for(x=-side;x<=side;x+=side){
                        for(y=-side;y<=side;y+=side){
                            for(z=-side;z<=side;z+=side){
                                Point Pt = Pj;
                                Pt.x += x;
                                Pt.y += y;
                                Pt.z += z;
                                distTemp  =  euclidean(Pi, Pt);
                                if(distTemp<dist)
                                    dist = distTemp;
                            }
                        }
                    }
                }
                else {
                    dist  =  euclidean(Pi, Pj);
                }
                outDistMat[i*NoPts+j] = dist ;
                outDistMat[j*NoPts+i] = dist ;
            }
        }
        start++;
    }
    
    return;
}

double euclidean(Point a, Point b){
    return sqrt((a.x-b.x)*(a.x-b.x)
    + (a.y-b.y)*(a.y-b.y)
    + (a.z-b.z)*(a.z-b.z));
}
