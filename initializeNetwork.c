#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "mex.h"
#include "matrix.h"
#include "math.h"

# define DEBUG

/*stamatiad.st@gmail.com*/

int decideConnection(double dist ,double* bias, double *CBins, double *CProb, int CNumel, double *RBins, double *RProb, int RNumel);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    /*Declerations:*/
    double *distMat, *connMat, *connBins, *connProb, *recipBins,  *recipProb;
    double *distBias;
    int numel , start, IR;
    int i,j,k, connType, connNumel, recipNumel;
    int incomming, outgoing;
	
    distMat = mxGetPr(prhs[0]);
    numel = (int)mxGetN(prhs[0]);
	/*Connection probabilities*/
    distBias = mxGetPr(prhs[1]);
    connBins = mxGetPr(prhs[2]);
	connProb = mxGetPr(prhs[3]);
    connNumel = (int)mxGetN(prhs[3]);
	/*reciprocal probabilities*/
	recipBins = mxGetPr(prhs[4]);
	recipProb = mxGetPr(prhs[5]);
    recipNumel = (int)mxGetN(prhs[5]);
	
	plhs[0] = mxCreateDoubleMatrix(numel,numel,mxREAL);
    connMat = mxGetPr(plhs[0]);
	
	srand(time(NULL));
	/*For every unique pair*/
	start=1;
	for(i=0;i<numel;i++){
		for(j=start;j<numel;j++){
            /*Decide connection type (no, single, reciprocal) based on distance and imported probabilities */
                connType = decideConnection(distMat[i*numel+j],distBias,connBins,connProb,connNumel,recipBins,recipProb,recipNumel);
            /* connection types: 0=No connection, 1,2=Single (side determined as in Perin et al. 2013), 3=Reciprocal */
			if(connType == 0){
				connMat[i*numel+j] = 0;
				connMat[j*numel+i] = 0;
			}
			if(connType == 1){
					connMat[i*numel+j] = 1;
					connMat[j*numel+i] = 0;}
			if(connType == 2){
					connMat[i*numel+j] = 0;
					connMat[j*numel+i] = 1;}
			
			if(connType == 3){
				connMat[i*numel+j] = 1;
				connMat[j*numel+i] = 1;
			}
		}
		start++;
	}
	
    return;
}

int decideConnection(double dist,double* bias, double *CBins, double *CProb, int CNumel, double *RBins, double *RProb, int RNumel){
	int iC, iR;
    double RND, CSp, CRp;
    RND = (double)rand()/RAND_MAX;
    
	iC=0;
	while((dist>CBins[iC]) && (iC<CNumel)){iC++;}
    
    iR=0;
	while((dist>RBins[iR]) && (iR<RNumel)){iR++;}
    
    CSp = CProb[iC] / bias[iC];
    CRp = RProb[iR] / bias[iR];

	if( RND > (CSp+CRp) ){
        return 0; /*No connection*/
    }else if(RND > ((CSp/2)+CRp)){
         return 1;/* one way */
    }else if (RND > (CRp)){
        return 2; /* the other way */
    }else {
        return 3;
    }
}
