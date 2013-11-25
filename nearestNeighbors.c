#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "mex.h"
#include "matrix.h"
#include "math.h"


/*stamatiad.st@gmail.com*/

void commonIncomming(int I,int J,double* connMat,int Numel, double *CommNeighbors);
void commonOutgoing(int I,int J,double* connMat,int Numel, double *CommNeighbors);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	/*Declerations:*/
	double *connMat;
	double* CommNeighbors;
	int i,j,numel;

	connMat = mxGetPr(prhs[0]);
	numel = (int)mxGetN(prhs[0]);

    plhs[0] = mxCreateDoubleMatrix(numel,numel,mxREAL);
	CommNeighbors = mxGetPr(plhs[0]);

	srand(time(NULL));
	for(i=0;i<numel;i++){
		for(j=0;j<numel;j++){
			if( i!=j){
				commonIncomming(i,j,connMat,numel,CommNeighbors);
                commonOutgoing(i,j,connMat,numel,CommNeighbors) ;
			}
		}
	}
	return;
}

/*Calculate connection Probability based on No of incomminb common neighborgs*/
void commonIncomming(int I,int J,double* connMat,int Numel, double *CommNeighbors){
	int i,commonNo;
	commonNo=0;
	for(i=0;i<Numel;i++){
		if(connMat[i*Numel+I] == 1 && connMat[i*Numel+J] == 1)
			commonNo++;
	}
    CommNeighbors[I*Numel+J] = commonNo;
	return ;
}

/*Calculate connection Probability based on No of outgoing common neighborgs*/
void commonOutgoing(int I,int J,double* connMat,int Numel, double* CommNeighbors){
	int i,commonNo;
	commonNo=0;
	for(i=0;i<Numel;i++){
		if(connMat[I*Numel+i] == 1 && connMat[J*Numel+i] == 1)
			commonNo++;
	}
    CommNeighbors[I*Numel+J] = commonNo;
	return ;
}


