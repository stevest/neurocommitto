#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "mex.h"
#include "matrix.h"
#include "math.h"


/*stamatiad.st@gmail.com*/
#ifndef _MSC_VER
#define max(a,b) a>b?a:b
#endif

double commonIncomming(int I,int J,double* connMat,int Numel, double* CBins, double* CProb, int connNumel, double *CommNeighbors);
double commonOutgoing(int I,int J,double* connMat,int Numel, double* CBins, double* CProb, int connNumel, double *CommNeighbors);
int decideConnection(double dist, double CProb, double* bias, double *RBins, double *RProb, int RNumel);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	/*Declerations:*/
	double *distMat, *connMat,*tempMat, *outConnMat,*connBins,*incomingProb, *outgoingProb,*recipBins,*recipProb;
	double* distBias, *CommNeighbors;
	int numel , start, connNumel, recipNumel;
	int i,j, connType, totalInput;
	double L2norm, *probsMat, TempVal;

	distMat = mxGetPr(prhs[0]);
	connMat = mxGetPr(prhs[1]);
	numel = (int)mxGetN(prhs[1]);
	distBias = mxGetPr(prhs[2]);

	/*Connection probabilities*/
	connBins = mxGetPr(prhs[3]);
	incomingProb = mxGetPr(prhs[4]);
	outgoingProb = mxGetPr(prhs[5]);
	connNumel = (int)mxGetN(prhs[5]);
	/*reciprocal probabilities*/
	recipBins = mxGetPr(prhs[6]);
	recipProb = mxGetPr(prhs[7]);
	recipNumel = (int)mxGetN(prhs[7]);
	totalInput = (int)(*mxGetPr(prhs[8]));

	plhs[0] = mxCreateDoubleMatrix(numel,numel,mxREAL);
	outConnMat = mxGetPr(plhs[0]);
	plhs[1] = mxCreateDoubleMatrix(numel,numel,mxREAL);
	probsMat = mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(numel,numel,mxREAL);
	CommNeighbors = mxGetPr(plhs[2]);

	tempMat = (double*)mxCalloc(numel*numel,sizeof(double));

	srand(time(NULL));
	/*For every unique pair*/
	start = 1;
	for(i=0;i<numel;i++){
		for(j=start;j<numel;j++){
			/*individual probabilities(?)*/
			/*Probabilities map is based in pair incomming conns only (from common neighbors) multiplied by*/
			/*pair outgoing conns only (from common neighbors)*/
			if( i!=j){
				tempMat[i*numel+j] = (commonIncomming(i,j,connMat,numel,connBins,incomingProb,connNumel,CommNeighbors) +
					commonOutgoing(i,j,connMat,numel,connBins,outgoingProb,connNumel,CommNeighbors) ) ; /*was /5 (?)*/
				if(tempMat[i*numel+j] < 1e-5){
					tempMat[i*numel+j]  = 0;
				}
			}
		}
		start++;
	}

	start = 0;
	for(i=0;i<numel;i++){
		for(j=0;j<start;j++){
			tempMat[i*numel+j] = tempMat[j*numel+i];
		}
		start++;
	}


	/*Normalize probability:*/
	for(j=0;j<numel;j++){
		L2norm = 0;
		for(i=0;i<numel;i++){
			L2norm += tempMat[i*numel+j] * tempMat[i*numel+j];
		}
		L2norm = sqrt(L2norm) / totalInput ;
		for(i=0;i<numel;i++){
			tempMat[i*numel+j] = tempMat[i*numel+j] / L2norm ;
			probsMat[i*numel+j] = tempMat[i*numel+j];
		}
	}

	start=1;
	for(i=0;i<numel;i++){
		for(j=start;j<numel;j++){
			/*Using Probabilities map similar to InitializeNetwork: Probabilities map is the connection probabilities and*/
			/*the reciprocal probabilities are imported from Perin et al. Fig1 */
			connType = decideConnection(distMat[i*numel+j],tempMat[i*numel+j],distBias,recipBins,recipProb,recipNumel);

			/* connection types: 0=No connection, 1,2=Single (side determined as in Perin et al. 2013), 3=Reciprocal */
			if(connType == 0){
				outConnMat[i*numel+j] = 0;
				outConnMat[j*numel+i] = 0;
			}
			if(connType == 1){
				/*POSSIBLY WRONG if single connection exists, do not reverse it!! */
				/*if(!((outConnMat[i*numel+j]==0 && outConnMat[j*numel+i] == 1) ||
				(outConnMat[i*numel+j]==1 && outConnMat[j*numel+i] == 0))){*/
				outConnMat[i*numel+j] = 1;
				outConnMat[j*numel+i] = 0;}
			/*}*/
			if(connType == 2){
				/*POSSIBLY WRONG if single connection exists, do not reverse it!! */
				/*if(!((outConnMat[i*numel+j]==0 && outConnMat[j*numel+i] == 1) ||
				(outConnMat[i*numel+j]==1 && outConnMat[j*numel+i] == 0))){*/
				outConnMat[i*numel+j] = 0;
				outConnMat[j*numel+i] = 1;}
			/*}*/
			if(connType == 3){
				outConnMat[i*numel+j] = 1;
				outConnMat[j*numel+i] = 1;
			}


		}
		start++;
	}

	mxFree(tempMat);
	return;
}


/*Calculate connection Probability based on No of incomminb common neighborgs*/
double commonIncomming(int I,int J,double* connMat,int Numel, double* CBins, double* CProb, int connNumel, double *CommNeighbors){
	int i,commonNo, iC;
	commonNo=0;
	for(i=0;i<Numel;i++){
		if(connMat[i*Numel+I] == 1 && connMat[i*Numel+J] == 1)
			commonNo++;
	}
    CommNeighbors[I*Numel+J] = commonNo;
	iC=0;
	while((commonNo>CBins[iC]) && (iC<connNumel-1) ){iC++;}
	return CProb[iC];
}

/*Calculate connection Probability based on No of outgoing common neighborgs*/
double commonOutgoing(int I,int J,double* connMat,int Numel, double* CBins, double* CProb, int connNumel, double* CommNeighbors){
	int i,commonNo, iC;
	commonNo=0;
	for(i=0;i<Numel;i++){
		if(connMat[I*Numel+i] == 1 && connMat[J*Numel+i] == 1)
			commonNo++;
	}
    CommNeighbors[I*Numel+J] = commonNo;
	iC=0;
	while((commonNo>CBins[iC]) && (iC<connNumel-1) ){iC++;}
	return CProb[iC];
}

int decideConnection(double dist, double CProb, double* bias, double *RBins, double *RProb, int RNumel){
	int iR;
	double RND, CSp, CRp;
	RND = (double)rand()/RAND_MAX;

	iR=0;
	while((dist>RBins[iR]) && (iR<RNumel-1) ){iR++;}

	CRp = RProb[iR] ;/*/ bias[iR];*/
    CSp = CProb + CRp ;/*/ bias[iR];*/
	/*printf("@ bin %d disconnection prob is %f\n",iR, (CProb + RProb[iR]) );*/

	if( RND > (CSp+CRp) ){
		return 0; /*No connection*/
	}else if(RND > (CSp/2+CRp)){
		return 1;/* one way */
	}else if (RND > (CRp)){
		return 2; /* the other way */
	}else {
		return 3;
	}
}
