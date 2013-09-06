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

double commonIncomming(int I,int J,double* connMat,int Numel, double* CBins, double* CProb, int connNumel);
double commonOutgoing(int I,int J,double* connMat,int Numel, double* CBins, double* CProb, int connNumel);
int decideConnection(double dist, double CProb, double *RBins, double *RProb, int RNumel);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    /*Declerations:*/
    double *distMat, *connMat,*tempMat, *outConnMat,*connBins,*incomingProb, *outgoingProb,*recipBins,*recipProb;
    int numel , start, connNumel, recipNumel;
    int i,j, connType;
    double maxProb, *probsMat, TempVal;
    
    distMat = mxGetPr(prhs[0]);
    connMat = mxGetPr(prhs[1]);
    numel = (int)mxGetN(prhs[1]);
    
    /*Connection probabilities*/
    connBins = mxGetPr(prhs[2]);
    incomingProb = mxGetPr(prhs[3]);
    outgoingProb = mxGetPr(prhs[4]);
    connNumel = (int)mxGetN(prhs[4]);
    /*reciprocal probabilities*/
    recipBins = mxGetPr(prhs[5]);
    recipProb = mxGetPr(prhs[6]);
    recipNumel = (int)mxGetN(prhs[6]);
    
    plhs[0] = mxCreateDoubleMatrix(numel,numel,mxREAL);
    outConnMat = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(numel,numel,mxREAL);
    probsMat = mxGetPr(plhs[1]);
    
    tempMat = (double*)mxCalloc(numel*numel,sizeof(double));
    
    srand(time(NULL));
    /*For every unique pair*/
    for(i=0;i<numel;i++){
        for(j=0;j<numel;j++){
            /*individual probabilities(?)*/
			/*Probabilities map is based in pair incomming conns only (from common neighbors) multiplied by*/
			/*pair outgoing conns only (from common neighbors)*/
            tempMat[i*numel+j] = (commonIncomming(i,j,connMat,numel,connBins,incomingProb,connNumel) *
                    commonOutgoing(i,j,connMat,numel,connBins,outgoingProb,connNumel) ) ; /*was /5 (?)*/
        }

    }
    /*Normalize probability:*/
    for(i=0;i<numel;i++){
        for(j=0;j<numel;j++){
            tempMat[i*numel+j] = tempMat[i*numel+j] / (numel-6) ;
            probsMat[i*numel+j] = tempMat[i*numel+j];
        }
    }
    
    start=1;
    for(i=0;i<numel;i++){
        for(j=start;j<numel;j++){
		/*Using Probabilities map similar to InitializeNetwork: Probabilities map is the connection probabilities and*/
		/*the reciprocal probabilities are imported from Perin et al. Fig1 */
            connType = decideConnection(distMat[i*numel+j],tempMat[i*numel+j],recipBins,recipProb,recipNumel);
            
            /* connection types: 0=No connection, 1,2=Single (side determined as in Perin et al. 2013), 3=Reciprocal */
            if(connType == 0){
                outConnMat[i*numel+j] = 0;
                outConnMat[j*numel+i] = 0;
            }
            if(connType == 1){
			/*POSSIBLY WRONG if single connection exists, do not reverse it!! */
                if(!((outConnMat[i*numel+j]==0 && outConnMat[j*numel+i] == 1) ||
                        (outConnMat[i*numel+j]==1 && outConnMat[j*numel+i] == 0))){
                    outConnMat[i*numel+j] = 1;
                    outConnMat[j*numel+i] = 0;}
            }
            if(connType == 2){
			/*POSSIBLY WRONG if single connection exists, do not reverse it!! */
                if(!((outConnMat[i*numel+j]==0 && outConnMat[j*numel+i] == 1) ||
                        (outConnMat[i*numel+j]==1 && outConnMat[j*numel+i] == 0))){
                    outConnMat[i*numel+j] = 0;
                    outConnMat[j*numel+i] = 1;}
            }
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
double commonIncomming(int I,int J,double* connMat,int Numel, double* CBins, double* CProb, int connNumel){
    int i,commonNo, iC;
    commonNo=0;
    for(i=0;i<Numel;i++){
        if(connMat[i*Numel+I] == 1 && connMat[i*Numel+J] == 1)
            commonNo++;
    }
    
    iC=0;
    while((commonNo>CBins[iC]) && (iC<connNumel) ){iC++;}
    return CProb[iC];
}

/*Calculate connection Probability based on No of outgoing common neighborgs*/
double commonOutgoing(int I,int J,double* connMat,int Numel, double* CBins, double* CProb, int connNumel){
    int i,commonNo, iC;
    commonNo=0;
    for(i=0;i<Numel;i++){
        if(connMat[I*Numel+i] == 1 && connMat[J*Numel+i] == 1)
            commonNo++;
    }
    
    iC=0;
    while((commonNo>CBins[iC]) && (iC<connNumel) ){iC++;}
    return CProb[iC];
}

int decideConnection(double dist, double CProb, double *RBins, double *RProb, int RNumel){
    int iR;
    double RND;
    RND = (double)rand()/RAND_MAX;
    
    iR=0;
    while((dist>RBins[iR]) && (iR<RNumel) ){iR++;}
    
    /*printf("@ bin %d disconnection prob is %f\n",iR, (CProb + RProb[iR]) );*/
    
    
    if( RND > (CProb*2+RProb[iR]) ){
        return 0; /*No connection*/
    }else if(RND > (CProb+RProb[iR])){
        return 1;/* one way */
    }else if (RND > (RProb[iR])){
        return 2; /* the other way */
    }else {
        return 3;
    }
}
