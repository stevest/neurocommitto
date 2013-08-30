#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "mex.h"
#include "matrix.h"
#include "math.h"


/*stefanou@imbb.forth.gr*/

double connProbability(int I,int J,double* connMat,int Numel, double* CBins, double* CProb, int connNumel);
double commonIncomming(int I,int J,double* connMat,int Numel, double* CBins, double* CProb, int connNumel);
double commonOutgoing(int I,int J,double* connMat,int Numel, double* CBins, double* CProb, int connNumel);
int decideConnection(double dist, double CProb, double *RBins, double *RProb, int RNumel);
/*double max(double a, double b);*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    /*Declerations:*/
    double *distMat, *connMat, *outConnMat,*connBins,*incomingProb, *outgoingProb,*recipBins,*recipProb;
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
    
    double* tempMat = (double*)mxCalloc(numel*numel,sizeof(double));
    
    srand(time(NULL));
    /*For every unique pair*/
    for(i=0;i<numel;i++){
        for(j=0;j<numel;j++){
            /*individual probabilities(?)*/
            tempMat[i*numel+j] = (commonIncomming(i,j,connMat,numel,connBins,incomingProb,connNumel) +
                    commonOutgoing(i,j,connMat,numel,connBins,outgoingProb,connNumel) ) ; /*was /5 (?)*/
        }
        /*printf("@ i %d j %d %f\n",i,j,tempMat[i*numel+j] );*/
        /*printf("@ i %d j %d %f\n",i,j,commonIncomming(i,j,connMat,numel,connBins,incomingProb,connNumel) );*/
    }
    /*Normalize probability:(NOT IN THE PAPER)*/
    maxProb = 0;
    for(i=0;i<numel;i++){
        for(j=0;j<numel;j++){
           if( tempMat[i*numel+j] > maxProb)
               maxProb = tempMat[i*numel+j];
        }
    }
    printf("Max probability is: %f\n",maxProb);
    TempVal = ((20*20)/(double)numel) / maxProb;
    printf("Temp is: %f\n",TempVal);
    for(i=0;i<numel;i++){
        for(j=0;j<numel;j++){
            tempMat[i*numel+j] = tempMat[i*numel+j] * TempVal ;
            probsMat[i*numel+j] = tempMat[i*numel+j];
        }
    }
    
    /*
     * start=1;
     * for(i=0;i<numel;i++){
     * for(j=start;j<numel;j++){
     * tempMat[i*numel+j] *= 30/(N-1);
     * }
     * start++;
     * }*/
    
    
    start=1;
    for(i=0;i<numel;i++){
        for(j=start;j<numel;j++){
            /*printf("@ i %d j %d %f\n",i,j,tempMat[i*numel+j] );*/
            connType = decideConnection(distMat[i*numel+j],tempMat[i*numel+j],recipBins,recipProb,recipNumel);
            
            /*not using start*/
            /*
            if(connType == 0){
                outConnMat[i*numel+j] = 0;
            }
            if(connType == 1){
                outConnMat[i*numel+j] = 1;
            }
            if(connType == 2){
                outConnMat[i*numel+j] = 1;
                outConnMat[j*numel+i] = 1;
            }
            */
            
            /*using start*/
            if(connType == 0){
                outConnMat[i*numel+j] = 0;
                outConnMat[j*numel+i] = 0;
            }
            if(connType == 1){
                if(!((outConnMat[i*numel+j]==0 && outConnMat[j*numel+i] == 1) ||
                        (outConnMat[i*numel+j]==1 && outConnMat[j*numel+i] == 0)))
                {
                if(((double)rand()/RAND_MAX) <= 0.5){
                    outConnMat[i*numel+j] = 1;
                    outConnMat[j*numel+i] = 0;}
                else{
                    outConnMat[i*numel+j] = 0;
                    outConnMat[j*numel+i] = 1;}
                }
            }
            if(connType == 2){
                outConnMat[i*numel+j] = 1;
                outConnMat[j*numel+i] = 1;
            }
            
            
        }
        start++;
    }
    
    mxFree(tempMat);
    return;
}

/*Calculate connection Probability based on No of common neighborgs*/
double connProbability(int I,int J,double* connMat,int Numel, double* CBins, double* CProb, int connNumel){
    int i,commonNo, iC;
    commonNo=0;
    for(i=0;i<Numel;i++){
        if(connMat[I*Numel+i] == 1 && connMat[i*Numel+J] == 1)/*crap*/
            commonNo++;
    }
    
    iC=0;
    while((commonNo>CBins[iC]) && (iC<connNumel) ){iC++;}
    return CProb[iC];
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
        
            
    if( RND > CProb  ){
        return 0;}
    else if( RND > RProb[iR] * CProb){
        return 1;}
    else {return 2;}
}
/*
double max(double a, double b){
    if (a >= b)
        return a;
    else 
        return b;
}*/