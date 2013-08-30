#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "mex.h"
#include "matrix.h"
#include "math.h"

# define DEBUG

/*stefanou@imbb.forth.gr*/

int decideConnection(double dist, double *CBins, double *CProb, int CNumel, double *RBins, double *RProb, int RNumel);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    /*Declerations:*/
    double *distMat, *connMat, *connBins, *connProb, *recipBins,  *recipProb;
    int numel , start, IR;
    int i,j,k, connType, connNumel, recipNumel;
    int incomming, outgoing;
	
    distMat = mxGetPr(prhs[0]);
    numel = (int)mxGetN(prhs[0]);
	/*Connection probabilities*/
	connBins = mxGetPr(prhs[1]);
	connProb = mxGetPr(prhs[2]);
    connNumel = (int)mxGetN(prhs[2]);
	/*reciprocal probabilities*/
	recipBins = mxGetPr(prhs[3]);
	recipProb = mxGetPr(prhs[4]);
    recipNumel = (int)mxGetN(prhs[4]);
    /*IR = (int)mxGetScalar(prhs[5]);*/
	
	plhs[0] = mxCreateDoubleMatrix(numel,numel,mxREAL);
    connMat = mxGetPr(plhs[0]);
	
	srand(time(NULL));
	/*For every unique pair*/
	start=1;
	for(i=0;i<numel;i++){
		for(j=start;j<numel;j++){
            /*restrict connectivity!*/
            incomming = 0;
            outgoing = 0;
            for(k=0;k<numel;k++){
                incomming += connMat[k*numel+j];
                outgoing  += connMat[i*numel+k];
            }
            /*If current pair can not project / receive skip it*/
            /*if ( (incomming>6) || (outgoing>5))
                continue;*/
            
            /*chose reciprocal probabilities independend or not!*/
           /* if(IR){
                connType = decideConnectionIR(distMat[i*numel+j],connBins,connProb,connNumel,recipBins,recipProb,recipNumel);
            }else{*/
                connType = decideConnection(distMat[i*numel+j],connBins,connProb,connNumel,recipBins,recipProb,recipNumel);
           /* }*/
            
            
			if(connType == 0){
				connMat[i*numel+j] = 0;
				connMat[j*numel+i] = 0;
			}
			if(connType == 1){
				if(((double)rand()/RAND_MAX) <= 0.5){
					connMat[i*numel+j] = 1;
					connMat[j*numel+i] = 0;}
				else{
					connMat[i*numel+j] = 0;
					connMat[j*numel+i] = 1;}
			}
			if(connType == 2){
				connMat[i*numel+j] = 1;
				connMat[j*numel+i] = 1;
			}
		}
		start++;
	}
	
    return;
}

int decideConnection(double dist, double *CBins, double *CProb, int CNumel, double *RBins, double *RProb, int RNumel){
	int iC, iR;
    double RND;
    RND = (double)rand()/RAND_MAX;

	iC=0;
	while((dist>CBins[iC]) && (iC<CNumel)){iC++;}
    
    iR=0;
	while((dist>RBins[iR]) && (iR<RNumel)){iR++;}
    
	if( RND > (CProb[iC]+RProb[iR]) ){
        return 0;
    }else if(RND > RProb[iR]){
         return 1;
    }else {
        return 2;
    }
}
/*
int decideConnectionIR(double dist, double *CBins, double *CProb, int CNumel, double *RBins, double *RProb, int RNumel){
	int iC, iR, Type;
	iC=0;
	while(dist>CBins[iC]){iC++;}
	if( ((double)rand()/RAND_MAX) > CProb[iC] ){
		Type = 0;}
    else{ Type = 1;}
    
	iR=0;
	while(dist>RBins[iR]){iR++;}
	if( ((double)rand()/RAND_MAX) <= RProb[iR] ) 
		{Type = 2;}
    
    return Type;

}
*/
        
/*int decideConnection(double dist, double *CBins, double *CProb, int CNumel, double *RBins, double *RProb, int RNumel){
	int iC, iR;
	iC=0;
	while(dist>CBins[iC]){iC++;}
	if( ((double)rand()/RAND_MAX) > CProb[iC] )
		return 0;
		
	iR=0;
	while(dist>RBins[iR]){iR++;}
	if( ((double)rand()/RAND_MAX) <= RProb[iR] ) 
		{return 2;}
	else
		{return 1;}
}*/