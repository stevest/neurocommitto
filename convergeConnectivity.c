#include <stdio.h>
#include <stdlib.h>


#include "mex.h"
#include "matrix.h"
#include "math.h"

/*Function to analyze a connected network of neurons*/
/*Algorithm performs the below:*/
/*1.Extract all possible pairs from population: n*(n-1) */
/*2.Extract eucledean distances between pairs*/
/*3.Create a histogram of frequency of connections in a distance range*/
/*4.Normalize the histogram to get propabilities as a function of intersomatic distance*/
/*stamatiad.st@gmail.com*/

#ifndef _MSC_VER
#define min(a,b) a<b?a:b
#endif


double flip(double x);
int returnBin(double el, double *bins);
double L2Norm(double *vector, int len);
void normalizeVector(double *vector, int len);
double simpleHistDiff(double* hist1, double* hist2, int len);
double absolut(double x);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    /*Declerations:*/
    double *inputConns ,*bins, *probsS, *probsR, *histS, *histR;
    double *connMat, cvgS, cvgR, t_cvgS, t_cvgR;
    double *improvementS, *improvementR, *inputDist;
    int steps, pntsPerTime;
    int bin, rpts, cntr;
    int i, j,*x,*y, start,NoPts, numBins;
    
    inputConns = mxGetPr(prhs[0]);/*pointer of connection Matrix*/
    NoPts = mxGetN(prhs[0]);
    
    inputDist = mxGetPr(prhs[1]);/*pointer of distance Matrix*/
    
    pntsPerTime = (int)(*(mxGetPr(prhs[2])));/*points flipped per step*/
    
    bins = mxGetPr(prhs[3]);
    numBins = mxGetN(prhs[3]);
    
    probsS = mxGetPr(prhs[4]);
    probsR = mxGetPr(prhs[5]);
    
    steps = (int)(*(mxGetPr(prhs[6]))); /*total steps*/
    
    plhs[0] = mxCreateDoubleMatrix(NoPts,NoPts,mxREAL);
    connMat = mxGetPr(plhs[0]);
    
    plhs[1] = mxCreateDoubleMatrix(1,steps,mxREAL);
    improvementS = mxGetPr(plhs[1]);
    
    plhs[2] = mxCreateDoubleMatrix(1,steps,mxREAL);
    improvementR = mxGetPr(plhs[2]);
        
    histS = (double*)mxCalloc(numBins , sizeof(double));
    histR = (double*)mxCalloc(numBins , sizeof(double));
    
    x = (double*)mxCalloc(pntsPerTime , sizeof(int));
    y = (double*)mxCalloc(pntsPerTime , sizeof(int));
    
    printf("Running for steps=%d and points=%d\n", steps, pntsPerTime);
    
    /*copy to output*/
    for(i=0;i<NoPts*NoPts;i++){
            connMat[i] = inputConns[i];
    }
    
    /*Histogram function*/
    start = 1;
    for(i=0;i<NoPts;i++){
        for(j=start;j<NoPts;j++){
            int bin=0;
            if(connMat[i*NoPts+j] > 0){
                bin = returnBin(inputDist[i*NoPts+j],bins);
                histS[bin]++;
                if(connMat[j*NoPts+i] > 0){
                    histR[bin]++;
                }
            }
        }
        start++;
    }
    
    normalizeVector(histS, numBins);
    normalizeVector(histR, numBins);
    
    cvgS = simpleHistDiff(probsS, histS, numBins);
    cvgR = simpleHistDiff(probsR, histR, numBins);
    
    /*begin iterative convergence*/
    srand(time(NULL));
    cntr = 0;
    
    while(cntr < steps){
       /* printf("@%d\n", cntr);*/
        
        /*flip random connections*/
        printf("@%d resetting rnd points\n", cntr);
        memset(x, 0, pntsPerTime*sizeof(int));
        memset(y, 0, pntsPerTime*sizeof(int));
        printf("@%d Picking points\n", cntr);
        for(rpts=0;rpts<pntsPerTime;rpts++){
            y[rpts] = (int)(((double)rand()/RAND_MAX) * NoPts);
            x[rpts] = (int)(((double)rand()/RAND_MAX) * NoPts);
            if(y[rpts]>=NoPts){y[rpts]=NoPts-1;}
            if(x[rpts]>=NoPts){x[rpts]=NoPts-1;}
            connMat[y[rpts]*NoPts+x[rpts]] = flip(connMat[y[rpts]*NoPts+x[rpts]]);
        }
        
        /*reset histogram*/
        printf("@%d Resetting histos\n", cntr);
        memset(histS, 0, numBins*sizeof(double));
        memset(histR, 0, numBins*sizeof(double));
        
        /*Histogram function*/
        printf("@%d Calculating histo\n", cntr);
        start = 1;
        for(i=0;i<NoPts;i++){
            for(j=start;j<NoPts;j++){
                int bin=0;
                if(connMat[i*NoPts+j] > 0){
                    bin = returnBin(inputDist[i*NoPts+j],bins);
                    histS[bin]++;
                    if(connMat[j*NoPts+i] > 0){
                        histR[bin]++;
                    }
                }
            }
            start++;
        }
        printf("@%d Normalizing vectors\n", cntr);
        normalizeVector(histS, numBins);
        normalizeVector(histR, numBins);
        printf("@%d Computing histo differences\n", cntr);
        t_cvgS = simpleHistDiff(probsS, histS, numBins);
        t_cvgR = simpleHistDiff(probsR, histR, numBins);
        
        if( (t_cvgS > cvgS) || (t_cvgR > cvgR) ){
            for(rpts=0;rpts<pntsPerTime;rpts++){
                connMat[y[rpts]*NoPts+x[rpts]] = flip(connMat[y[rpts]*NoPts+x[rpts]]);
            }
        }else{
            printf("@%d We have improved!\n", cntr);
            improvementS[cntr] = cvgS - t_cvgS;
            improvementR[cntr] = cvgR - t_cvgR;
            
            cvgS = t_cvgS;
            cvgR = t_cvgR;
            
            cntr++;
        }
    } /*END for times*/
    
    
    printf("Deallocate memmory\n");
    mxFree(histS);
    mxFree(histR);
    mxFree(x);
    mxFree(y);
    
    return;
}


double flip(double x){
    if (x > 0){
        return 0;
    }else{
        return 1;
    }
}


int returnBin(double el, double *bins){
    int c=0;
    while (el > (bins[c]) ) {c++;}
    return c;
}

double L2Norm(double *vector, int len){
    int i=0;
    double cumulate=0;
    for(i=0;i<len;i++){
        cumulate+= vector[i] * vector[i];
    }
    return sqrt(cumulate);
}

void normalizeVector(double *vector, int len){
    int i=0;
    double l2norm;
    
    l2norm = L2Norm(vector,len);
    
    for (i=0;i<len;i++){
        vector[i]=vector[i]/l2norm;
    }
    return;
}

double simpleHistDiff(double* hist1, double* hist2, int len){
    int i=0;
    double cumulate=0;
    for(i=0;i<len;i++){
        cumulate += (hist1[i] - hist2[i]);
    }
    return absolut(cumulate);
}

double absolut(double x){
    if(x<0){
        return -x;
    }else{
        return x;
    }
}