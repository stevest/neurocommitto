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

int returnBin(double el, double *bins);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    /*Declerations:*/
    double *inputConns, *inputDist ,*allCombs, *histMinVal, *histMaxVal, histRange, *histSingle. *histReciprocal, distSum;
    double avrIntersomatic;
    int bin, maxBins;
    int i, j, start,uniqueConns,CombsM,CombsN,NoPts, numBins, totalConnections;
    int cluster, *currCluster,*histoCtr;
    maxBins = 0;
    /*Input: the 3d points' coordinates*/
    /*double *inputPts = mxGetPr(prhs[0]); //concatenated column-wise (Matlab does this..)*/
    /*int NoPts = mxGetN(prhs[0]);*/
    
    inputConns = mxGetPr(prhs[0]);/*pointer of connection Matrix*/
    NoPts = mxGetN(prhs[0]);

	bins = mxGetPr(prhs[1]);
	numBins = mxGetN(prhs[1]);
	
    plhs[0] = mxCreateDoubleMatrix(1,numBins,mxREAL);
    histSingle = mxGetPr(plhs[0]);
	
	plhs[1] = mxCreateDoubleMatrix(1,numBins,mxREAL);
    histReciprocal = mxGetPr(plhs[1]);
    
    /*unsigned int *histo[numBins]; high resolution bin*/
   
        start = 1;
        for(i=0;i<NoPts;i++){
            for(j=start;j<NoPts;j++){
				int bin=0;
				if(inputConns[i*NoPts+j] > 0){
					bin = returnBin(inputConns[i*NoPts+j],bins);
					histSingle[bin]++;
					if(inputConns[j*NoPts+i] > 0)){
						bin = returnBin(inputConns[j*NoPts+i],bins);
						histReciprocal[bin]++;
					}
				}
            }
            start++;
        }
    return;
}


int returnBin(double el, double *bins){
	int c=0;
	while(el<bins[c]) c++;
	return c;
}