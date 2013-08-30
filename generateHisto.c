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
/*stefanou@imbb.forth.gr*/

double min(double a, double b);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    /*Declerations:*/
    double *inputConns, *inputDist ,*allCombs, *histMinVal, *histMaxVal, histRange, *histo, distSum;
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
    
    inputDist = mxGetPr(prhs[1]);/*pointer of distance matrix*/
    
    allCombs = mxGetPr(prhs[2]);/*pointer of distance matrix*/
    CombsM = mxGetM(prhs[2]);
    CombsN = mxGetN(prhs[2]);
    
    histMinVal = mxGetPr(prhs[3]); /*min intersoma distance*/
    histMaxVal = mxGetPr(prhs[4]);  /*max distance*/
    
    histRange = *histMaxVal - *histMinVal;
    /*printf("Histo range is: %f\n",histRange);*/
    
    numBins = (int)(histRange + 1); /*something like ceil..*/
    /*printf("No of bins is: %d\n",numBins);*/
    plhs[0] = mxCreateDoubleMatrix(1,numBins,mxREAL);
    histo = mxGetPr(plhs[0]);
    
    histoCtr = (int*)mxCalloc(numBins,sizeof(int));
    /*unsigned int *histo[numBins]; high resolution bin*/
    for(i=0;i<numBins;i++)
        histo[i]=0; /*initialize histogram to zero*/
    
    
    /*For all n cell clusters, measure the average intersomatic distance*/
    currCluster = (int*)mxCalloc(CombsM,sizeof(int));
    for(cluster=0;cluster<CombsN;cluster++){
        /*pick cluster*/
       /* printf("cluster: ");*/
        for(i=0;i<CombsM;i++){
            currCluster[i] = allCombs[cluster*CombsM+i];
            /*printf("%d ",currCluster[i]);*/
        }
       /* printf("\n");*/
        /*measure intersomatic distance*/
        start = 1;
        distSum = 0;
        uniqueConns = 0;
        for(i=0;i<CombsM;i++){
            for(j=start;j<CombsM;j++){
                distSum += inputDist[currCluster[i]*NoPts+currCluster[j]];
                uniqueConns++;
            }
            start++;
        }
        avrIntersomatic = distSum / uniqueConns; /*for this cluster*/
        /*printf("intersomatic: %f \n",avrIntersomatic);*/
        
        /*measure total connections in cluster*/
        totalConnections = 0;
        for(i=0;i<CombsM;i++){
            for(j=0;j<CombsM;j++){
                if(inputConns[currCluster[i]*NoPts+currCluster[j]]==1 ||
                        inputConns[currCluster[j]*NoPts+currCluster[i]]==1){
                    totalConnections++ ;
                }}}
        
        bin = -1;
        if((histRange * numBins)>0){
            bin = min( numBins-1 , (*histMaxVal - avrIntersomatic / histRange * numBins) );
            /*printf("Connections: %d \n",totalConnections);*/
            
        }
        if(bin>=0){
            /*printf("Bin: %d \n",bin);*/
            histoCtr[bin]++;
            histo[bin] += totalConnections;
        }
    }
    
    /*Average connections*/
    for(i=0;i<numBins;i++){
        histo[i] /= histoCtr[i];
    }
    
    
    
    
    /*Form histogram from connected pairs:*/
    /*totalConnections = 0;
     * for(i=0;i<NoPts;i++)
     * for(j=0;j<NoPts;j++)
     * if(inputConns[i*NoPts+j]>0){
     * bin = min( numBins-1 , (*histMaxVal - (inputDist[i*NoPts+j]) / histRange * numBins) );
     * histo[bin]++;
     * totalConnections++;
     * }
     */
    
    
    
    /*Normalize histogram to get probabilities:*/
    /* for(i=0;i<numBins;i++)
     * if(histo[i] > maxBins)
     * maxBins = histo[i];*/
    /*
     * for(i=0;i<numBins;i++)
     * histo[i] /= totalConnections;
     */
    mxFree(currCluster);
    return;
}

double min(double a, double b){
    if(a<b)
        return a;
    else
        return b;
}
