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

void CopyArray(double A[],double B[], int n);
void BottomUpSort(int n, int A[], int B[], int C[]);
void BottomUpMerge(int A[], int iLeft, int iRight, int iEnd, int B[], int C[]);
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
    double *vect_Dist, *vectIDX, *t_vect_Dist;
    int steps, pntsPerTime;
    int bin, rpts, cntr, vectLen;
    int i, j,*x,*y, start,NoPts, numBins;
    int *vect_Ids, *vect_Conn;
    
    inputConns = mxGetPr(prhs[0]);/*pointer of connection Matrix*/
    NoPts = mxGetN(prhs[0]);
    
    pntsPerTime = (int)(*(mxGetPr(prhs[1])));/*points flipped per step*/
    
    bins = mxGetPr(prhs[2]);
    numBins = mxGetN(prhs[2]);
    
    probsS = mxGetPr(prhs[3]);
    probsR = mxGetPr(prhs[4]);
    
    steps = (int)(*(mxGetPr(prhs[5]))); /*total steps*/
    
    plhs[0] = mxCreateDoubleMatrix(NoPts,NoPts,mxREAL);
    connMat = mxGetPr(plhs[0]);
    
    plhs[1] = mxCreateDoubleMatrix(1,steps,mxREAL);
    improvementS = mxGetPr(plhs[1]);
    
    plhs[2] = mxCreateDoubleMatrix(1,steps,mxREAL);
    improvementR = mxGetPr(plhs[2]);
    
    plhs[3] = mxCreateDoubleMatrix(1,NoPts*(NoPts-1),mxREAL);
    vectIDX = mxGetPr(plhs[3]);
    
    /*histS = (double*)mxCalloc(numBins , sizeof(double));
    histR = (double*)mxCalloc(numBins , sizeof(double));
    
    x = (int*)mxCalloc(pntsPerTime , sizeof(int));
    y = (int*)mxCalloc(pntsPerTime , sizeof(int));*/
    
    vect_Dist = (double*)mxCalloc(NoPts*(NoPts-1) , sizeof(double));
    vect_Conn = (int*)mxCalloc(NoPts*(NoPts-1) , sizeof(int));
    
    
    printf("Running for steps=%d and points=%d\n", steps, pntsPerTime);
    
    /*copy to intermediate vectors for sorting later*/
    cntr = 0;
    start = 1;
    for(i=0;i<NoPts;i++){
        for(j=start;j<NoPts;j++){
            if(inputConns[i*NoPts+j] > 0){ /*if single connection*/
                vect_Dist[cntr] = inputConns[i*NoPts+j];
                if(inputConns[j*NoPts+i] > 0){ /*if reciprocal connection*/
                    vect_Conn[cntr] = 2;
                }else{
                    vect_Conn[cntr] = 1;
                }
                cntr++;
            }
        }
        start++;
    }
    
    vectLen = cntr ;
    mxRealloc(vect_Dist, vectLen * sizeof(double));
    mxRealloc(vect_Conn, vectLen * sizeof(int));
    vect_Ids = (int*)mxCalloc(vectLen , sizeof(int));
    
    /*perform mergesort based on distances*/
    t_vect_Dist = (double*)mxCalloc(pntsPerTime , sizeof(double));
    BottomUpSort(vectLen, vect_Dist, t_vect_Dist, vect_Ids);
    
    
    

    
    
    printf("Deallocate memmory\n");
    mxFree(t_vect_Dist);
    mxFree(vect_Dist);
    mxFree(vect_Conn);
    mxFree(vect_Ids);
    /*mxFree(histS);
    mxFree(histR);
    mxFree(x);
    mxFree(y);*/
    
    return;
}

/* array A[] has the items to sort; array B[] is a work array */
void BottomUpSort(int n, int A[], int B[], int C[])
{
  int width;
 
  /* Each 1-element run in A is already "sorted". */
 
  /* Make successively longer sorted runs of length 2, 4, 8, 16... until whole array is sorted. */
  for (width = 1; width < n; width = 2 * width)
    {
      int i;
 
      /* Array A is full of runs of length width. */
      for (i = 0; i < n; i = i + 2 * width)
        {
          /* Merge two runs: A[i:i+width-1] and A[i+width:i+2*width-1] to B[] */
          /* or copy A[i:n-1] to B[] ( if(i+width >= n) ) */
          BottomUpMerge(A, i, min(i+width, n), min(i+2*width, n), B, C);
        }
 
      /* Now work array B is full of runs of length 2*width. */
      /* Copy array B to array A for next iteration. */
      /* A more efficient implementation would swap the roles of A and B */
      CopyArray(A, B, n);
      /* Now array A is full of runs of length 2*width. */
    }
}
 
void BottomUpMerge(int A[], int iLeft, int iRight, int iEnd, int B[], int C[])
{
  int i0 = iLeft;
  int i1 = iRight;
  int j;
 
  /* While there are elements in the left or right lists */
  for (j = iLeft; j < iEnd; j++)
    {
      /* If left list head exists and is <= existing right list head */
      if (i0 < iRight && (i1 >= iEnd || A[i0] <= A[i1]))
        {
          B[j] = A[i0];
          C[j] = i0;
          i0 = i0 + 1;
        }
      else
        {
          B[j] = A[i1];
          C[j] = i1;
          i1 = i1 + 1;
        }
    }
}

void CopyArray(double A[],double B[], int n){
    int i=0;
    for(i=0;i<n;i++){
        A[i] = B[i];
    }
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