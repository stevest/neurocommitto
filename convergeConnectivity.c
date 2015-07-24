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

int compare_function(const void *name1, const void *name2);
double flip(double x);
int returnBin(double el, double *bins);
double L2Norm(double *vector, int len);
void normHisto2one(double *histo, int len);
void simpleHistDiff(double* hist1, double* hist2, double* diff, int len);
double absolut(double x);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	/*Declerations:*/
	double *inputConns ,*bins, *probsS, *probsR, *histS, *histR, *histC;
	double *connMat, cvgS, cvgR, t_cvgS, t_cvgR;
	double *improvementS, *improvementR, *inputDist;
	double *vect_Dist, **p_vect_Dist,*vect_ptr;
	double* diffS, *diffR, rndf;
	int steps, pntsPerTime, *org_IDX, rndi, *vectIDX;
	int bin, cbin,rpts, cntr, vectLen;
	int i, j,ic,jc,x,t, start,NoPts, numBins;
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

	/*plhs[1] = mxCreateDoubleMatrix(1,steps,mxREAL);
	improvementS = mxGetPr(plhs[1]);

	plhs[2] = mxCreateDoubleMatrix(1,steps,mxREAL);
	improvementR = mxGetPr(plhs[2]);

	plhs[3] = mxCreateDoubleMatrix(1,NoPts*(NoPts-1),mxREAL);
	vectIDX = mxGetPr(plhs[3]);*/

	org_IDX = (int*)mxCalloc(NoPts*(NoPts-1)/2 , sizeof(int));
	
	histS = (double*)mxCalloc(numBins , sizeof(double));
	histR = (double*)mxCalloc(numBins , sizeof(double));
	histC = (double*)mxCalloc(numBins , sizeof(double));

	diffS = (double*)mxCalloc(numBins , sizeof(double));
	diffR = (double*)mxCalloc(numBins , sizeof(double));

	/*x = (int*)mxCalloc(pntsPerTime , sizeof(int));
	y = (int*)mxCalloc(pntsPerTime , sizeof(int));*/

	vect_Dist = (double*)mxCalloc(NoPts*(NoPts-1)/2 , sizeof(double));
	vect_Conn = (int*)mxCalloc(NoPts*(NoPts-1)/2 , sizeof(int));


	/*printf("Running for steps=%d and points=%d\n", steps, pntsPerTime);*/

	/*copy to intermediate vectors for sorting later*/
	/*UNSORTED DATA*/
	cntr = 0;
	start = 1;
	for(i=0;i<NoPts;i++){
		for(j=start;j<NoPts;j++){
			if((inputConns[i*NoPts+j] > 0) && (inputConns[j*NoPts+i] > 0)){ /*if reciprocal connection*/
				vect_Dist[cntr] = inputConns[i*NoPts+j];
				org_IDX[cntr] = i*NoPts+j; /*upper diagonal index*/
				vect_Conn[cntr] = 2;
				cntr++;}
			else if((inputConns[i*NoPts+j] > 0) || (inputConns[j*NoPts+i] > 0)){
				vect_Dist[cntr] = inputConns[i*NoPts+j];
				org_IDX[cntr] = i*NoPts+j;
				vect_Conn[cntr] = 1;
				cntr++;}
		}
		start++;
	}

	vectLen = cntr ;
	vectIDX = (int*)mxCalloc(vectLen , sizeof(int));
	/*mxRealloc(vect_Dist, vectLen * sizeof(double));
	mxRealloc(vect_Conn, vectLen * sizeof(int));
	mxRealloc(org_IDX, vectLen * sizeof(int));*/
	/*vect_Ids = (int*)mxCalloc(vectLen , sizeof(int));*/

	/*perform sort based on distances*/
	p_vect_Dist = (double**)mxCalloc(vectLen , sizeof(double*));
	vect_ptr = (double*)vect_Dist;
	for(i=0;i<vectLen;i++){
		p_vect_Dist[i] = vect_ptr+i;/*kanw swsta to incrementing?*/
	}
	qsort(p_vect_Dist, vectLen, sizeof(double*), compare_function);

	/*transform pointers to indices*/
	for(i=0;i<vectLen;i++){
		vectIDX[i] =  (p_vect_Dist[i] - vect_Dist) ;
	}



	srand(time(NULL));
	/*Begin iterative process*/
	for(t=0;t<steps;t++){
		/*Clear histos*/
		memset(histS, 0, numBins*sizeof(double));
		memset(histR, 0, numBins*sizeof(double));
		memset(diffS, 0, numBins*sizeof(double));
		memset(diffR, 0, numBins*sizeof(double));
		memset(histC, 0, numBins*sizeof(double));

		/*Create histo using sorted elements (ascending!) */
		cbin = 0;
		for(i=0;i<vectLen;i++){
			if(vect_Dist[vectIDX[i]] > 0){
				if(vect_Dist[vectIDX[i]] > bins[cbin]){
					cbin++;
					i--;
					continue;
				}
				if(vect_Conn[i] == 2){
					histR[cbin]++;}
				if(vect_Conn[i] > 0){
					histS[cbin]++;}
			}
		}

		/*Create cumulative histo*/
		histC[0] = histS[0];
		for(i=1;i<numBins;i++){
			histC[i] = histC[i-1] + histS[i];
		}
		/*compare histograms with reference ones*/
		simpleHistDiff(probsS, histS, diffS, numBins);
		simpleHistDiff(probsR, histR, diffR, numBins);
		/*the score must be normalized per bin (summing to one), to be used as probability*/
		normHisto2one(diffS, numBins);/*cumulative probability!*/
		normHisto2one(diffR, numBins);
		/*pick bins to remove from, based in above probability*/
		for(rpts=0;rpts<pntsPerTime;rpts++){
			rndf = ((double)rand()/RAND_MAX);
			x = numBins-1;/*bin idx to remove from*/
			if(rpts%2){
				while( (rndf < diffS[x]) && (x>=0) ){x--;}
				x++;
				/*Mark points to be removed by setting negative distance*/
				/*do{*/
				rndi = (int)((((double)rand()/RAND_MAX)*(histS[x] - histC[x-1])) + histC[x-1] );
				/*}while(vect_Conn[vectIDX[rndi]] == 0);*/
			}else{
				while( (rndf < diffS[x]) && (x>=0) ){x--;}
				x++;
				/*Mark points to be removed by setting negative distance*/
				/*do{*/
				rndi = (int)((((double)rand()/RAND_MAX)*(histS[x] - histC[x-1])) + histC[x-1] );
			}
			
				/*remove random connections using above probability*/
				if( vect_Conn[vectIDX[rndi]] == 2  ){
					vect_Conn[vectIDX[rndi]] = 1;
				}else{
					vect_Conn[vectIDX[rndi]] = 0;/*double set it*/
					vect_Dist[vectIDX[rndi]] = -1;
				}
		}

	}

	/*save to out connection matrix and quit*/
	start = 1;
	for(i=0;i<vectLen;i++){
		if( vect_Dist[i] > 0){ 
			if( vect_Conn[i] == 2) {
				connMat[org_IDX[i]] = 1;
				/*inverse the index to column-wise*/
				ic = (int)(org_IDX[i] % NoPts);
				jc = (int)((org_IDX[i]-ic) / NoPts);
				connMat[ic*NoPts+jc] = 1;
			}else if(vect_Conn[i] == 1){
				if( 0.5 > ((double)rand()/RAND_MAX) ){
					connMat[org_IDX[i]] = 1;}
				else{
					/*inverse the index to column-wise*/
					ic = (int)(org_IDX[i] % NoPts);
					jc = (int)((org_IDX[i]-ic) / NoPts);
					connMat[ic*NoPts+jc] = 1;}
			}
		}
	}

	/* printf("Deallocate memmory\n");*/
	mxFree(p_vect_Dist);
	mxFree(vect_Dist);
	mxFree(vect_Conn);
	/*mxFree(vect_Ids);*/
	mxFree(histS);
	mxFree(histR);
	mxFree(histC);
	mxFree(diffS);
	mxFree(diffR);
	mxFree(org_IDX);
	mxFree(vectIDX);
	/*mxFree(x);
	mxFree(y);*/

	return;
}


int compare_function(const void *name1, const void *name2)
{
	double name1_ = *(*((const double **)name1));
	double name2_ = *(*((const double **)name2));
	if(name1_ < name2_){
		return -1;}
	if(name1_ > name2_){
		return 1;}
	if(name1_ == name2_){
		return 0; }
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

void normHisto2one(double *histo, int len){
	int i=0;
	double cumulate=0;
	double factor, l2norm;

	l2norm = L2Norm(histo,len);

	for (i=0;i<len;i++){
		histo[i]=histo[i]/l2norm;
		cumulate += histo[i];
	}

	factor = 1 / cumulate;

	for (i=0;i<len;i++){
		histo[i]=histo[i]*factor;
	}
	/*cumulative probability*/
	for (i=1;i<len;i++){
		histo[i] += histo[i-1];
	}
	return;
}

void simpleHistDiff(double* hist1, double* hist2, double* diff, int len){
	int i=0;
	for(i=0;i<len;i++){
		diff[i] = absolut((hist1[i] - hist2[i]));
	}
	return;
}

double absolut(double x){
	if(x<0){
		return -x;
	}else{
		return x;
	}
}