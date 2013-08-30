#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mex.h"
#include "matrix.h"
#include "math.h"

/*As in  Duncan J. Watts and Steven Strogatz, 1998 (wikipedia)*/
/*stefanou@imbb.forth.gr*/


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    /*Declerations:*/
    double *connMat, *clustCoeff;
    int numel , nodesNo, linksNo, start, *nodesVect;
    int i,j, n,c, counter;
    
    connMat = mxGetPr(prhs[0]);
    numel = (int)mxGetN(prhs[0]);

    plhs[0] = mxCreateDoubleMatrix(numel,1,mxREAL);
    clustCoeff = mxGetPr(plhs[0]);
    
    nodesVect = (int*)mxCalloc(numel,sizeof(int));

    /*For every neuron calculate culster coefficient*/
    /*(the current vertex is not checked for directed connectivity; its neighbours yes)*/
    for(n=0;n<numel;n++){
        nodesNo = 0;
        linksNo = 0;
        memset(nodesVect, 0, sizeof(int)*numel );
        for(c=0;c<numel;c++){
            if (connMat[n*numel+c] || connMat[n+numel*c]) {/*incoming or outgoing*/
                nodesVect[nodesNo] = c;
                nodesNo++;
            }
        }
        /*printf("NodesNo = %d\n",nodesNo);*/
        /*extract possible pairs between neighbours*/
        /*(neighbours are checked for directed connectivity; the current vertex not)*/
        start = 1;
        for(i=0;i<nodesNo;i++){
            for(j=start;j<nodesNo;j++){
                if( connMat[nodesVect[i]*numel+nodesVect[j]] )
                        linksNo++; 
                if( connMat[nodesVect[j]*numel+nodesVect[i]] )
                        linksNo++;
            }
            start++;
        }
        /*printf("linksNo = %d\n",linksNo);*/
        /*from Barabasi Network biology understands the cell's functional organization,2004*/
         clustCoeff[n] = (double)linksNo / ( (double)nodesNo*((double)nodesNo-1) );
    }
    
    mxFree(nodesVect);
    return;
}


