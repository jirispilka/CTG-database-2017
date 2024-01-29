/* fastmovingaverage.c calculate moving average, fast

   Compile: mex fastmovingaverage.c
   Syntax: y = fastmovingaverage(x, w);

*/

/*************************************************************************/
/*
% Example
% Data
T = 250; N = 100; x = cumsum(randn(T, N));
y = fastmovingaverage(x, ceil(T / 2));
*/
/*************************************************************************/
/*************************************************************************/

/*************************************************************************/
/*                                                                       */
/*            Author: Francesco Pozzi                                    */
/*            E-Mail: francesco.pozzi@anu.edu.au                         */
/*            Date: 11 December 2013                                     */
/*                                                                       */
/*************************************************************************/

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    mwSignedIndex hh, ii;            /*  counters                       */
    mwSize m, n, m2;                 /*  size of matrices               */
    double *x, *y;                   /*  input and output matrices      */
    mwSize w;                        /*  third-order correlations       */
    /*  matrices needed */
    mxArray *xArray, *yArray;

/*  Now we need to get the data */
    if (nrhs == 2) {
        x = mxGetPr(prhs[0]);
        m = (mwSize) mxGetM(prhs[0]);
        n = (mwSize) mxGetN(prhs[0]);
        w = (mwSize) mxGetScalar(prhs[1]);
        m2 = m - w + 1;
    }

/*  Then build the matrices for the output */
    yArray = mxCreateDoubleMatrix(m2, n, mxREAL);
    y = mxGetPr(yArray);
    plhs[0] = yArray;

/* Calculate first element of each column */
    for(ii=0; ii<n; ii++) {
        for(hh=0; hh<w; hh++) {
            y[ii * m2] += x[ii * m + hh];
        }
    }

/* Calculate all other elements for each column */
    for(ii=0; ii<n; ii++) {
        for(hh=1; hh<m2; hh++) {
            y[ii * m2 + hh] = y[ii * m2 + hh - 1] - x[ii * m + hh - 1] + x[ii * m + hh - 1 + w];
        }
    }
/* Divide by length of window */
    for(ii=0; ii<n; ii++) {
        for(hh=0; hh<m2; hh++) {
            y[ii * m2 + hh] /= w;
        }
    }
}
