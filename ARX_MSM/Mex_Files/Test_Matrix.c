#include <math.h>
#include <limits.h>
#include "mex.h"
#include "matrix.h"


/* To compile in matlab, type mex ARX_MSM_2.c */

void mat_times(double *pmat, double *A, double *pout,size_t rowp, size_t colp,size_t rowA, size_t colA){

    size_t tr, tc;

    for (tc = 0; tc < colA; tc++)
    {
        pout[tc]=0.0000;
        for (tr = 0; tr < rowA; tr++)
        {
            pout[tc] += A[tc*rowA + tr]*pmat[tr];       
        }

    }
        for (tc = 0; tc < rowA; tc++)
    {
        for (tr = 0; tr < colA; tr++)
        {
            pout[tc]=A[tc*rowA+tr];
            mexPrintf("T=%f\n", A[tc*rowA+tr]);    
        }

    }
    
}

/* Begin the Gateway function  */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]){

    double *pmat, *A;           /*The inputs */
    double *pout;                  /*The outputs */ 
    size_t rowp, rowA, colp,colA;

  /*  Create a pointer to the input matrices . */
   pmat  = mxGetPr(prhs[0]);
   A     = mxGetPr(prhs[1]);
 	
 /*  Get the dimensions of the input matricec. */
   rowp = mxGetM(prhs[0]); 
   colp = mxGetN(prhs[0]);
   rowA = mxGetM(prhs[1]); 
   colA = mxGetN(prhs[1]);

   /*  Set the output pointer to the output matrix. */
    plhs[0] = mxCreateDoubleMatrix((mwSize)rowp, (mwSize)colA, mxREAL);
      
    
    /*  Create a C pointer to a copy of the output matrix. */
    pout = mxGetPr(plhs[0]);
       
    /*  Call the C subroutine. */
    mat_times(pmat, A,pout,rowp,colp,rowA,colA);
}