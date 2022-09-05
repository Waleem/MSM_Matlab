#include <math.h>
#include <limits.h>
#include "mex.h"
#include "matrix.h"


/* To compile in matlab, type mex ARX_MSM_2.c */

void msmcore(double *pmat, double *A, double *y, size_t T, size_t k, double *pout, double *LLs, double *C, double *pa, double *pt, double *f){

    size_t t, tr, tk, offSet;
    int i=1;
    int zero=0;

	 for (tk = 0; tk < k; tk++)
	 {
		 pt[tk] = pmat[tk*(T+1)+0];  
	
     }

    for (t = 1; t < (T+1); t++)
    {
        offSet = t-i;
		f[0] =0.0000;
        
        for (tk = 0; tk < k; tk++)
        {
			pa[tk]=0.0000;
            for(tr=0; tr<k; tr++)
            {
                pa[tk] += A[tk*k+tr]*pt[tr]; /*Code works up to this line*/  
				/*mexPrintf("pa[k]=%f\n", pa[tk]);*/
            }
        }
        
        for(tr = 0; tr<k; tr++)
        {
         C[tr] = y[tr*T + (t-1)] * pa[tr]; /*Code works up to this line*/ 
         f[0] += C[tr];
        }
        
        if(f[0]==0)
           pt[0]=1;
        else
            for(tk = 0; tk<k; tk++)
			{
				pt[tk]             = C[tk]/f[0];
                pout[tk*(T+1) + t] = C[tk]/f[0];
			}
                   
       LLs[t-1]=log(f[0]);
        
    }
    
}

/* Begin the Gateway function  */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]){

    double *pmat, *A, *y;           /*The inputs */
    double *pout, *LLs, *C, *pa, *pt, *f;                  /*The outputs */ 
    size_t T, k;

  /*  Create a pointer to the input matrices . */
   pmat  = mxGetPr(prhs[0]);
   A     = mxGetPr(prhs[1]);
   y     = mxGetPr(prhs[2]);
	
 /*  Get the dimensions of the input matricec. */
   T = mxGetM(prhs[2]); 
   k = mxGetN(prhs[2]);

   /*  Set the output pointer to the output matrix. */
    plhs[0] = mxCreateDoubleMatrix((mwSize)(T+1), (mwSize)k, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, (mwSize)T, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1, (mwSize)k, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1, (mwSize)k, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(1, (mwSize)k, mxREAL);
	plhs[5] = mxCreateDoubleMatrix(1, 1, mxREAL);
    
    /*  Create a C pointer to a copy of the output matrix. */
    pout = mxGetPr(plhs[0]);
    LLs  = mxGetPr(plhs[1]);
    C    = mxGetPr(plhs[2]);
    pa   = mxGetPr(plhs[3]);
    pt   = mxGetPr(plhs[4]);
	f    = mxGetPr(plhs[5]);
    
    /*  Call the C subroutine. */
    msmcore(pmat, A, y, T, k, pout,LLs, C, pa, pt,f);
}