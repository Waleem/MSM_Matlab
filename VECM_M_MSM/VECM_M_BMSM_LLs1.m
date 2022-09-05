    function LLs=VECM_M_BMSM_LLs1(para1, kbar, r, S, F, A,n)
  
        [LL,LLs1, LLs2] =VECM_M_BMSM_likelihood1(para1, kbar, r, S, F, A,n);
        
        LLs=LLs1+LLs2;
                
    end