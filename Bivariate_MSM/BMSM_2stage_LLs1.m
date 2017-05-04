    function LLs=BMSM_2stage_LLs1(para1, data, kbar, n, A)
  
        [~, LLs1, LLs2] = BMSM_2stage_likelihood1(para1, kbar, data, A, n);
        
        LLs=LLs1+LLs2;
                
    end