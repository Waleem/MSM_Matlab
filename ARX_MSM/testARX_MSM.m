function [LLs,p_mat] = testARX_MSM(p_mat, w_t,A)

T=size(w_t,1);
for t=2:T+1   
    
    piA = (p_mat(t-1,:)*A);
    C   = (w_t(t-1,:).*piA); 
    ft  = sum(C);
    if ft == 0                      %This stop div by zero if probs are too low
        p_mat(t,1) = 1;   
    else
        p_mat(t,:) = C / ft; 
    end
    
 LLs(t-1) = log(dot(w_t(t-1,:),piA));
end
