function [LL, LLs, A, p_mat, M, g_m, et] = MSM2_likelihood(para, kbar, et,n, para2)
       

if length(para)==1
    para=[para para2];
end
%Gather MSM paremeters
b          = para(1);
m0         = para(2);
gamma_kbar = para(3);
sigma      = para(4)/sqrt(n);

k2 =2^kbar;

A = ARX_MSM_Transition_Mat(b, gamma_kbar, kbar);
% Get state vector:
[M, ~,g_m] = Univ_MSM_states(m0,kbar);

T = size(et,1); 
LLs = zeros(1,T);

p_mat      = zeros(T+1,k2);
p_mat(1,:) = (1/k2)*ones(1,k2);


%*----------------------------------------------------------------------*
%*                        Likelihood algorithm  begins                        *
%*----------------------------------------------------------------------*

s   = repmat(sigma*g_m,T,1);
w_t = repmat(et,1,k2);
w_t = ((2*pi)^-0.5) * exp( - 0.5.*((w_t./s).^2))./s; 
w_t = w_t + 1e-16;
         

pmat=p_mat;   
%This is a mex function. Please compile before you use. Else, uncomment and
%use the loop below.
[p_mat, LLs]=ARX_MSM_Core(pmat,A,w_t);  
% for t=2:T+1          
%     piA = (p_mat(t-1,:)*A);
%     C   = (w_t(t-1,:).*piA); 
%     ft  = sum(C);
%     if ft == 0                      %This stop div by zero if probs are too low
%         p_mat(t,1) = 1;   
%     else
%         p_mat(t,:) = C / ft; 
%     end
%     
%     LLs(t-1) = log(dot(w_t(t-1,:),piA));
% end

LL = -sum(LLs);
if ~ isfinite(LL)
    disp('Log-likelihood is inf. Probably due to all zeros in pi_mat.')
end
end
