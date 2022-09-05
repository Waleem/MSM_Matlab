function [ect,as_m, af_m,X1,X2] = VECM_M_BMSM_ect(para,r,S,F,kbar,n)
    % This algorithm computes the error term matrix and the time varying speeds of 
    % adjustment from the Vector Error Correction Bivariate Markov Switching Multifractal Model 
    % with Time Varying Speed of Long-Run Adjustment (VECM(M)-MSM).
    % VECM(M) Bivariate MSM:
    %                                         R_s = b_s + a_s((Ms_1 -1 + Ms_2 -1 +. ...+ Ms_kbar-1))*(S-F) + e_s
    %                                         R_f = b_f + a_f((Mf_1 -1 + Mf_2 -1 +. ...+ Mf_kbar-1)-1)*(S-F) + e_f
    %                                         e_s = sigma_s*(Ms_1*Ms_2*. ...*Ms_kbar)
    %                                         e_f = sigma_f*(Mf_1*Mf_2*. ...*Mf_kbar)
    
    % Input:  para     = 10-by-1 vector of parameters [m01, m02, sigma1,sigma2, gamma_k, b, bs, bf, as, af]
    %         kbar        -  Number of volatility components (scalar)
    %         r           -  T-by-2 matrix of returns [log_spot_return log_futures_return]
    %         S           -  T-by-1 vector of log_spot prices
    %         F           -  T-by-1 vector of log_futures prices
    %         n           -  The number of trading periods in 1 year. E.g
    %                        n= 252 for the number of business days in a year. n=12 for monthly data

    % Output: ect         -  T-by-2 matrix of residuals matrix for the spot and futures returns, respectively
    %         as_m        -  T-by-1 vector of time varying speed of adjustment for spot returns
    %         as_f        -  T-by-1 vector of time varying speed of adjustment for futures returns
   
    
    T=size(r,1);
    %Get probability matrix
    A_template = T_mat_template(kbar);
    [~, ~, ~, p1, p2,A] = VECM_M_BMSM_likelihood1(para, kbar, r, S, F, A_template,n);
    
    %Get BMSM state matrix for each series
    Mmat1 = Univ_MSM_states(para(1),kbar);
    Mmat2 = Univ_MSM_states(para(2),kbar);
    
    %Get the marginal probability for each M_k i.e prob[M_k = m01]
    ph_s = MSM_marginals(p1(2:end,:),para(1),Mmat1,kbar);
    ph_f = MSM_marginals(p2(2:end,:),para(2),Mmat2,kbar);
    
    %Get expected value of Mk, for each Mk
    em_s = para(1)*ph_s + (2-para(1))*(1-ph_s);
    em_f = para(2)*ph_f + (2-para(2))*(1-ph_f);
    
    as_m = para(9)*sum(em_s - 1,2);
    af_m = para(10)*sum(em_f - 1,2);
    
    % Get the T-by-1 error correction terms from the VECM model
    X1 = para(7)*ones(T,1) + as_m.*(S-F);
    X2 = para(8)*ones(T,1) + af_m.*(S-F);
    % Remove VECM component to get the residual terms matrix i.e e_s and e_f.
    ect = [r(:,1) - X1  r(:,2) - X2];
    
   % Now, we need to get conditional mean estimates, using the smoothed probability.
   % Just repeat the same process but with smoothed probability.
   
   
    % Get smoothed probability
    p1=kim_smooth(A,p1);
    p2=kim_smooth(A,p2);
    
    %Get the marginal probability for each M_k i.e prob[M_k = m01]
    ph_s = MSM_marginals(p1(2:end,:),para(1),Mmat1,kbar);
    ph_f = MSM_marginals(p2(2:end,:),para(2),Mmat2,kbar);
    
    %Get expected value of Mk, for each Mk
    em_s = para(1)*ph_s + (2-para(1))*(1-ph_s);
    em_f = para(2)*ph_f + (2-para(2))*(1-ph_f);
    
    as_m = para(9)*sum(em_s - 1,2);
    af_m = para(10)*sum(em_f - 1,2);
    
    % Get the T-by-1 conditional mean terms for spot and futures equation, from the VECM model
    X1 = para(7)*ones(T,1) + as_m.*(S-F);
    X2 = para(8)*ones(T,1) + af_m.*(S-F);
%----------------------------------------------------------------------------------------------------    
    
    
    function A = T_mat_template(kbar)  
        A = zeros((2^kbar),(2^kbar));
        for i =0:2^kbar-1       
            for j = i:(2^kbar-1)-i  
                A(i+1,j+1) = bitxor(i,j);
            end
        end
    end

end