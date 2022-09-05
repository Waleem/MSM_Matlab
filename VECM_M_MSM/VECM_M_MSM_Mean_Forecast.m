function [E_rs,E_rf]=VECM_M_MSM_Mean_Forecast(para,p1,p2,lnS,lnF,kbar)
 % This is a helper function used to compute E(Rs_t+1) and E(Rf_t+1) for VECM(M) - MSM.
 % Input:   para    -   A 1-by-12 vector of VECM(M) - MSM parameter estimates:
 %                         [m01, m02, sigma1,sigma2, gamma_k, b, bs, bf, as, af,  rho_e, lambda]
 %          p1      -   A kbar^2-by-1 state probability matrix forecast for spot equation.
 %          p2      -   A kbar^2-by-1 state probability matrix forecast for futures equation.
 %          lnS     -   Scalar - log-spot price.
 %          lnF     -   Scalar - log-futures price.
 %          KBAR    -   The number of volatility components.

 % Output:  
 %          E_rs    -   E_t[R_st+1] = E_t[ a_s((Ms_1 -1 + Ms_2 -1 +....+Ms_kbar-1))*(ln(St - lnFt))] - Scalar
 %          E_rf    -   E_t[R_ft+1] = E_t[ a_f((Ms_1 -1 + Ms_2 -1 +....+Ms_kbar-1))*(ln(St - lnFt))] - Scalar

    
    %Get BMSM state matrix for each series
    Mmat1 = Univ_MSM_states(para(1),kbar);
    Mmat2 = Univ_MSM_states(para(2),kbar);
    
    %Get the marginal probability for each M_k i.e prob[M_k = m01]
    ph_s = MSM_marginals(p1,para(1),Mmat1,kbar);
    ph_f = MSM_marginals(p2,para(2),Mmat2,kbar);
    
    %Get expected value of Mk, for each Mk
    em_s = para(1)*ph_s + (2-para(1))*(1-ph_s);
    em_f = para(2)*ph_f + (2-para(2))*(1-ph_f);
    
    as_m = para(9)*sum(em_s - 1,2);
    af_m = para(10)*sum(em_f - 1,2);
    
    % Get the T-by-1 error correction terms from the VECM model
    E_rs = para(7) + as_m.*(lnS-lnF);
    E_rf = para(8) + af_m.*(lnS-lnF);




end