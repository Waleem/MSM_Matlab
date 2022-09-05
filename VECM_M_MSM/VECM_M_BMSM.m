function [out,L,states,diagnostics] = VECM_M_BMSM(r,S,F,kbar,n,startingvals, options)
    % -------------------------------------------------------------------------
    %         Vector Error Correction Bivariate Markov Switching Multifractal Model with Time Varying Speed 
    %          with Time Varying Speed of Long-Run Adjustment (VECM(M)-BMSM)                   
    %                      Maximum likelihood estimation                       
    %                                  v1.0                                      
    %                   Copyright © 2013 Waleem Alausa                
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    % Vector Error Correction Bivariate Markov Switching Multifractal Model with Time Varying Speed 
    % of Long-Run Adjustment (VECM(M)-MSM), assuming that rho_m (the correlation between M^alpha and M^beta) is 1.
    % VECM(M) Bivariate MSM:
    %                                         R_s = b_s + a_s((Ms_1 -1 + Ms_2 -1 +. ...+ Ms_kbar-1))*(S-F) + e_s
    %                                         R_f = b_f + a_f((Mf_1 -1 + Mf_2 -1 +. ...+ Mf_kbar-1)-1)*(S-F) + e_f
    %                                         e_s = sigma_s*(Ms_1*Ms_2*. ...*Ms_kbar)
    %                                         e_f = sigma_f*(Mf_1*Mf_2*. ...*Mf_kbar)
    %
    % The estimation is done in two stages: First stage estimates [m01,m02,sigma1,sigma2,gamma_k,b, bs, bf, as, af]
    %                                       Second stage estimates [rho_e, lambda]
    %
    % INPUTS:
    %    R                   -   T-by-2 matrix of returns [log_spot_return log_futures_return]
    %    S                   -   T-by-1 vector of log_spot prices
    %    F                   -   T-by-1 vector of log_futures prices
    %    KBAR                -   The number of volatility components.
    %    n                   -   The number of trading periods in 1 year. E.g
    %                            n= 252 for the number of business days in a year. n=12 for monthly data.
    %    STARTINGVALS        -   [OPTIONAL] Starting values for optimization
    %                            [m01, m02, sigma1,sigma2, gamma_k, b, bs, bf, as, af, rho_e, lambda]
    %                               b       - (1,inf) 
    %                               m01     - (1,2]
    %                               m02     - (1,2]
    %                               gamma_k - (0,1) 
    %                               sigma1  - [0,inf)
    %                               sigma2  - [0,inf)
    %                               rho_e   - [-1,1]
    %                               lambda  - [0,1]
    %    OPTIONS             -   {OPTIONAL] User provided options structure
    %
    % OUTPUTS:
    %    output.parameters   -   1-by-12 vector of parameters 
    %                              [m01, m02, sigma1,sigma2, gamma_k, b, bs, bf, as, af,  rho_e, lambda] 
    %    output.se           -   1-by-12 row vector of parameter standard errors
    %    out.ect             -   T-by-2 matrix of residuals after removing the VECM terms.
    %    out.as              -   T-by-1 vector of time-varying speed of adjustment to long-run equilibrium 
    %                            for the spot equation.
    %    out.af              -   T-by-1 vector of time-varying speed of adjustment to long-run equilibrium
    %                            for the futures equation.
    %    out.ect             -   T-by-2 matrix of residuals after removing the VECM terms.
    %    L.LL1               -   The first-stage log-likelihood at the optimum(scalar)
    %    L.LL2               -   The second-stage log-likelihood at the optimum(scalar)
    %    L.LLs1              -   Individual first-stage daily log-likelihoods at optimum for data^alpha (T-by-1)
    %    L.LLs2              -   Individual first-stage daily log-likelihoods at optimum for data^beta (T-by-1)
    %    states.pi1          -   T-by-4^kbar state probability for data_alpha
    %    states.pi2          -   T-by-4^kbar state probability for data_beta 
    %    states.A1           -   2^kbar-by-2^kbar transition matrix for stage-1
    %    diagnostics         -   Structure of optimization output information.Useful for checking convergence problems
    %
    % REFERENCES:
    %    [1] Calvet, L., Adlai Fisher (2004). "How to Forecast long-run 
    %        volatility: regime-switching and the estimation of multifractal 
    %        processes". Journal of Financial Econometrics 2: 49–83.
    %    [2] Calvet, L., Adlai Fisher (2008). "Multifractal Volatility: Theory, 
    %        Forecasting and Pricing". Elsevier - Academic Press..
    %    [3] Calvet et al (2006). Volatility comovement: a multifrequency approach. 
    %        Journal of Econometrics, 131(1-2):179{215.
    % ------------------------------------------------------------------------- 
    switch nargin
        case 5
            [r, S, F, kbar, startingvals, LB1, UB1, LB2, UB2,n, options] = ...
                            VECM_BMSM_parameter_check(r,S,F,kbar,n);
        case 6
            [r, S, F, kbar, startingvals, LB1, UB1, LB2, UB2,n, options] = ...
                            VECM_BMSM_parameter_check(r,S,F,kbar,n,startingvals);
        case 7
            [r, S, F, kbar, startingvals, LB1, UB1, LB2, UB2,n, options] = ...
                            VECM_BMSM_parameter_check(r,S,F,kbar,n, startingvals, options);
        otherwise
            error('Number of inputs must be between 5 and 7');
    end


    % First stage maximum likelihood estimation for [m01,m02,sigma1,sigma2,gamma_k,b]
    % Set up template for transition matrix to save initializing at each iteration
    
    A_template = T_mat_template(kbar);
    
    % Get starting values (if none are supplied by user)
    startingvals = VECM_BMSM_starting_values(r, S, F, startingvals, kbar, A_template, n);
    start1=startingvals(1:10);
    start2=startingvals(11:12);
    
    [para1,LL1,exitflag1,output1]=fmincon('VECM_M_BMSM_likelihood1',start1,...
        [],[],[],[],LB1,UB1,[], options, kbar, r, S, F, A_template, n);
    
    % Remove the conditional mean (VECM) term and get residuals for estimating rho_e and lambda
    [ect_resid,as_m, af_m,mu_s, mu_f] = VECM_M_BMSM_ect(para1,r,S,F,kbar,n);
        
    [para2,LL2,exitflag2,output2]=fmincon('RES_BMSM_2stage_likelihood2',start2,...
        [],[],[],[],LB2,UB2,[],options ,kbar, ect_resid, para1(1:6),n);
    


    %Store parameter estimates and standard errors
    if isrow(para1) ==0
        para1=para1';
    end
    if isrow(para2) ==0
        para2=para2';
    end
    
    %Estimate standard errors using GMM approach
   se = VECM_M_BMSM_std_err(r, S, F, [para1 para2],kbar,n);

    out.para    =[para1 para2];
    out.para(3) =out.para(3)/sqrt(n);
    out.para(4) =out.para(4)/sqrt(n);
    out.se      =se;
    out.se(3)   =out.se(3)/sqrt(n);
    out.se(4)   =out.se(4)/sqrt(n); 
    out.ect     =ect_resid;
    out.as      =as_m;
    out.af      =af_m;
    out.mu_s    =mu_s;
    out.mu_f    =mu_f;

    L.LL1       =-LL1;
    L.LL2       =-LL2;
    
    % Compute log-likelihood and daily likelihood values
    if nargout>2
       [~, LLs1, LLs2, pi_mat1, pi_mat2, A1] = VECM_M_BMSM_likelihood1(para1, kbar, r, S, F, A_template, n);
       
       L.LLs1=LLs1;
       L.LLs2=LLs2;
       states.pi1=pi_mat1;
       states.pi2=pi_mat2;
       states.A1=A1;

    end
    
   
    % Set up diagnostics for output
    diagnostics.EXITFLAG1=exitflag1;
    diagnostics.EXITFLAG2=exitflag2;
    diagnostics.ITERATIONS1=output1.iterations;
    diagnostics.ITERATIONS2=output2.iterations;
    diagnostics.FUNCCOUNT1=output1.funcCount;
    diagnostics.FUNCCOUNT2=output2.funcCount;
    diagnostics.MESSAGE1=output1.message;
    diagnostics.MESSAGE2=output2.message;
    out.diagnostics=diagnostics;

    % Save a template for the transition matrix to prevent initialization at
    % each log-likelihood iteration
    function A = T_mat_template(kbar)  
        A = zeros((2^kbar),(2^kbar));
        for i =0:2^kbar-1       
            for j = i:(2^kbar-1)-i  
                A(i+1,j+1) = bitxor(i,j);
            end
        end
    end

end