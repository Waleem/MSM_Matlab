function [output,L,states,diagnostics] = RES_BMSM_2stage(data,kbar,n,startingvals, options)
    % -------------------------------------------------------------------------
    %                   Bivariate Markov Switching Multifractal (MSM)                   
    %                      Maximum likelihood estimation                       
% -------------------------------------------------------------------------
% Author:
% Waleem Alausa, alausa.babs@gmail.com
% Date:     19/Feb/2014: Initial development
%
%--------------------------------------------------------------------------
    %Estimates 2-stage bivariate MSM, with the assumption that rho_m (the correlation between M^alpha and M^beta) is 1. 
    %The estimation is done in two stages: First stage estimates [m01,m02,sigma1,sigma2,gamma_k,b]
    %                                      Second stage estimates [rho_e, lambda]
    % USAGE
    %    [output]                      = RES_BMSM_2stage(data,kbar,n)
    %    [output,L,states,diagnostics] = RES_BMSM_2stage(data,kbar,n,startingvalues)
    %    [output,L,states,diagnostics] = RES_BMSM_2stage(data,kbar,n,startingvalues,options)
    %
    % INPUTS:
    %    DATA                -   2 columns (or row) of mean zero data
    %    KBAR                -   The number of volatility components.
    %    n                   -   The number of trading periods in 1 year. E.g
    %                            n= 252 for the number of business days in a year. n=12 for monthly data.
    %    STARTINGVALS        -   [OPTIONAL] Starting values for optimization
    %                            [m01, m02, sigma1, sigma2, gamma_k, b, rho_e, lambda]
    %                               b        - (1,inf) 
    %                               m01      - [1,2]
    %                               m02      - [1,2]
    %                               gamma_k  -  [0,1] 
    %                               sigma1   - [0,inf)
    %                               sigma2   - [0,inf)
    %                               rho_e    - [-1,1]
    %                               lambda   - [0,1]
    %    OPTIONS             -   {OPTIONAL] User provided options structure
    %
    % OUTPUTS:
    %    output.parameters   -   A 1-by-8 row vector of parameters 
    %                            [m01, m02, sigma1, sigma2, gamma_kbar, b, rho_e, lambda]
    %    output.se           -   A 1-by-8 row vector of parameter standard errors
    %    L.LL1=LLstep1       -   The first-stage log-likelihood at the optimum(scalar)
    %    L.LL2=LLstep2       -   The second-stage log-likelihood at the optimum(scalar)
    %    L.LLs=LLs           -   Individual second-stage daily log-likelihoods at optimum (T-by-1)
    %    L.LLs1=LLs1         -   Individual first-stage daily log-likelihoods at optimum for data^alpha (T-by-1)
    %    L.LLs2=LLs2         -   Individual first-stage daily log-likelihoods at optimum for data^beta (T-by-1)
    %    states.pi1=pi_mat1  -   T-by-4^kbar state probability for data_alpha
    %    states.pi2=pi_mat2  -   T-by-4^kbar state probability for data_beta 
    %    states.A1=A1        -   2^kbar-by-2^kbar transition matrix for stage-1
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
        case 3
            [data, kbar, startingvals, LB1, UB1, LB2, UB2,n, options] = ...
                RES_BMSM_2stage_parameter_check(data, kbar,n);
        case 4
            [data, kbar, startingvals, LB1, UB1, LB2, UB2,n, options] = ...
               RES_BMSM_2stage_parameter_check(data, kbar,n,startingvals);
        case 5
            [data, kbar, startingvals, LB1, UB1, LB2, UB2,n, options] = ...
               RES_BMSM_2stage_parameter_check(data, kbar,n, startingvals, options);
        otherwise
            error('Number of inputs must be between 3 and 5');
    end


    % First stage maximum likelihood estimation for [m01,m02,sigma1,sigma2,gamma_k,b]
    
    
    A_template = T_mat_template(kbar);
    
    % Get starting values (if none are supplied by user)
    startingvals = RES_BMSM_2stage_starting_values(data,startingvals,kbar,A_template,n);
    start1=startingvals(1:6);
    start2=startingvals(7:8);
    
    [para1,LL1,exitflag1,output1]=fmincon('BMSM_2stage_likelihood1',start1,...
        [],[],[],[],LB1,UB1,[],options,kbar,data,A_template,n);

    [para2,LL2,exitflag2,output2]=fmincon('RES_BMSM_2stage_likelihood2',start2,...
        [],[],[],[],LB2,UB2,[],options,kbar,data,para1,n);
    


    %Store parameter estimates and standard errors
    if isrow(para1) ==0
        para1=para1';
    end
    if isrow(para2) ==0
        para2=para2';
    end
    
%     %Estimate standard errors using GMM approach
     se = RES_BMSM_2stage_std_err(data,[para1 para2],kbar,n);

    output.para    =[para1 para2];
    output.para(3) =out.para(3)/sqrt(n);
    output.para(4) =out.para(4)/sqrt(n);
    output.se      =se;
    output.se(3)   =out.se(3)/sqrt(n);
    output.se(4)   =out.se(4)/sqrt(n);    

    L.LL1       =LL1;
    L.LL2       =LL2;
    
    % Compute log-likelihood and daily likelihood values
    if nargout>2
       [~,LLs1, LLs2,pi_mat1,pi_mat2,A1] = BMSM_2stage_likelihood1(para1,kbar,data,A_template,n);

       L.LLs1=LLs1;
       L.LLs2=LLs2;

       states.pi1=pi_mat1;
       states.pi2=pi_mat2;
       states.A1=A1;
    
       [~,LLs] = RES_BMSM_2stage_likelihood2(para2,kbar,data,para1,n);
       L.LLs=LLs;

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