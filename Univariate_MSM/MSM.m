function [output,LL,LLs,states,diagnostics] = MSM(data, kbar,n, startingvals, options)
% -------------------------------------------------------------------------
%                   Markov Switching Multifractal (MSM)                   
%                      Maximum likelihood estimation                       
%                                  v1.0                                      
%                Copyright © 2010 Multifractal-finance.com                 
% -------------------------------------------------------------------------
%                    Original code modified by Waleem Alausa
%                        This version: 31/Aug/2014
%--------------------------------------------------------------------------
% USAGE
%    [output] = MSM(DATA, K, n)
%    [output, LL, LLs,states, DIAGNOSTICS] = MSM(DATA, K, n, STARTING_VALUES)
%    [output, LL, LLs,states, DIAGNOSTICS] = MSM(DATA, K, n, STARTING_VALUES, OPTIONS)
%
% INPUTS:
%    DATA                -   A column (or row) of mean zero data
%    KBAR                -   The number of frequency components
%    n                   -   A scalar for the number of trading periods in 1 year. E.g
%                            n= 252 for the number of business days in a year. n=12 for monthly data.
%    STARTINGVALS        -   [OPTIONAL] Starting values for optimization
%                            [b, m0, gamma_k, sigma]
%                               b       - (1,inf) 
%                               m0      - (1,2] 
%                               gamma_k - (0,1) 
%                               sigma   - [0,inf)
%    OPTIONS             -   {OPTIONAL] User provided options structure
%
% OUTPUTS:
%    output.parameters   -   A 4x1 column vector of parameters 
%                            [b, m0, gamma_k, sigma]
%    output.se           -   A 4x1 column vector of parameter standard errors
%    LL                  -   The log-likelihood at the optimum
%    LLs                 -   Individual daily log-likelihoods at optimum
%    states              -   Structure of states quantities:A=transition
%                            matrix, state_vals= possible state values,
%                            state_prob=probability of states.
%    diagnostics         -   Structure of optimization output information.
%                            Useful for checking convergence problems
%
% ASSOCIATED FILES:
%    MSM_likelihood.m, MSM_parameter_check.m, MSM_starting_values.m
%
% REFERENCES:
%    [1] Calvet, L., Adlai Fisher (2004). "How to Forecast long-run 
%        volatility: regime-switching and the estimation of multifractal 
%        processes". Journal of Financial Econometrics 2: 49–83.
%    [2] Calvet, L., Adlai Fisher (2008). "Multifractal Volatility: Theory, 
%        Forecasting and Pricing". Elsevier - Academic Press..
% ------------------------------------------------------------------------- 
switch nargin
    case 3
        [data, kbar, startingvals, LB, UB, options] = ...
            MSM_parameter_check(data, kbar, n);
    case 4
        [data, kbar, startingvals, LB, UB, options] = ...
            MSM_parameter_check(data, kbar, n, startingvals);
    case 5
        [data, kbar, startingvals, LB, UB, options] = ...
            MSM_parameter_check(data, kbar, n, startingvals, options);
    otherwise
        error('Number of inputs must be between 3 and 5');
end

% Set up template for transition matrix to save initializing each iteration
A_template = T_mat_template(kbar);

% Get starting values (if none are supplied by user)
[startingvals,LLs] = MSM_starting_values(data,startingvals,kbar,n,A_template);

% Start seraching for minimum negative likelihood value
[para,LL,exitflag,diagnosis]=fmincon('MSM_likelihood',startingvals,[],[],[],[],LB,UB,[],options,kbar,n, data,A_template);

%Store parameter estimates and standard errors
if size(para,2)>1
    para=para';
end
output.para=para;
output.para(4)=output.para(4)/sqrt(n);


output.se=MSM_std_err(para,kbar, n, data, A_template);
output.se(4)=output.se(4)/sqrt(n);
LL=-LL;

% Compute log-likelihood and daily likelihood values
if nargout>2
    [LL,LLs,A,g_m,pi_mat] = MSM_likelihood(para,kbar,n ,data,A_template);
    LL = -LL;
    states.A=A;
    states.state_vals=g_m;
    states.state_prob=pi_mat;
end

% Set up diagnostics for output
diagnostics.EXITFLAG=exitflag;
diagnostics.ITERATIONS=diagnosis.iterations;
diagnostics.FUNCCOUNT=diagnosis.funcCount;
diagnostics.MESSAGE=diagnosis.message;

% Save a template for the transition matrix
function A = T_mat_template(kbar)  
    A = zeros((2^kbar),(2^kbar));
    for i =0:2^kbar-1       
        for j = i:(2^kbar-1)-i  
            A(i+1,j+1) = bitxor(i,j);
        end
    end
