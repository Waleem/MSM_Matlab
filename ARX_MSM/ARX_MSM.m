function [para, LL, se, LLs, states, summary, diagnostics] = ARX_MSM(y, p, kbar, model, X, n, nw, startingvals, options)
% -------------------------------------------------------------------------
%                   Autoregressive Markov Switching Multifractal (ARX-MSM) Model                 
%                               Maximum likelihood estimation                       
%                  
% -------------------------------------------------------------------------
% Author:
% Waleem Alausa, alausa.babs@gmail.com
% Date:     19/Feb/2014: Initial development
%
%--------------------------------------------------------------------------
% USAGE:
%    [para, LL, se, LLs, states, summary, diagnostics] = ARX_MSM(y, p, kbar, model, X, n)
%    [para, LL, se, LLs, states, summary, diagnostics] = ARX_MSM(y, p, kbar, model, X, n, nw)
%    [para, LL, se, LLs, states, summary, diagnostics] = ARX_MSM(y, p, kbar, model, X, n, nw, startingvals)
%    [para, LL, se, LLs, states, summary, diagnostics] = ARX_MSM(y, p, kbar, model, X, n, nw, startingvals, options)
%
% INPUTS:
%    Y                   -   A (T x 1) column (or row) of price or returns
%    P                   -   A positive scalar or vector of autoregressive order(s)
%    KBAR                -   A positive scalar for number of frequency/volatility components to be included in MSM(kbar)
%    Model               -   A scalar for model type: 1 = ARX-MSM, 2 = ARX(M)-MSM, 3 = ARX-MSM-M
%    X                   -   A (T x N) vector of factors (regressors) for
%                            the mean process. But it can also be empty forthe case of AR-MSM
%    n                   -   A scalar for the number of trading periods in 1 year. E.g
%                            n= 252 for the number of business days in a year. n=12 for monthly data.
%                            This is needed to estimate unconditional volatility in MSM
%    NW                  -   [OPTIONAL] A positive scalar for the Newey-West lag order to use in robust standard errors
%    STARTINGVALS        -   [OPTIONAL] Starting values for optimization
%                            [b  m0  gamma_k  sigma b0 b1 b2 ... b_AR a_x1 a_x2 ... a_xN  x_sigma(if model=3)]
%                               b       - (1,inf) 
%                               m0      - (1,2] 
%                               gamma_k - [0,1] 
%                               sigma   - [0,inf)
%                             Other parameters can be any value between +inf and -inf. But choose wisely.
%                             If you don't know what values to use, better to leave STARTINGVALS empty.
%    OPTIONS             -   [OPTIONAL] User provided options structure
%
% OUTPUTS:
%    Para                -   A 4+numel(p)+size(X,2)+1-by-1 column vector of parameters 
%                            [b  m0  gamma_k  sigma b0 b1 b2 ... b_AR a_x1 a_x2 ... a_xN ], for model=1 or 2, or
%                            A 4+numel(p)+size(X,2)+2-by-1 column vector of parameters 
%                            [b  m0  gamma_k  sigma b0 b1 b2 ... b_AR a_x1 a_x2 ... a_xN a_sigma], for model=3
%    LL                  -   The log-likelihood at the optimum
%    se                  -   A 4+numel(p)+1-by-1 column vector of parameter standard errors
%    LLs                 -   Individual daily log-likelihoods at optimum
%    states.p_mat        -   T-by-2^kbar state probability matrix 
%    states.A            -   2^kbar-by-2^kbar transition matrix
%    states.M            -   kbar-by-2^kbar matrix of volatility components
%    states.g_m          -   1-by-2^kbar vector of state values
%    Summary             -   Summary of results, including r_sq, adjusted r_sq, AIC, BIC, residuals(et) and yhat
%    diagnostics         -   Structure of optimization output information.Useful for checking convergence problems
%
% ASSOCIATED FILES:
%    ARX_MSM_parameter_check.m, ARX_MSM_starting_values.m, ARX_MSM_likelihood.m, ARX_MSM_std_err.m 
%
% REFERENCES:
%    [1] Calvet, L., Adlai Fisher (2004). "How to Forecast long-run 
%        volatility: regime-switching and the estimation of multifractal 
%        processes". Journal of Financial Econometrics 2: 49–83.
%    [2] Calvet, L., Adlai Fisher (2008). "Multifractal Volatility: Theory, 
%        Forecasting and Pricing". Elsevier - Academic Press.


switch nargin
    case 6
         [Y, xmat, kbar, startingvals, nw, Lb, Ub, A, b, options] =...
            ARX_MSM_parameter_check(y, p, kbar, model, X, n,[]);
    case 7
         [Y, xmat, kbar, startingvals, nw, Lb, Ub, A, b, options] =...
            ARX_MSM_parameter_check(y, p, kbar, model, X, n,nw);
    case 8
        [Y, xmat, kbar, startingvals, nw, Lb, Ub, A, b, options] =...
            ARX_MSM_parameter_check(y, p, kbar, model, X, n, nw, startingvals);
    case 9
        [Y, xmat, kbar, startingvals, nw, Lb, Ub, A, b, options] =...
            ARX_MSM_parameter_check(y, p, kbar, model, X, n,nw, startingvals, options);
    otherwise
        error('Number of inputs must be between 6 and 9');
end


% Get starting values if none are supplied by user:
if isempty(startingvals)
   disp('No starting values entered: Using grid-search')
   startingvals = ARX_MSM_starting_values(Y, xmat, kbar, n, model);
end

% Start seraching for minimum negative likelihood value
[parameters, LL, exitflag, diagnosis,~,~,hessian]=fmincon('ARX_MSM_likelihood', startingvals, A, b, [],[],Lb,Ub,[],...
                                           options, kbar, Y, xmat, n, model, p);



%Store parameter estimates and standard errors
para=parameters;
if size(para,2) >1
    para = para';
end
para(4)=  para(4)/sqrt(n);
LL     = -LL;

if nargout>2
    %Estimate standard errors using GMM approach
    if nw~=0
        disp('Estimating Robust Standard Errors')
    else
        disp('Estimating Standard Errors')
    end
    se = ARX_MSM_std_err(parameters, kbar, Y, xmat, n, model, p, nw);
end

if nargout > 3
    [~, LLs, A, p_mat, M, g_m, et, yhat] = ARX_MSM_likelihood(parameters, kbar, Y, xmat, n, model, p);
    states.p_mat = p_mat;
    states.A     = A;
    states.M     = M;
    states.g_m   = g_m; 

    T=size(Y,1);
    summary.Rsquared = (1-sum(et.^2)/sum((Y-mean(Y)).^2));
    summary.AdjRsquared = 1-((T-1)/(T-length(parameters))*(1-summary.Rsquared));
    summary.LLF = LL;
    summary.AIC = -2*LL+2*length(parameters);
    summary.BIC = -2*LL+length(parameters)*log(size(Y,1)); 
    %summary.et  = et;
    %summary.yhat= yhat;
    summary.y   = Y;
    summary.x   =xmat;
    %summary.LLs =LLs;
    %summary.se2 =sqrt(diag(inv(hessian)));

    % Set up diagnostics for output
    diagnostics.EXITFLAG   = exitflag;
    diagnostics.ITERATIONS = diagnosis.iterations;
    diagnostics.FUNCCOUNT  = diagnosis.funcCount;
    diagnostics.MESSAGE    = diagnosis.message;


end

end

