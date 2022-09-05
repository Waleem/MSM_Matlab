function [data, S, F, kbar, startvals,LB1,UB1,LB2,UB2,n,options] = ...
                                 VECM_BMSM_parameter_check(data,S,F,kbar,n, startingvals, options)
%Parameter check for Vector Error Correction Bivariate MSM
% -------------------------------------------------------------------------
%     data    
% -------------------------------------------------------------------------
if size(data,2) > 2
    
    data = data';
    
end

if size(data,2) > 2 
    error('data series must be a T-by-2 matrix')
elseif isempty(data(:,:))
    error('data is empty')
end

if size(S,2) > 1
    
    S = S';
    
end

if size(S,2) > 1 
    error('S must be a T-by-1 vector')
elseif isempty(S(:))
    error('S is empty')
end

if size(F,2) > 1
    
    F = F';
    
end

if size(F,2) > 1 
    error('F must be a T-by-1 vector')
elseif isempty(F(:))
    error('F is empty')
end

%-------------------------------------------------------------------------
%     Check starting values    
% -------------------------------------------------------------------------
if (kbar < 1) 
    error('k must be a positive integer')
end
kbar = floor(kbar);
if (n < 1) 
    error('n must be a positive integer')
end
n = floor(n);

if nargin>5
    
    if (startingvals(1) < 1) || (startingvals(1) > 1.99)
        error('m01 must be between (1,1.99]')
    end
    if (startingvals(2) < 1) || (startingvals(2) > 1.99)
        error('m02 must be between (1,1.99]')
    end
    
    if (startingvals(3) <= 0)
        error('sigma1 must be a positive (non-zero) number')
    end
    if (startingvals(4) <= 0)
        error('sigma2 must be a positive (non-zero) number')
    end
    
    if (startingvals(5) < 0.00001) || (startingvals(5) > 0.99999)
        error('gamma_k must be between [0,1]')
    end
    
    if (startingvals(6) < 1)
    error('b must be greater than 1')
    end
    
    if (startingvals(11) < -0.99999) || (startingvals(11) > 0.99999)
        error('rho_e must be between [-1,1]')
    end
    
    if (startingvals(12) < 0.00001) || (startingvals(12) > 0.99999)
        error('lambda must be between [0,1]')
    end
   
    startvals = startingvals;
   
else
    startvals = [];
   
end
% -------------------------------------------------------------------------
%     Upper and lower bounds    
% -------------------------------------------------------------------------
%     [m01,   m02,     sigma1, sigma2, gamma_k, b,   bs,    bf,   as,   af]  
LB1 = [1,      1,      0.0001, 0.0001, 0.0001,   2,  -inf,  -inf, -10, -10];
UB1 = [1.9999, 1.9999, 100,    100,    0.9999,  500, inf,   inf,  10,  10];
%     [rho_e,    lambda] 
LB2 = [-0.9999, 0.0001];
UB2 = [ 0.9999, 0.9999];

% -------------------------------------------------------------------------
%     options    
% -------------------------------------------------------------------------
if nargin>6 && ~isempty(options)
    try
        optimset(options);
    catch
        error('options is not a valid minimization option structure');
    end
else
    %Setup the options in case of none provided
    options  =  optimset('fmincon');
    options  =  optimset(options , 'TolFun'      , 1e-005);
    options  =  optimset(options , 'TolX'        , 1e-005);
    options  =  optimset(options , 'Display'     , 'iter');
    options  =  optimset(options , 'Diagnostics' , 'on');
    options  =  optimset(options , 'LargeScale'  , 'on');
    options  =  optimset(options , 'MaxFunEvals' , 1000);
    options  =  optimset(options , 'MaxSQPIter'  , 1000);
    options  =  optimset(options , 'Algorithm'   ,'sqp');
    %options  = optimset(options , 'FinDiffType' ,'central');
    options  = optimset(options , 'UseParallel' , 'Always');
 
end