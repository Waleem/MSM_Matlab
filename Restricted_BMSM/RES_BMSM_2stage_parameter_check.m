function [data, kbar, startvals,LB1,UB1,LB2,UB2,n,options] = RES_BMSM_2stage_parameter_check(data, kbar,n,startingvals, options)
%Parameter check for 2-stage restricted  bivariate MSM
% -------------------------------------------------------------------------
%     data    
% -------------------------------------------------------------------------
if size(data,2) > 2
    
    data = data';
    
end

if size(data,2) > 2 
    error('data series must be a T-by-2 matrix')
elseif isempty(data)
    error('data is empty')
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

if nargin>3
    
    if (startingvals(1) < 1) || (startingvals(1) > 1.9999)
        error('m01 must be between (1,1.99]')
    end
    if (startingvals(2) < 1) || (startingvals(2) > 1.9999)
        error('m02 must be between (1,1.99]')
    end
    
    if (startingvals(3) <= 0)
        error('sigma1 must be a positive (non-zero) number')
    end
    if (startingvals(4) <= 0)
        error('sigma2 must be a positive (non-zero) number')
    end
    
    if (startingvals(5) < 0.00001) || (startingvals(5) > 0.9999)
        error('gamma_k must be between [0,1]')
    end
    
    if (startingvals(6) < 1)
    error('b must be greater than 1')
    end
    
    if (startingvals(7) < -0.9999) || (startingvals(7) > 0.9999)
        error('rho_e must be between [-1,1]')
    end
    
    if (startingvals(8) < 0.00001) || (startingvals(8) > 0.9999)
        error('lambda must be between [0,1]')
    end
   
    startvals = startingvals;
   
else
    startvals = [];
   
end
% -------------------------------------------------------------------------
%     Upper and lower bounds    
% -------------------------------------------------------------------------
LB1 = [1,     1,      0.0001, 0.0001, 0.001,  1];
UB1 = [1.9999,1.9999, 100,    100,    0.9999, 100];
LB2 = [-0.9999, 0.0001]; 
UB2 = [0.9999,  0.9999];




% -------------------------------------------------------------------------
%     options    
% -------------------------------------------------------------------------
if nargin>4 && ~isempty(options)
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
    options  =  optimset(options , 'Algorithm'   ,'sqp');     % Note that sqp is only available for later versions of matlab.
    %                                                                         % If you got the error that sqp is not an option, that means you have 
    %                                                                         % an older version of matlab. So change this to either 'Interior-point' or 'Active-set'.
    options  =  optimset(options , 'UseParallel' , 'Always'); % If you have a multicore system, this may speed up estimation.
 
end