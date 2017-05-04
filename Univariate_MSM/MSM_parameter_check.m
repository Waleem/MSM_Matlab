function [data, kbar, startingvals, LB, UB, options] =...
    MSM_parameter_check(data, kbar,n, in, options)

% -------------------------------------------------------------------------
%     data    
% -------------------------------------------------------------------------
if size(data,2) > 1
    data = data';
end
if size(data,2) > 1 || length(data) == 1
    error('data series must be a vector')
elseif isempty(data)
    error('data is empty')
end

%-------------------------------------------------------------------------
%     starting values    
% -------------------------------------------------------------------------
if (kbar < 1) 
    error('k must be a positive integer')
end
kbar = floor(kbar);

if nargin>3
    if (in(1) < 1)
    error('b must be greater than 1')
    end
    if (in(2) < 1) || (in(2) > 1.99)
        error('m0 must be between (1,1.99]')
    end
    if (in(3) < 0.00001) || (in(3) > 0.99999)
        error('gamma_k be between [0,1]')
    end
    if (in(4) < 0.00001)
        error('sigma must be a positive (non-zero) number')
    end
    %if (in(4) > 1)
        %warning('Sigma value is very large - consider using smaller value')
    %end
    startingvals = in;
else
    startingvals = [];
end
% -------------------------------------------------------------------------
%     Upper and lower bounds    
% -------------------------------------------------------------------------

%LB = [1, 1, 0.001, 0.0001];
%UB = [50, 1.99, 0.99999, 50];
LB = [1, 1, 0.0001, 0.0001];
UB = [50, 1.9999, 0.9999, 50];

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
    options  =  optimset(options , 'LargeScale'  , 'off');
    options  =  optimset(options , 'MaxFunEvals' , 2000);
    options  =  optimset(options , 'MaxSQPIter' , 2000);
    options  =  optimset(options , 'Algorithm'   ,'sqp'); % Note that sqp is only available for later versions of matlab.
    %                                                         % If you got the error that sqp is not an option, that means you have 
    %                                                         % an older version of matlab. So change this to either 'Interior-point' or 'Active-set'.
    options  =  optimset(options , 'UseParallel' , 'Always'); % If you have a multicore system, this may speed up estimation.

end