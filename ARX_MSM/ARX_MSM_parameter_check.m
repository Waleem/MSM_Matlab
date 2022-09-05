function [Y, xmat, kbar, startingvals, nw, Lb, Ub, A, b, options] =...
    ARX_MSM_parameter_check(y, p, kbar, model, X, n, nw, startingvals, options)


if size(y,2) > 1
    y = y';
end

if size(y,1) < 1 || size(y,2) >1
    error('Dependent variable must be a vector')
elseif isempty(y)
    error('Dependent variable is empty')
end

if isempty(X)==0
   if size(X,1) ~= size(y,1)
       error('Depndent variable and factor mayrix must have the same number of rows')
   end
end


if (length(kbar) > 1) || kbar < 1 
    error('k must be a positive integer')
end
kbar = floor(kbar);

if p < 1 
    error('AR order should be a positive integer or a vector of positive integer orders')
end
p=floor(p);


if ismember(model, [1 2 3]) ==0
    error('Model must be one of 1 = ARX-MSM, 2 = ARX(M)-MSM, 3 = ARX-MSM-M')
end

if n < 1 
    error('n must be a positive integer. E.g n = 252 for daily data, n=12 for monthly data')
end
n=floor(n);

if isempty(nw)==0
    if nw < 1 
        error('nw (Newey-West lag order) must be a positive scalar.')
    end
else
    nw=0;
end


%-------------------------------------------------------------------------
%     starting values    
% -------------------------------------------------------------------------
if nargin>7
    if (startingvals(1) < 1)
    error('b must be greater than 1')
    end
    if (startingvals(2) < 1) || (startingvals(2) > 1.9999)
        error('m0 must be between (1,1.99]')
    end
    if (startingvals(3) < 0.00001) || (startingvals(3) > 0.9999)
        error('gamma_k must be between [0,1]')
    end
    if (startingvals(4) < 0.00001)
        error('sigma must be a positive (non-zero) number')
    end
    %if (in(4) > 1)
        %warning('Sigma value is very large - consider using smaller value')
    %end
    if length(startingvals) ~= 4+1+numel(p)+size(X,2) && model~=3
       error('The number of initial parameters is incorrect')
    elseif length(startingvals) ~= 4+2+numel(p)+size(X,2) && model==3
        error('The number of initial parameters is incorrect')
    end
else
    startingvals = [];
end

% -------------------------------------------------------------------------
%     Upper and lower bounds 
%     Remember abs(bi)<1 for covariance stationarity in the AR process
% -------------------------------------------------------------------------

%         [b m0 gamma_kbar sigma]
Lb_MSM = [1, 1, 0.0001, 0.0001];
Ub_MSM = [50, 1.999, 0.999, 50]; % Remember the parameter restrictions imposed by MSM, should you decide to change the MSM bounds.
Lb_AR  = -0.999*ones(1, numel(p));
Ub_AR  =  0.999*ones(1, numel(p));
Lb_X   = -inf*ones(1, size(X,2));
Ub_X   =  inf*ones(1, size(X,2));

if model ==3
    Lb     = [Lb_MSM -inf Lb_AR Lb_X -inf];
    Ub     = [Ub_MSM  inf Ub_AR Ub_X  inf];
else
    Lb     = [Lb_MSM -inf Lb_AR Lb_X];
    Ub     = [Ub_MSM  inf Ub_AR Ub_X];
end

% -------------------------------------------------------------------------
%     Linear equality constraints 
%     AR(p) => b1+b2+...+bp <1 => Covariance stationaity
% -------------------------------------------------------------------------
if model ==3
    A = [zeros(1,5) ones(1, numel(p)) zeros(1, size(X,2)) 0];
    b = 1;
else
    A = [zeros(1,5) ones(1, numel(p)) zeros(1, size(X,2))];
    b = 1;
end

% -------------------------------------------------------------------------
%     Construct dependent variable vector and independent variable matrix
% -------------------------------------------------------------------------
 
Y = y(max(p)+1 : end);
if isempty(X)==0
   X = X(max(p)+1 : end,:);
end
ylag=zeros(size(y,1)-max(p),1);

for i=1:numel(p)
    ylag(:,i) = y(p(end)-p(i)+1: end-p(i)); 
end

xmat = [ones(size(Y,1),1) ylag X];


% -------------------------------------------------------------------------
%     options    
% -------------------------------------------------------------------------
if nargin>8 && ~isempty(options)
    try
        optimset(options);
    catch
        error('options is not a valid minimization option structure');
    end
else
    %Setup the options in case none is provided
    options  =  optimset('fmincon');
    options  =  optimset(options , 'TolFun'      , 1e-005);
    options  =  optimset(options , 'TolX'        , 1e-005);
    options  =  optimset(options , 'Display'     , 'iter');
    options  =  optimset(options , 'Diagnostics' , 'on');
    options  =  optimset(options , 'LargeScale'  , 'off');
    options  =  optimset(options , 'MaxFunEvals' , 2000);
    options  =  optimset(options , 'MaxSQPIter'  , 2000);
    options  =  optimset(options , 'Algorithm'   ,'sqp');     % Note that sqp is only available for later versions of matlab.
    %                                                         % If you got the error that sqp is not an option, that means you have 
    %                                                         % an older version of matlab. So change this to either 'Interior-point' or 'Active-set'.
    options  =  optimset(options , 'UseParallel' , 'Always'); % If you have a multicore system, this may speed up estimation.
end