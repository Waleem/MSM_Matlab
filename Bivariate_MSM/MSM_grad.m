% =========================================================================
% numerical derivative function
% calculates the numerical derivative of function f at the parameter value
% bd (uses central difference)
% inputs: f, needs to be a function handle, ie if the function's name is
% test you need to pass @test 
% f can be (columns) vector valued, ie its output
% can be (q x 1)
% bd, this is the parameter value at which to evaluate the derivative
% this vector has dim (k x 1)
% y and x, further inputs into function f
% eg function a = test(b,y,x)
% outputs: der, a (q x k) vector of numerical derivatives
% =========================================================================

function der = MSM_grad(f,bd,varargin)

k = size(bd,1);
temp = repmat(feval(f,bd,varargin{:}),k,2);

% compute stepsize (as in GAUSS procedure gradp)
abd = abs(bd);
if bd ~= 0
    dabd = bd./abd;
else
    dabd = 1;
end

h1 = [abd (1e-2*ones(k,1))];
h = 1e-8 * max(h1,[],2) .* dabd;
temp = bd + h;
h = temp - bd;

for i = 1:k
    b_temp = repmat(bd,1,2); 
    b_temp(i,1) = b_temp(i,1) + h(i,1);
    b_temp(i,2) = b_temp(i,2) - h(i,1);
    temp1(:,i) = feval(f,b_temp(:,1),varargin{:});
    temp2(:,i) = feval(f,b_temp(:,2),varargin{:});
end

q = size(temp1,1);
der = (temp1 - temp2)./(2*repmat(h,1,q)');

end





