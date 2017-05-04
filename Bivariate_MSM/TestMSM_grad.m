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


load caret
output = MSM(caret, 2, 252)
para = output.para
para(4)=para(4)*sqrt(252)

kbar = 2; n=252;
    A = zeros((2^kbar),(2^kbar));
    for i =0:2^kbar-1       
        for j = i:(2^kbar-1)-i  
            A(i+1,j+1) = bitxor(i,j);
        end
    end
    
bd=para;
k = size(bd,1);



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
    temp1(:,i) = MSM_LLs(b_temp(:,1),kbar,n,data,A);
    temp2(:,i) = MSM_LLs(b_temp(:,2),kbar,n,data,A);
end

q = size(temp1,1);
der = (temp1 - temp2)./(2*repmat(h,1,q)');







