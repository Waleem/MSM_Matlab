function [yhat, para] = ARX_MSM_Fcast(y, p, k, model, X, n, nw, startingvals)

Xt1 = X(end,:);
p1  = p-1; 
yt  = y;
yt1 = [];

for i=1:numel(p1)
    yt1 = [yt1 y(end-p1(i))];
end

if isempty( startingvals)
    para = ARX_MSM(yt, p, k, model, X(1:end-1,:), n, nw);
else
    para = ARX_MSM(yt, p, k, model, X(1:end-1,:), n, nw,startingvals);
end

if size(para,2) >1
    para = para'; 
end

switch model
        case 1 
            
            yhat  = [1 yt1 Xt1]*para(5:end);
            
        case 2 
            
            para2    = para;
            para2(4) = para2(4)*sqrt(n);
            [Y, xmat] =   getxmat(yt, p, X(1:end-1,:));
            [~, ~, A, p_mat, M] = ARX_MSM_likelihood(para2, k, Y, xmat, n, model, p);
            phat = p_mat(end,:)*A;
            pm   = MSM_marginals(phat,para(2),M,k);
            em   = para(2)*pm + (2-para(2))*(1-pm);
            bm   = sum(em - 1,2)*para(6:6+numel(p)-1);
            yhat  = [1 yt1 Xt1]*[para(5); bm; para(numel(p)+6:end)];

        case 3
            
            para2    = para;
            para2(4) = para2(4)*sqrt(n);
            [Y, xmat] =   getxmat(yt, p, X(1:end-1,:));
            [~, ~, A, p_mat, ~, g_m] = ARX_MSM_likelihood(para2, k, Y, xmat, n, model, p);
            phat = p_mat(end,:)*A;
            ht   = para(4)*(phat*g_m');
            yhat = [1 yt1 Xt1 ht]*para(5:end);
end

end


% -------------------------------------------------------------------------
%     Construct dependent variable vector and independent variable matrix
% -------------------------------------------------------------------------
function [Y, xmat] =   getxmat(y, p, X) 
Y = y(max(p)+1 : end);
if isempty(X)==0
   X = X(max(p)+1 : end,:);
end
ylag=zeros(size(y,1)-max(p),1);

for i=1:numel(p)
    ylag(:,i) = y(p(end)-p(i)+1: end-p(i)); 
end

xmat = [ones(size(Y,1),1) ylag X];

end
