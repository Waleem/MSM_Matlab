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