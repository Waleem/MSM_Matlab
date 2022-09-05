function [se,VCV,A,B,scores,hess,gross_scores] = ARX_MSM_std_err(para, kbar, Y, xmat, n, model, p, nw)
 
% Compute robust standard errors using variance covariance matrix computed 
% numerically, including Newey-West style score covariance using 2-sided
% derivatives.

% This code benefits heavily from the robustvcv.m code of the MFEToolbox by 
% Kevin Sheppard (% kevin.sheppard@economics.ox.ac.uk)
           
   % T=size(Y,1);
    if size(para,2) >1
         para = para';
    end
    if kbar==1
       gt  = MSM_grad('get_ARX_MSM_LLs',para,kbar, Y, xmat, n, model, p);
       gt = gt(:,2:end); %Exclude b from the estimation of covariance matrix because b is not defined at k=1.
                           %Else the standard errors for all the parameters will be inaccurate or may even return Inf
        J1    = gt'*gt;
        se    = [nan; sqrt(diag(inv(J1)))];
        se(4) = se(4)/sqrt(n);
        
    else
        
        [VCV,A,B,scores,hess,gross_scores]=robustvcv('get_ARX_MSM_likelihood',para,nw,kbar, Y, xmat, n, model, p);
        
        se    = sqrt(diag(VCV));
        se(4) = se(4)/sqrt(n);
        
    end
       
end
