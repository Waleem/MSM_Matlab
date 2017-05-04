function se = MSM_std_err(para,kbar,n ,data, A_template)

%%This algorithm computes the standard error from univariate MSM using the
% second-derivative estimate of the  information matrix. See Hamilton, Time Series Analysis,page 143.

% This code benefits heavily from the robustvcv.m code of the MFEToolbox by 
% Kevin Sheppard (% kevin.sheppard@economics.ox.ac.uk)

%Input:
%      x    = MSM parameter estimates vector [b m gamma sigma]
%      k    = number of volatility components i.e in MSM(k)
%      data = data vector i.e returns

%Output: 
%      s= A 4x1 standard error vector


%Checks:
    if size(data,2) > 1
        data = data';
    end
    if size(data,1) == 1
        error('data series must be a vector')
    elseif isempty(data)
        error('data is empty')
    end
    
    if size(para,2) > 1
        para=para';
    end
    if size(para,1) == 1||size(para,1)>4
        error('para must be a 4-by-1 vector')
    end
    
    if (kbar < 1) 
    error('k must be a positive integer')
    end
    kbar = floor(kbar);
    
    T=size(data,1);
    
  
    if kbar==1
        gt1 = MSM_grad(@MSM_LLs,para,kbar,n,data,A_template);
        gt1 = gt1(:,2:end); %Exclude b from the estimation of covariance matrix because b is not defined at k=1.
                            %Else the standard errors for all the parameters will be inaccurate or may even return Inf
        J1    = gt1'*gt1;
        se    = [nan; sqrt(diag(inv(J1)))];    
    else
        
        VCV = robustvcv('get_MSM_likelihood',para,2,kbar,n, data, A_template);
        
        se   = sqrt(diag(VCV));
        
    end
    
    %s = se';
end

%%%%%%%%This code is fine. But the above code is more efficient%%%%
%     gt1       = MSM_grad(@MSM_LLs,para,kbar,data);
%     if kbar==1
%         
%         gt1=gt1(:,2:end); %Exclude b from the estimation of covariance matrix because b is not defined at k=1.
%                           %Else the standard errors for all the parameters will be inaccurate or may even return Inf
%         J1    = gt1'*gt1;
%         se    = [nan; sqrt(diag(inv(J1)))];    
%     else
%         J1    = gt1'*gt1;
%         H1        = MSM_hessian(@MSM_LL,para,kbar,data);
%         [~,I]     = chol(H1);
%         if I~=0 %if H1 (the negative of Hessian) is not positive-definite (this usually happens when m0 or gamma_k takes boundary values)
%                 % use the Outer Product Gradient for standard errors.
%             se    = sqrt(diag(inv(J1)));
%         else
%             se    = sqrt(diag(inv(H1/T)*J1/T^2*inv(H1/T)));
%         end
%         
%     end