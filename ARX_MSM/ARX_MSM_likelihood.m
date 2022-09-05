function [LL, LLs, A, p_mat, M, g_m, et, yhat] = ARX_MSM_likelihood(para, kbar, Y, xmat, n, model, p)
%%Likelihood calculation - Computes the log-likelihood for ARX_MSM models.

%Input:  
%         Para        -  3-by-1 vector of second-stage parameters [rho_e lambda rho_m]
%         Kbar        -  Number of volatility components (scalar)
%         Data        -  T-by-2 matrix of zero-mean returns
%         para1       -  First-stage parameter estimates [m01, m02, sigma1,sigma2, gamma_k, b]
%         n           -  The number of trading periods in 1 year. E.g
%                        n= 252 for the number of business days in a year. n=12 for monthly data
%Output:  LL          -  sum of log-likelihood for second stage of a 2-stage bivariate MSM
%         LLs         -  Individual daily log-likelihoods at optimum (T-by-1)
%         


%Gather MSM paremeters
b          = para(1);
m0         = para(2);
gamma_kbar = para(3);
sigma      = para(4)/sqrt(n);

k2 =2^kbar;


A = ARX_MSM_Transition_Mat(b, gamma_kbar, kbar);
% Get state vector:
[M, ~,g_m] = Univ_MSM_states(m0,kbar);

T = size(Y,1); 
LLs = zeros(1,T);

p_mat      = zeros(T+1,k2);
p_mat(1,:) = (1/k2)*ones(1,k2);


%*----------------------------------------------------------------------*
%*                        Likelihood algorithm  begins                        *
%*----------------------------------------------------------------------*

if size(para,2)>1
    para=para';
end
switch model
    case 1 %ARX-MSM
        
        mean_para  = para(5:end);
        et         = Y - xmat*mean_para;
        s   = repmat(sigma*g_m,T,1);
        w_t = repmat(et,1,k2);
        w_t = ((2*pi)^-0.5) * exp( - 0.5.*((w_t./s).^2))./s; 
        w_t = w_t + 1e-16;
        
    case 2 %ARX-MSM_M (Regime dependent AR parameters)
        bm = para(6:6+numel(p)-1)*sum(M-1,2)';
        et = repmat(Y - para(5).*xmat(:,1),1,k2);
        if size(xmat,2) == numel(p)+1
           ylag=xmat(:,2:end);
           for i = 1:numel(p)
               et = et - bsxfun(@times, ylag(:,i),bm(i,:));
           end
        else
           ylag=xmat(:,2:2+numel(p)-1);
           for i = 1:numel(p)
               et = et - bsxfun(@times, ylag(:,i),bm(i,:));
           end
           et = et - repmat(xmat(:,numel(p)+2:end)*para(numel(p)+6:end),1,k2);
        end
        s   = repmat(sigma*g_m,T,1);
        w_t = ((2*pi)^-0.5) * exp( - 0.5.*((et./s).^2))./s; 
        w_t = w_t + 1e-16;
        
    case 3
        et  = Y - xmat*para(5:end-1);
        s   = repmat(sigma*g_m,T,1);
        w_t = repmat(et,1,k2) - para(end)*s;
        w_t = ((2*pi)^-0.5) * exp( - 0.5.*((w_t./s).^2))./s; 
        w_t = w_t + 1e-16;
end

pmat=p_mat; 
%ARX_MSM_Core is a mex function. Please compile before you use. Else, uncomment and
%use the loop below. But the loop may be very slow for a very large data
%and kbar.
[p_mat, LLs]=ARX_MSM_Core(pmat,A,w_t);      
% for t=2:T+1          
%     piA = (p_mat(t-1,:)*A);
%     C   = (w_t(t-1,:).*piA); 
%     ft  = sum(C);
%     if ft == 0                      %This stop div by zero if probs are too low
%         p_mat(t,1) = 1;   
%     else
%         p_mat(t,:) = C / ft; 
%     end
%     
%     LLs(t-1) = log(dot(w_t(t-1,:),piA));
% end

LL = -sum(LLs);
if ~ isfinite(LL)
    disp('Log-likelihood is inf. Probably due to all zeros in pi_mat.')
end

if nargout >6
    switch model
        case 1 %ARX-MSM
            yhat       = xmat*mean_para;
            et         = Y - yhat;
        case 2 %ARX-MSM_M (Regime dependent AR parameters)
            %Get the marginal probability for each M_k i.e prob[M_k = m0]
            pm   = MSM_marginals(p_mat(2:end,:),m0,M,kbar);
            %Get expected value of Mk, for each Mk
            em   = m0*pm + (2-m0)*(1-pm);
            bm   = sum(em - 1,2)*para(6:6+numel(p)-1)';
           
            if size(xmat,2) == numel(p)+1
               ylag = xmat(:,2:end);
               yhat = para(5)*xmat(:,1) + sum(bm.*ylag,2);
            else
               ylag = xmat(:,2:2+numel(p)-1);
               yhat = para(5)*xmat(:,1) + sum(bm.*ylag,2) + xmat(:,numel(p)+2:end)*para(numel(p)+6:end);
            end
            et = Y-yhat;
        case 3
            mean_para  = para(5:end);
            ht   = sigma*(p_mat*g_m');
            ht   = ht(2:end,1);
            yhat = [xmat ht]*mean_para;
            et     = Y - yhat;
    end
end
end


% %Get switching probabilities for each Mk
% gamma = zeros(kbar,1);                          
% gamma(1) = 1-(1-gamma_kbar)^(1/(b^(kbar-1)));
% for i = 2:(kbar)
%     gamma(i,1) = 1-(1-gamma(1,1))^(b^(i-1));
% end
% The approach used to compute transition matrix here is different from the Calvet(2004) approach.
% This approach exploits the independence of the volatility components.
% Each Mk has its own 2-by-2 transition matrix Tk. It can be shown that the
% transition matrix for M=[M1;M2,..Mkbar] is equal to the kronecker product
% of all the transition matrices for all the Mk. 
% This approach is 100 times faster in my test, using Kbar=10. Combine
% this with the use of the KronProd toolbox, this approach also requires
% significantly less memory.

% % Compute the transition matrices for each Mk and gather them into opset array. There will be kbar transition matrices.
% opset = arrayfun(@(x) ARX_MSM_Transition(x),gamma,'UniformOutput', false);
% % Get kronecker product of all the transition matrices. Use KronProd because it
% % is faster and requires less memory.
% A = KronProd(opset,(kbar:-1:1),2*ones(kbar,1),1);

