function [LL,LLs1, LLs2,pi_mat1,pi_mat2,A] = BMSM_2stage_likelihood1(para,kbar,data,A_template,n)
    %%Likelihood calculation - Computes the first stage log-likelihood for the
    %%2-stage bivariate MSM.
    % The log-likelihhod computed is also equivalent to the log-likelihood for
    % the combined univariate MSM
    %
    %Input:  
    %         Para        -  6-by-1 vector of parameters [m01, m02, sigma1,sigma2, gamma_k, b]
    %         Kbar        -  Number of volatility components (scalar)
    %         Data        -  T-by-2 matrix of zero-mean returns
    %         A_template  -  Initial template for the univariate MSM transition matrix (2^kbar-by-2^kbar)
    %         n           -  The number of trading periods in 1 year. E.g
    %                        n= 252 for the number of business days in a year. n=12 for monthly data
    %Output:  LL          -  sum of log-likelihood for first stage of a 2-stage bivariate MSM
    %         LLs1        -  Individual daily log-likelihoods at optimum for d data_alpha (T-by-1)
    %         LLs2        -  Individual daily log-likelihoods at optimum for d data_beta (T-by-1)
    %         pi_mat1     -  Filtered state probabilities at optimum for data^alpha (T+1-by-2^kbar)
    %         pi_mat2     -  Filtered state probabilities at optimum for data^beta (T+1-by-2^kbar)
    %         A           -  Transition matrix at optimum (2^kbar-by-2^kbar)

    m01=para(1);
    m02=para(2);
    sigma1 = para(3)/sqrt(n);
    sigma2 = para(4)/sqrt(n);
    if para(5)>=1
       para(5)=0.9999; %This has nothing to do with the estimation algorithm itself, as the condition will always be satisfied.
                       % It is rather for the algorithm used to estimate the GMM standard errors. Sometimes the estimated
                       % gamma_k may take boundary, or close to boundary values such 0.9999. This is fine. But when
                       % computing the gradient or hessian of the objective function, the algorithm may attempt to evaluate 
                       % the objective function  at a point where gamma_k >1, which violates the MSM model constraint that gamma_k =(0,1). 
                       % That is when these lines of code kick in to stop such violation.
    end
    gamma_k=para(5);
    b=para(6);
    k =2^kbar;

    A = transition_mat(A_template,gamma_k,b,kbar);

    %g_m1 = gofm(m01,kbar);
    %g_m2 = gofm(m02,kbar);
    [~, ~,g_m1] = Univ_MSM_states(m01,kbar);
    [~, ~,g_m2] = Univ_MSM_states(m02,kbar);
    
    T = size(data,1);

    LLs1=zeros(1,T);
    LLs2=LLs1;
    
    pi_mat1 = zeros(T+1,k);  
    pi_mat1(1,:) = (1/k)*ones(1,k);
    pi_mat2 = pi_mat1; 

    %*----------------------------------------------------------------------*
    %*                        Likelihood algorithm                          *
    %*----------------------------------------------------------------------*
    pa = (2*pi)^-0.5;
    s1 = repmat(sigma1*g_m1,T,1);
    s2 = repmat(sigma2*g_m2,T,1);
    w_t1 = repmat(data(:,1),1,k);
    w_t1 = pa*exp( - 0.5.*((w_t1./s1).^2))./s1; 
    w_t1 = w_t1 + 1e-16;
    w_t2 = repmat(data(:,2),1,k);
    w_t2 = pa*exp( - 0.5.*((w_t2./s2).^2))./s2; 
    w_t2 = w_t2 + 1e-16;

    pmat1 = pi_mat1;
    pmat2 = pi_mat2;
    
    %ARX_MSM_Core is a mex function. The compiled copy included in the toolbox was compiled
    %on a 64-bit win8 pc. This may work for you. If the compiled copy does not work for you,
    %please compile the included c-code before you use. 
    %Else, uncomment and use the loop below. But the loop may be very slow for a very large data
    %and kbar.
    
    [pi_mat1, LLs1]=ARX_MSM_Core(pmat1,A,w_t1); 
    [pi_mat2, LLs2]=ARX_MSM_Core(pmat2,A,w_t2); 

    % for t=2:T+1          
        % piA1 = (pi_mat1(t-1,:)*A);
        % piA2 = (pi_mat2(t-1,:)*A);
        % C1 = (w_t1(t-1,:).*piA1); 
        % ft1 = sum(C1);
        % if ft1 == 0                      %This stop div by zero if probs are too low
            % pi_mat1(t,1) = 1;   
        % else
            % pi_mat1(t,:) = C1 / ft1; 
        % end
    % 
        % C2 = (w_t2(t-1,:).*piA2); 
        % ft2 = sum(C2);
        % if ft2 == 0                      %This stop div by zero if probs are too low
            % pi_mat2(t,1) = 1;   
        % else
            % pi_mat2(t,:) = C2 / ft2; 
        % end
    % 
        % LLs1(t-1) = log(dot(w_t1(t-1,:),piA1));
        % LLs2(t-1) = log(dot(w_t2(t-1,:),piA2));
    % end 
        LL1=-sum(LLs1);
        LL2=-sum(LLs2);
        LL=LL1+LL2;
        if ~isfinite(LL)
            disp('Log-likelihood is inf.')
        end

    
end

% Calculate the transition matrix using the template.
% This code is from the original MSM_MLE code by Calvet.
function A = transition_mat(A,gamma_kbar,b,kbar)
        
    gamma = zeros(kbar,1);                          
    gamma(1) = 1-(1-gamma_kbar)^(1/(b^(kbar-1)));
    for i = 2:(kbar)
        gamma(i,1) = 1-(1-gamma(1,1))^(b^(i-1));
    end
    gamma = gamma*0.5;
    gamma(:,2)=gamma(:,1);
    gamma(:,1) = 1 - gamma(:,1);  
    kbar1 = kbar +1;
    kbar2 = 2^kbar;
    prob = ones(kbar2,1);    
    
    for i=0:2^kbar-1    %Works out probability associated with each XOR number
        for m = 1:kbar  
            prob(i+1,1) = prob(i+1,1) * gamma(kbar1-m, (bitget(i,m)+1));
        end
    end
    for i =0:2^(kbar-1)-1   %Copies probabilities to the transition matrix
        for j = i:(2^(kbar-1)-1)  
            A(kbar2-i,j+1)     =  prob(kbar2-A(i+1,j+1),1);%Copies each probability to the other 8 symmetrical locations
            A(kbar2-j,i+1)     =  A(kbar2-i,j+1);
            A(j+1,kbar2-i)     =  A(kbar2-i,j+1);
            A(i+1,kbar2-j)     =  A(kbar2-i,j+1);    
            A(i+1,j+1)         = prob(A(i+1,j+1)+1,1);
            A(j+1,i+1)         = A(i+1,j+1);
            A(kbar2-j,kbar2-i) = A(i+1,j+1);
            A(kbar2-i,kbar2-j) = A(i+1,j+1);
        end
    end 
end

% Calculate all of the possible volatility states. This is the Calvet code included in the original MSM_MLE
% function g_m = gofm(m0,kbar)
% 
% m1=2-m0;
% kbar2 = 2^kbar;
% g_m1 = [0:(kbar2-1)];
% 
% for i = 1:(kbar2)
    % g=1;
    % for j = 0:(kbar-1)       
        % if(bitand(g_m1(i),(2^j))~=0)    %
            % g=g*m1;
        % else g=g*m0;
        % end
    % end
    % g_m1(i)=g;
% end
% 
% g_m=sqrt(g_m1);
% end

