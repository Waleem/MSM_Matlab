function [LL, LLs1, LLs2, pi_mat1, pi_mat2, A] = VECM_M_BMSM_likelihood1(para, kbar, r, S, F, A_template,n)
    %%Likelihood calculation - Computes the first stage log-likelihood for the
    %%Vector Error Correction Bivariate Markov Switching Multifractal Model with Time Varying Speed 
    % of Long-Run Adjustment (VECM(M)-MSM), assuming that rho_m (the correlation between M^alpha and M^beta) is 1.
    %                                         R_s = b_s + a_s((Ms_1 -1 + Ms_2 -1 +. ...+ Ms_kbar-1))*(S-F) + e_s
    %                                         R_f = b_f + a_f((Mf_1 -1 + Mf_2 -1 +. ...+ Mf_kbar-1)-1)*(S-F) + e_f
    %                                         e_s = sigma_s*(Ms_1*Ms_2*. ...*Ms_kbar)
    %                                         e_f = sigma_f*(Mf_1*Mf_2*. ...*Mf_kbar)
    %Input:  
    %         Para        -  10-by-1 vector of parameters [m01, m02, sigma1,sigma2, gamma_k, b, bs, bf, as, af]
    %         Kbar        -  Number of volatility components (scalar)
    %         r           -  T-by-2 matrix of returns [log_spot_return log_futures_return]
    %         S           -  T-by-1 vector of log_spot prices
    %         F           -  T-by-1 vector of log_futures prices
    %         A_template  -  Initial template for the univariate MSM transition matrix (2^kbar-by-2^kbar)
    %         n           -  The number of trading periods in 1 year. E.g
    %                        n= 252 for the number of business days in a year. n=12 for monthly data
    %
    %Output:  LL          -  sum of log-likelihood for first stage of a 2-stage  VECM Bivariate MSM
    %         LLs1        -  Individual daily log-likelihoods at optimum for spot (T-by-1)
    %         LLs2        -  Individual daily log-likelihoods at optimum for futures (T-by-1)
    %         pi_mat1     -  Filtered state probabilities at optimum for spot (T+1-by-2^kbar)
    %         pi_mat2     -  Filtered state probabilities at optimum for futures (T+1-by-2^kbar)
    %         A           -  Transition matrix at optimum (2^kbar-by-2^kbar)

    m01    =para(1);
    m02    =para(2);
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
    gamma_k =para(5);
    b       =para(6);
    bs      =para(7);
    bf      =para(8);
    as      =para(9);
    af      =para(10);
    
    T = size(r,1);
    k =2^kbar;
    
    %Get BMSM state vectors
    [Mmat1, ~,g_m1] = Univ_MSM_states(m01,kbar);
    [Mmat2, ~,g_m2] = Univ_MSM_states(m02,kbar);
    Mmat1=sum(Mmat1-1,2)';
    Mmat2=sum(Mmat2-1,2)';
   
    % Get the T-by-2^kbar error vectors from the VECM model
    X1 = bs*ones(T,k) + (S-F)*(as*(Mmat1));
    X2 = bf*ones(T,k) + (S-F)*(af*(Mmat2));
            
    %Begin the MSM model
    A = transition_mat(A_template,gamma_k,b,kbar);

    
    pi_mat1 = zeros(T+1,k);  
    pi_mat1(1,:) = (1/k)*ones(1,k);
    pi_mat2 = pi_mat1; 
    LLs1=zeros(1,T);
    LLs2=LLs1;

    %*----------------------------------------------------------------------*
    %*                        Likelihood algorithm                          *
    %*----------------------------------------------------------------------*
    pa = (2*pi)^-0.5;
    s1 = repmat(sigma1*g_m1,T,1);
    s2 = repmat(sigma2*g_m2,T,1);
    w_t1 = repmat(r(:,1),1,k)-X1;
    w_t1 = pa*exp( - 0.5.*((w_t1./s1).^2))./s1; 
    w_t1 = w_t1 + 1e-16;
    w_t2 = repmat(r(:,2),1,k)-X2;
    w_t2 = pa*exp( - 0.5.*((w_t2./s2).^2))./s2; 
    w_t2 = w_t2 + 1e-16;


    for t=2:T+1          
        piA1 = (pi_mat1(t-1,:)*A);
        piA2 = (pi_mat2(t-1,:)*A);
        C1 = (w_t1(t-1,:).*piA1); 
        ft1 = sum(C1);
        if ft1 == 0                      %This stop div by zero if probs are too low
            pi_mat1(t,1) = 1;   
        else
            pi_mat1(t,:) = C1 / ft1; 
        end

        C2 = (w_t2(t-1,:).*piA2); 
        ft2 = sum(C2);
        if ft2 == 0                      %This stop div by zero if probs are too low
            pi_mat2(t,1) = 1;   
        else
            pi_mat2(t,:) = C2 / ft2; 
        end

        LLs1(t-1) = log(dot(w_t1(t-1,:),piA1));
        LLs2(t-1) = log(dot(w_t2(t-1,:),piA2));
    end 
        LL1=-sum(LLs1);
        LL2=-sum(LLs2);
        LL=LL1+LL2;
        if ~isfinite(LL)
            disp('Log-likelihood is inf. Probably due to all zeros in pi_mat.')
        end

    
end

% Calculate the transition matrix using the template
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
            A(kbar2-i,j+1) = prob(kbar2-A(i+1,j+1),1);%Copies each probability to the other 8 symmetrical locations
            A(kbar2-j,i+1) =  A(kbar2-i,j+1);
            A(j+1,kbar2-i) =  A(kbar2-i,j+1);
            A(i+1,kbar2-j) =  A(kbar2-i,j+1);    
            A(i+1,j+1) = prob(A(i+1,j+1)+1,1);
            A(j+1,i+1) = A(i+1,j+1);
            A(kbar2-j,kbar2-i) = A(i+1,j+1);
            A(kbar2-i,kbar2-j) = A(i+1,j+1);
        end
    end 
end



