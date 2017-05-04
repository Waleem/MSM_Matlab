function [LL,LLs] = RES_BMSM_2stage_likelihood2(para,kbar,data,para1,n)
    %%Likelihood calculation - Computes the second stage log-likelihood for
    %%the 2-stage restricted bivariate MSM. In the restricted bivariate MSM
    %%model, rho_m is assumed to be 1.
    
    % This code benefits a great deal, from the KronProd toolbox of Matt J. 
    % The KronProd toolbox makes it possible to construct the transition matrix_type
    % without actually expanding it into a 4^kbar-by-4^kbar matrix.
    % For kbar >5, expanding the transition matrix will require a large amount of memory,
    % in gigabytes.
    % Without the toolbox, it will be difficult, if not impossible, to estimate
    % bivariate MSM with more than 5 volatility component.
    % The toolbox is included in the MSM toolbox, and is also available here
    % http://www.mathworks.com/matlabcentral/fileexchange/25969-efficient-object-oriented-kronecker-product-manipulation
    
    
    %Input:  
    %         Para        -  2-by-1 vector of second-stage parameters [rho_e lambda]
    %         Kbar        -  Number of volatility components (scalar)
    %         Data        -  T-by-2 matrix of zero-mean returns
    %         para1       -  First-stage parameter estimates [m01, m02, sigma1,sigma2, gamma_k, b]
    %         n           -  The number of trading periods in 1 year. E.g
    %                        n= 252 for the number of business days in a year. n=12 for monthly data
    %Output:  LL          -  sum of log-likelihood for second stage of a 2-stage bivariate MSM
    %         LLs         -  Individual daily log-likelihoods at optimum (T-by-1)
    %         

    
    rho_m=1;
    %Gather paremeters
    rho_e=para(1);
    lamda=para(2);
    
    if rho_e >=1
       rho_e=0.9999;   %This has nothing to do with the estimation algorithm itself, as the condition will always be satisfied.
                       % It is rather for the algorithm used to estimate the GMM standard errors. Sometimes the estimated
                       % rho_e may take boundary, or close to boundary values such as 0.9999. This is fine. But when
                       % computing the gradient or hessian of the objective function, the algorithm may attempt to evaluate 
                       % the objective function  at a point where rho_e >1, which violates the MSM model constraint that rho_e =[0,1]. 
                       % That is when these lines of code kick in to stop such violation.
    end
    if lamda>=1
       lamda=0.9999; 
    end
    
    k2 =4^kbar;
    %Parameters from first-stage estimation
    m01=para1(1);
    m02=para1(2);
    sigma1 = para1(3)/sqrt(n);
    sigma2 = para1(4)/sqrt(n);
    gamma_kbar = para1(5);
    b = para1(6);
    
    %Get switching probabilities
    gamma = zeros(kbar,1);                          
    gamma(1) = 1-(1-gamma_kbar)^(1/(b^(kbar-1)));
    for i = 2:(kbar)
        gamma(i,1) = 1-(1-gamma(1,1))^(b^(i-1));
    end
    %Compute transition matrices and gather them into opset array. There will be kbar transition matrices.
    opset=arrayfun(@(x,y,z) BMSM_Transition(x,y,z),gamma,lamda*ones(kbar,1),rho_m*ones(kbar,1),'UniformOutput', false);
    %Get kronecker product of all the transition matrices.
    A = KronProd(opset,1:kbar,4*ones(kbar,1),1);
    % A is a KronProd object. To see the content of A, use
    % full(A). Caution!!! You may run out of memory

    T = size(data,1); 
    LLs = zeros(1,T);

    g_m = BMSM_states(m01,m02,kbar);
    
%     May be used, but not stable in my testing
%     pi_mat = zeros(T+1,k2); 
%     pi_mat(1,:)=BMSM_ergodic(gamma,lamda,rho_m,kbar);
%     matObj = matfile('myBigData.mat','Writable',true);
%     matObj.pi_mat=pi_mat;
%     clear pi_mat
        
    piA=BMSM_ergodic(gamma,lamda,rho_m,kbar)*A;
    pa = 1/(2*pi*sqrt(1-rho_e^2));

    %*----------------------------------------------------------------------*
    %*                        Likelihood algorithm  begins                  *
    %*----------------------------------------------------------------------*
for t=2:T+1          
       
    s1=(data(t-1,1)./(sigma1*g_m(1,:))).^2;
    s2=(data(t-1,2)./(sigma2*g_m(2,:))).^2;
    s=(data(t-1,1)./(sigma1*g_m(1,:))).*(data(t-1,2)./(sigma2*g_m(2,:)));
    z= s1+s2-(2*rho_e*s);
    w_t_1=pa*exp(-z/(2*(1-rho_e^2)))./((sigma1*g_m(1,:)).*(sigma2*g_m(2,:))) +  1e-16;
    
    C=w_t_1.*piA;
    ft = sum(C);
%     Requires lines 48 to 53
%     if ft == 0
%         matObj.pi_mat(t,1)=1;
%     else
%         matObj.pi_mat(t,:)=C / ft;
%     end

    LLs(t-1) = log(dot(w_t_1,piA));
    if ft == 0
        piA=[1 zeros(k2-1)]*A;
    else
        piA=(C / ft)*A;
    end

end 
    LL=-sum(LLs);
    if ~isfinite(LL)
        disp('Log-likelihood is inf. Probably due to all zeros in pi_mat.')
    end
end
