function se = VECM_M_BMSM_std_err(r, S, F, para, kbar, n)
    %INPUTS:  
    %         Para        -  12-by-1 vector of parameters 
    %                              [m01, m02, sigma1,sigma2, gamma_k, b, bs, bf, as, af, rho_e, lambda]
    %         Kbar        -  Number of volatility components (scalar)
    %         r           -  T-by-2 matrix of returns [log_spot_return log_futures_return]
    %         S           -  T-by-1 vector of log_spot prices
    %         F           -  T-by-1 vector of log_futures prices
    %         A_template  -  Initial template for the univariate MSM transition matrix (2^kbar-by-2^kbar)
    %         n           -  The number of trading periods in 1 year. E.g
    %                        n= 252 for the number of business days in a year. n=12 for monthly data
    % OUTPUTS:
    %    se               -  A 12x1 row vector of standard errors for parameters 
    %                            [m01, m02, sigma1,sigma2, gamma_k, b, bs, bf, as, af, rho_e, lambda]
    
    %Note that the hassian matrix for the 2 stage VECM bivariate MSM is a block
    %diagonal matrix 12-by-12 of [H1 0_2]
    %                          [0_10 H2], where H1 is a 10-by-10 matrix, H2 is
    %                          a 2-by-2 matrix, 0_2 is a 10-by-2 matrix and
    %                          0_6 is a 2-by-10 matrix. Therefor, the
    %                          standard errors for the first 10 parameters
    %                          [m01, m02, sigma1, sigma2, gamma_kbar, b, bs, bf, as, af] are estimated separately and 
    %                          the standard errors for the last 2
    %                          parameters [rho_e, lambda] are also estimated separately.
    %
    %References:  Volatility comovement: a multifrequency approach. 
    %             Journal of Econometrics, 131(1-2):179{215.
    % ------------------------------------------------------------------------- 
    
    %Checks:
    if size(r,2) > 2
        r = r';
    end
    if size(r,1) == 1
        error('return matrix must be T-by-2')
    elseif isempty(r)
        error('return matrix is empty')
    end
    
    if size(para,2) > 1
        para=para';
    end
    if size(para,1) == 1 ||size(para,1)>12
        error('para must be a 12-by-1 vector')
    end
    
    if (kbar < 1) 
    error('k must be a positive integer')
    end
    kbar = floor(kbar);
    
    T=size(r,1);
    
    para1=para(1:10);
    para2=para(11:12);
    
    %Stage 1 parameter standard errors:
    A = zeros((2^kbar),(2^kbar));
    for i =0:2^kbar-1       
        for j = i:(2^kbar-1)-i  
            A(i+1,j+1) = bitxor(i,j);
        end
    end
    
    if kbar==1
        %Exclude b from the estimation of covariance matrix because b is not defined at k=1.
        %Else the standard errors for all the stage 1 parameters will be inaccurate or may even return Inf
        
        gt1       = MSM_grad(@VECM_M_BMSM_LLs1, para1, kbar, r, S, F, A, n);
        gt1(:,6)  =[] ;
        J1        = gt1'*gt1;
        H1        = MSM_hessian(@VECM_M_BMSM_likelihood1, para1, kbar, r, S, F, A, n);
        H1(6,:)   = [];
        H1(:,6)   = [];
        [~,I]     = chol(H1);
                    %if H1 (the negative of Hessian) is not positive-definite (this usually happens when m0 or gamma_k takes boundary values)
                    % use the Outer Product Gradient for standard errors.
        if I~=0
            se1   = sqrt(diag(inv(J1)));
        else
            se1   = sqrt(diag(inv(H1/T)*J1/T^2*inv(H1/T)));
        end
        se1       = [se1(1:5);nan;se1(6:end)];
    else
        gt1       = MSM_grad(@VECM_BMSM_LLs1, para1, kbar, r, S, F, A, n);
        J1        = gt1'*gt1;
        H1        = MSM_hessian(@VECM_BMSM_likelihood1, para1, kbar, r, S, F, A, n);
        [~,I]     = chol(H1);
        if I~=0
            se1   = sqrt(diag(inv(J1)));
        else
            se1   = sqrt(diag(inv(H1/T)*J1/T^2*inv(H1/T)));
        end
    end
       
    %Stage 2 parameter standard errors:
    
    % Get error vector from the VECM model and  for estimating rho_e and lambda
    
 % Get the error correction term from the VECM model for estimating rho_e and lambda
    et = VECM_M_BMSM_ect(para1,r,S,F,kbar,n);
        
    gt2       = MSM_grad(@RES_BMSM_2stage_LLs2,para2,para1,et,kbar,n);
    J2        = gt2'*gt2;
    H2        = MSM_hessian(@RES_BMSM_2stage_likelihood2,para2,kbar,et,para1,n);
    [~,I]     = chol(H2);
    if I~=0
        se2   = sqrt(diag(inv(J2)));
    else
        se2   = sqrt(diag(inv(H2/T)*J2/T^2*inv(H2/T)));
    end
        
  
    se=[se1' se2'];
   
end
