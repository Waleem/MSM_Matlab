function startingvals = VECM_BMSM_starting_values(r, S, F, startingvals, kbar, A_template, n)

options=optimset('fminbnd');
options.TolX=1e-3;


if isempty(startingvals)
    disp('No starting values entered: Using grid-search')
    % A grid search used to find the best set of starting values in the
    % event none are supplied by the user
          
    % Step 1, fit Ordinary Least Square as follows: R_s = b_s + a_s(S-F) + e_s
    %                                               R_f = b_f + a_f(S-F) + e_f
    rs   = r(:,1);
    rf   = r(:,2);
    x    = [ones(size(rs,1),1) S-F];    
    Bs   = (x'*x)\x'*rs;
    Bf   = (x'*x)\x'*rf;
    
    et = [(rs - x*Bs)  (rf - x*Bf)];
    
    %Step2, use MSM grid search for the MSM parameters
    
    b=[1.5 3 6 20];
    lb=length(b);
    g=[.1 .5 .9  ];
    lg=length(g);
    sigma1 = std(et(:,1))*sqrt(n);
    sigma2 = std(et(:,2))*sqrt(n);
    output_parameters=zeros(lb*lg,6);
    LLs1=zeros(lb*lg,1);
    LLs2=zeros(lb*lg,1);
    m0_lower = 1.2;         m0_upper = 1.8; 

    index=1;
    for i=1:lb         
        for j=1:lg
            
            [m01,LL1]=fminbnd('MSM_likelihood',m0_lower,m0_upper,options,kbar,et(:,1),A_template,[b(i);g(j);sigma1]);
            [m02,LL2]=fminbnd('MSM_likelihood',m0_lower,m0_upper,options,kbar,et(:,2),A_template,[b(i);g(j);sigma2]);
            
            output_parameters(index,:)=[m01 m02 sigma1 sigma2 g(j) b(i)];
            LLs1(index)=LL1;
            LLs2(index)=LL2;
            index=index+1;
        end
    end
    [LLs1,index1]=sort(LLs1);
    [LLs2,index2]=sort(LLs2);
    bg1=output_parameters(index1(1),5:6);
    bg2=output_parameters(index2(1),5:6);
    bg=(bg1+bg2)/2;
    
    m01 = output_parameters(index1(1),1);
    m02 = output_parameters(index2(1),2);
    c=corr(et);
    
    startingvals=[m01 m02 sigma1 sigma2 bg Bs(1) Bf(1) Bs(2) Bf(2) c(1,2) c(1,2)];
         
end

end