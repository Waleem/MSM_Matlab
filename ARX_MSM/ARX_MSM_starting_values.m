function startingvals = ARX_MSM_starting_values(Y, xmat, kbar, n, model)

options1=optimset('fminbnd');
options1.TolX=1e-3;


    
% A grid search used to find the best set of starting values in the
% event none are supplied by the user
% A clever way to get starting values. Estimate ARX to get errors. Estimate MSM(k) to get MSM parameters.
% If model is MSM-in-mean, then use ht from MSM to estimate ARX again to get approximate MSM-in-mean

% Step 1, fit Ordinary Least Square and get errors as follows: 
 mean_parameters  = (xmat'*xmat)\xmat'*Y;
 et               = Y - xmat*mean_parameters;

if model ==3
    options2 = get_options();
    msm_para = fmincon('MSM2_likelihood', [5 1.5 .5 std(Y)*sqrt(n)], [], [], [],[],...
              [1, 1, 0.0001, 0.0001],[50, 1.9999, 0.9999, 50],[],options2, kbar, et, n);
    
    [~,~, ~, p_mat, ~, g_m] = MSM2_likelihood(msm_para, kbar, et,n);
    ht               = (msm_para(4)/sqrt(n))*(p_mat*g_m');
    ht=ht(2:end,1);
    mean_parameters  = ([xmat ht]'*[xmat ht])\[xmat ht]'*Y;
    et               = Y - [xmat ht]*mean_parameters;
end
    
    %Step2, use grid search for the MSM parameters
    
    b     = [1.5 3 6 20];
    lb    = length(b);
    g     = [.1 .5 .9  ];
    lg    = length(g);
    sigma = std(et)*sqrt(n);

    msm_parameters = zeros(lb*lg,4);
    LLs            = zeros(lb*lg,1);
    m0_lower       = 1.2;         
    m0_upper       = 1.8; 

    index=1;
    for i=1:lb         
        for j=1:lg
            
            [m0,LL]=fminbnd('MSM2_likelihood',m0_lower,m0_upper,options1,kbar,et,n,[b(i) g(j) sigma]);
            msm_parameters(index,:)=[b(i) m0 g(j) sigma];
            LLs(index) = LL;
            index      = index+1;
        end
    end
    [~,index1]=sort(LLs);

   
    startingvals=[msm_parameters(index1(1),:) mean_parameters'];
         

    function options = get_options()
        options  =  optimset('fmincon');
        options  =  optimset(options , 'TolFun'      , 1e-005);
        options  =  optimset(options , 'TolX'        , 1e-005);
        options  =  optimset(options , 'Display'     , 'off');
        options  =  optimset(options , 'Diagnostics' , 'on');
        options  =  optimset(options , 'LargeScale'  , 'off');
        options  =  optimset(options , 'MaxFunEvals' , 2000);
        options  =  optimset(options , 'MaxSQPIter' , 2000);
        options  =  optimset(options , 'Algorithm'   ,'sqp');
        options  =  optimset(options , 'UseParallel' , 'Always');
    end

end