function [startingvals,LLs,output_parameters] =...
    MSM_starting_values(data,startingvals,kbar,n ,A_template)

options=optimset('fminbnd');
options.TolX=1e-3;


LLs=[];
output_parameters=[];
if isempty(startingvals)
    disp('No starting values entered: Using grid-search')
    % A grid search used to find the best set of starting values in the
    % event none are supplied by the user
    
    b=[1.5 3 6 20];
    lb=length(b);
    g=[.1 .5 .9 ];
    lg=length(g);
   
    sigma = std(data)*sqrt(n);
    output_parameters=zeros(lb*lg,3);
    LLs=zeros(lb*lg,1);
    m0_lower = 1.2;         m0_upper = 1.8;     

    index=1;
    for i=1:lb         
            for j=1:lg
                temp_in = [b(i);1.5;g(j);sigma];
                [m0,LL]=fminbnd('MSM_likelihood',m0_lower,m0_upper,options,kbar,n,data,A_template,[b(i);g(j);sigma]);
                parameters=[b(i);m0;g(j)];
                output_parameters(index,:)=parameters';
                LLs(index)=LL;
                index=index+1;
            end
    end
    [LLs,index]=sort(LLs);
    startingvals=[output_parameters(index(1),:),sigma]';
    output_parameters=output_parameters(index,:);
end