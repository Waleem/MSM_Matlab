function A = ARX_MSM_Transition_Mat(b, gamma_kbar, kbar)
% Calculate the transition matrix  for univariate MSM
    A = zeros((2^kbar),(2^kbar));
    for i =0:2^kbar-1       
        for j = i:(2^kbar-1)-i  
            A(i+1,j+1) = bitxor(i,j);
        end
    end
    
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