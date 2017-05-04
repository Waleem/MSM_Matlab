%--------------------------------------------------------------------------
% Original code modified by Waleem Alausa
% This version: 25/Nov/2011
%--------------------------------------------------------------------------
%Likelihood calculation - uses method from Calvet&Fisher 2004

function [LL,LLs,A,g_m,pi_mat] = MSM_likelihood(input,kbar,n, data,A_template,estim_flag)

if length(input) ==1
    input = [estim_flag(1);input;estim_flag(2);estim_flag(3)];
end

sigma = input(4)/sqrt(n);
k2 =2^kbar;

A = transition_mat(A_template,input,kbar);
g_m = gofm(input,kbar);       
T = length(data);                       
    
LLs = zeros(1,T);
pi_mat = zeros(T+1,k2); 
pi_mat(1,:) = (1/k2)*ones(1,k2);


%*----------------------------------------------------------------------*
%*                        Likelihood algorithm                          *
%*----------------------------------------------------------------------*
pa = (2*pi)^-0.5;
s = repmat(sigma*g_m,T,1);
w_t = repmat(data,1,k2);
w_t = pa*exp( - 0.5.*((w_t./s).^2))./s; 
w_t = w_t + 1e-16;

pmat = pi_mat;

%ARX_MSM_Core is a mex function. Please compile before you use. Else, uncomment and
%use the loop below. But the loop may be very slow for a very large data
%and big kbar.
[pi_mat, LLs]=ARX_MSM_Core(pmat,A,w_t);  
% for t=2:T+1          
%     piA = (pi_mat(t-1,:)*A);
%     C = (w_t(t-1,:).*piA); ft = sum(C);
%     if ft == 0                      %This stop div by zero if probs are too low
%         pi_mat(t,1) = 1;   
%     else
%         pi_mat(t,:) = C / ft; 
%     end
%     
%     LLs(t-1) = log(dot(w_t(t-1,:),piA));
%  end 
    LL=-sum(LLs);
    if ~isfinite(LL)
        disp('Log-likelihood is inf. Probably due to all zeros in pi_mat.')
    end
end



% Calculate the transition matrix using the template
function A = transition_mat(A,input,kbar)
    b = input(1);
    gamma_kbar = input(3);
    
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

% calculate all of the possible volatility states
function g_m = gofm(input,kbar)
    m0 = input(2); 
    m1=2-m0;
    kbar2 = 2^kbar;
    g_m1 = [0:(kbar2-1)];

    for i = 1:(kbar2)
        g=1;
        for j = 0:(kbar-1)       
            if(bitand(g_m1(i),(2^j))~=0)    %
                g=g*m1;
            else g=g*m0;
            end
        end
        g_m1(i)=g;
    end
    
    g_m=sqrt(g_m1);

end



