clc;
clear;

kbar=10;
b=2.7;
gamma_kbar=.959;
tic
A2 = zeros((2^kbar),(2^kbar));
    for i =0:2^kbar-1       
        for j = i:(2^kbar-1)-i  
            A2(i+1,j+1) = bitxor(i,j);
        end
    end
tic  
[A,g] = A_test(A2,b,gamma_kbar,kbar);
toc

tic
gamma = zeros(kbar,1);                          
gamma(1) = 1-(1-gamma_kbar)^(1/(b^(kbar-1)));
for i = 2:(kbar)
    gamma(i,1) = 1-(1-gamma(1,1))^(b^(i-1));
end
 %Compute transition matrices and gather them into opset array. There will be kbar transition matrices.
 opset=arrayfun(@(x) ARX_MSM_Transition(x),gamma,'UniformOutput', false);
 %Get kronecker product of all the transition matrices.
T = KronProd(opset,(kbar:-1:1),2*ones(kbar,1),1);
toc
% Elapsed time is 1.550789 seconds.
% Elapsed time is 0.016491 seconds. KronProd is about 100 times faster
%%
clc;
clear;
m0=1.5;
kbar=10;
tic
g_m = g_test(m0,kbar);
toc
tic
[Mmat, mvec,sqrt_mvec] = Univ_MSM_states(m0,kbar);
toc
% Elapsed time is 0.026769 seconds.
% Elapsed time is 0.006139 seconds.KronProd is about 5 times faster

%  min(sqrt_mvec==g_m)
% 
% ans =
% 
%      1
% 
% max(sqrt_mvec==g_m)
% 
% ans =
% 
%      1
