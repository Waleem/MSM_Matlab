function T=ARX_MSM_Transition(gamma)
%Computes the 2-by-2 transition matrix of M_k for univariate MSM
%Input: 
%        gamma  - Scalar. The switching probability for volatility component, M_k,t

%Output: 
%        T      - 2-by-2 symmetric transition matrix with T_i,j = P(M_t+1 =m^j| M_t=m^i)


    T=[1-gamma + (.5*gamma)  .5*gamma;...
        .5*gamma             1-gamma + (.5*gamma) ];