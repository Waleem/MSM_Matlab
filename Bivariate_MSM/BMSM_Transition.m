function T=BMSM_Transition(gamma,lamda,rhom)
%Computes the 4-by-4 transition matrix of M_k for 2-stage bivariate MSM
%Input: 
%        gamma  - Scalar. The switching probability for volatility component, M_k,t
%        lambda - the correlation of arrival between M^alpha and M^beta
%        rhom   - the correlation between M^alpha and M^beta
%Output: 
%        T      - 4-by-4 transition matrix with t_i,j = P(M_t+1 =m^j| M_t=m^i)

p=1-gamma+(gamma*((1-lamda)*gamma+lamda))*((1+rhom)/4);
    q=1-gamma+(gamma*((1-lamda)*gamma+lamda))*((1-rhom)/4);
    T=[p 1-(gamma/2)-p 1-(gamma/2)-p gamma-1+p;...
        1-(gamma/2)-q q gamma-1+q 1-(gamma/2)-q;...
        1-(gamma/2)-q gamma-1+q q 1-(gamma/2)-q;...
        gamma-1+p 1-(gamma/2)-p 1-(gamma/2)-p p];