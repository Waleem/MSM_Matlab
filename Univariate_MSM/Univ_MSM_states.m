function [Mmat, mvec,sqrt_mvec] = Univ_MSM_states(m0,kbar)

%Computes all 2^kbar possible volatility states for a univariate MSM, with kbar volatility components.
%It is also faster than the method of Calvet included in the MSM_likelihood function
%Input:   m01       -  high state for M_k^alpha (scalar)
%         kbar      -  number of volatility components (scalar)
%Output:  Mmat      -  2^kbar-by-kbar matrix of volatility components possible combinations
%         mvec      -  a vector of all possible volatility states. See Calvet & Fischer (2004), pg 8
%         sqrt_mvec -  a vector of sqrt(mvec)

m1=2-m0;
m=[m0 m1];

for n = 1:kbar
    eval(['P', int2str(n), '= m;']);
end
s=[]; 
for d=1:kbar
    if d<kbar
   s = [s 'P' int2str(d) ','];
    else s = [s 'P' int2str(d)];
    end
end
b=['allcomb(' s ')'];
Mmat=eval(b);
mvec=(prod(Mmat,2))';
sqrt_mvec=sqrt(mvec);

end