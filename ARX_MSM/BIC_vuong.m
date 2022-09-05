% This algorithm computes vuong statistics, with BIC correction (i.e corrected for number of parameters)
% and the coresponding t-value and p-value
%                                   Author:
%                                Waleem Alausa
%                            alausa.babs@gmail.com
%--------------------------------------------------------------------------
% Input:             l1= T-by-1 vector of expected log-likelihood for base model i.e model 1
%                    l2= T-by-1 vector of expected log-likelihood for comparison model i.e model 2
%                    H0: model1 and model2 are equivalent
%                    H1: model1 is worse than model2
%                    See Vuong(1989) pg.11. 
%                    l1-l2<0 means model 1 is worse than model2
%                    b1=number of parameters in model 1
%                    b2=number of parameters in model 2
%Output:             t=t-statistics
%                    p=lower-tail p-value
%--------------------------------------------------------------------------
%Interpreting results: Note that positive t-value automatically means the
%null CANNOT be rejected. For negative t-values, check the p-values along
%the columns. E.g comparing all other models against MSM(10), i.e H0: MSM(k) is equivalent to MSM(10). You should be
%looking at p-values on column 10. For 5% significance levels, any
%p-value<=0.05 means H0 is rejected i.e MSM(k) is worse than MSM(10) (or loosely speaking, MSM(10) is better than MSM(k) )
%Any p-value >0.05 means the null cannot be rejected.
function[t,p]=BIC_vuong(l1,l2,b1, b2)
  if size(l1,2)>1 %both l1 and l2 must each be column vectors
      l1=l1';
  end
  if size(l2,2)>1
      l2=l2';
  end
  if length(l1)~=length(l2)
      error('Vector of expected log-likelihood must be of the same length for both models')
  end
  T  = size(l1,1);
  l  = l1-l2;
  v1 = var(l);
  
  t=(sum(l) - ((b1-b2)/2)*log(T))/(sqrt(T)*v1); % One-tail p-value of Z-test
  p = tdis_prb(abs(t),T-1)/2;

%p=normcdf(-abs(t),0,1);
end
