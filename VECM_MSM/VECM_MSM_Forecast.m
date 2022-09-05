function [covt,vars,varf,phat,std1,std2]=VECM_MSM_Forecast(input,k)
% This is a helper function used to compute variance and covariance forcast
% for VECM_M_BMSM.
% Input:
%      Input      -   A 1-by-1 output structure returned by VECM_M_MSM
%       K         -   The number of volatility components in VECM_M_MSM
% Output:         - 
%        covt    -   a T-by1 vector of conditional covariance
%        vars    -   a T-by1 vector of forecasted variance for spot series
%        varf    -   a T-by1 vector of forecasted variance for futures series
%        stds    -   a T-by1 vector of forecasted std. deviation for spot series
%        stdf    -   a T-by1 vector of forecasted std. deviation for futures series


para=input.para;
%Remove the VECM parameters. We need only the RES_BMSM parameters
% Remember to add rho_m. Even though it is assumed to be 1 in VECM_M_BMSM,
% BMSM_2stage_pimat expects it.
para1 = [para(1:6) para(11:12) 1];

[p, A2] = BMSM_2stage_pimat(para1,k,input.resid);
phat   = p(end,:)*A2;
[~,~,~,covt,vars,varf,std1,std2]=BMSM_2stage_comvt(phat,para1,k);

end