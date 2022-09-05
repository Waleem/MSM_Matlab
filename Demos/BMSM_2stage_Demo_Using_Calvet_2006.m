% NB: to replicate Calvet's results on table 7 (DM JA), don't multiply
% returns by 100 i.e r= lnP-lnP_1, only. The results of JA-UK are also
% replicated pefectly as follows [ 1.6311    1.5752    0.0070    0.0043 0.6949   13.5759   -0.4390    0.5255].
% The only difference with the results of Calvet 2006 is that they reported
% the estimate of rho_e as 0.439 instead of -0.439. Their estimate of rho_e
% is wrong because UK is negatively correlated with all other currencies.
% See simple correlation matrix below:
%   UKRET  1.0000   -0.4440   -0.4047   -0.1503
%   DMRET -0.4440    1.0000    0.3823    0.0683
%   JARET -0.4047    0.3823    1.0000    0.0783
%   CARET -0.1503    0.0683    0.0783    1.0000

%DM-UK = [1.4715    1.5737    0.0078    0.0043    0.6814   11.3663   -0.5859    0.8899]
%JA-UK = [ 1.6311    1.5752    0.0070    0.0043 0.6949   13.5759   -0.4390    0.5255]
%DM-JA = [ 1.4449    1.5784    0.0050    0.0048    0.8438    9.1229    0.5169    0.6723]
%----------------------------------------------------------------------------------------------------------------
clc;
clear;
load Calvet_2006_Data %[UK DM YEN CAN]


[out1,L1] = RES_BMSM_2stage([Return(:,2) Return(:,3)],5,252); %DM-JA
[out2,L2] = BMSM_2stage([Return(:,3) Return(:,1)],2,252); %JA-UK
% out2 =       para: [1.6311 1.5752 0.7021 0.4318 0.6949 13.5770 -0.4390 0.5256 0.9999]
%              se: [0.0126 0.0129 0.0570 0.0541 0.0704 2.9989 0.0107 0.0938 0.1345]

