clc;
clear all;
load caret
caret = caret*100;
[output,LL] = MSM(caret, 2,252);
%   b         m0       gamma_k   sigma
% 4.7564    1.4739    0.1288    0.2929
% 1.6105    0.0186    0.0755    0.0297
