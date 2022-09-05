clc
clear all
close all

load Aeso_forecast

% In-sample period: Jan/01/2011 - Dec/31/2011 => 8448 - 16918
% Forecast period : Jan/01/2012 - Jan/30/2012 => 16919- 17638 i.e 720 observations
in_data = data(1:1000,:);

P = in_data(:,3);
P(P<0.01) = 0.01;
lnp    = log(P);
lnload = log(in_data(:,6));
X=[lnload in_data(:,[8,10])];

tout = 100;
T    = size(lnp,1);
T_in = T-tout;
n = 366*24;  p=[1 2 24]; nw=2;

%matlabpool(8)
kbar  = 5;
model = 3;

output = struct('para',cell(1));
output = repmat(output,kbar,1);
yhat   = zeros(T,kbar);


parfor t=T_in:T-1
        [yhat(t+1,1), para(:,t)] = ARX_MSM_Fcast(lnp(1:t), p, 1, model, X(1:t+1,:), n, nw,[])
end
output(1).para=para;
startingvals = para;



for k=2:kbar
    parfor t=T_in:T-1
        [yhat(t+1,k), para(:,t)] = ARX_MSM_Fcast(lnp(1:t), p, k, model, X(1:t+1,:), n, nw,startingvals(:,t))
    end
    output(k).para=para;
    startingvals = para;
 
    
end

y=lnp(901:1000,:);
yhat=yhat(901:1000,:);
mean(abs(repmat(y,1,5)-yhat))

matlabpool close


