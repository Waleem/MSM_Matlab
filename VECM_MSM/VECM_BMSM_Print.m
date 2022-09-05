% This algorithm prints bivariate MSM results
%                                   Author:
%                                Waleem Alausa
%                            alausa.babs@gmail.com
%--------------------------------------------------------------------------

function b = VECM_BMSM_Print(LL,para,se)
lnl=repmat(LL,12,1);

%VECM_BMSM parameters are in the following order:
%                                  [m01, m02, sigma1,sigma2, gamma_k, b, bs, bf, as, af, rho_e, lambda]
%Rearrange as                      [bs, bf, as, af, m01, m02, sigma1,sigma2, gamma_k, b,  rho_e, lambda]
para2 = [para(7:10,:); para(1:6,:); para(11:12,:)];
se2   = [se(7:10,:);   se(1:6,:);   se(11:12,:)];
xx=['$beta_s$  ';'          ';'          ';'          ';...
    '$beta_f$  ';'          ';'          ';'          ';...
    '$alpha_s$ ';'          ';'          ';'          ';...
    '$alpha_f$ ';'          ';'          ';'          ';...
    '$m_s$     ';'          ';'          ';'          ';...
    '$m_f$     ';'          ';'          ';'          ';...
    '$\sigma_s$';'          ';'          ';'          ';...
    '$\sigma_f$';'          ';'          ';'          ';...
    '$\gamma$  ';'          ';'          ';'          ';...
    '$b$       ';'          ';'          ';'          ';...
    '$\rho_e$  ';'          ';'          ';'          ';...
    '$\lambda$ ';'          ';'$lnL$     '];
info2.ldum = 1;
a=mprint3(para2,se2,lnl,info2);

b=[xx a];
b(3:4,:)=[];
b(5:6,:)=[];
b(7:8,:)=[];
b(9:10,:)=[];
b(11:12,:)=[];
b(13:14,:)=[];
b(15:16,:)=[];
b(17:18,:)=[];
b(19:20,:)=[];
b(21:22,:)=[];
b(23:24,:)=[];

end
