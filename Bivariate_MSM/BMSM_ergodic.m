function ergo = BMSM_ergodic(gamma,lamda,rho_m,kbar)
%Computes the 1-by-4^kabr ergodic distribution for bivariate MSM. 

pihh=zeros(kbar,1);
pihl=pihh;
for k=1:kbar
    pihh(k)=(1/4)*((1-((1-rho_m)*((1-lamda)*gamma(k)+lamda))/2)/(1-((1-lamda)*gamma(k)+lamda)/2));
    pihl(k)=(1/2)-pihh(k);
end
    pill=pihh;
    pilh=pihl;
    pik=[pihh pihl pilh pill];
         
ergo=pik(1,:);
if kbar>1
    for i=2:kbar
        ergo=kron(pik(i,:),ergo);
    end
end

end