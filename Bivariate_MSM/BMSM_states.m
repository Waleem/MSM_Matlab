function [gm,vol1,vol2] = BMSM_states(m01,m02,kbar)
%Computes all 4^kbar possible volatility states for a bivariate MSM, with
%kbar volatility components.
%Input:   m01  -  high state for M_k^alpha (scalar)
%         m02  -  high state for M_k^beta  (scalar)
%         kbar -  number of volatility components (scalar)
%Output:  gm   -  All 4^kar-by-2 possible volatility states
%         vol1 -  4^kar-by-kbar M_k values for M_alpha
%         vol1 -  4^kar-by-kbar M_k values for M_beta

%--------------------------------------------------------------------------------------------------------------
%This code was tested with Calvet_2006 data (DM-JA, JA-UK) and was found to
%be accurate. Results were replicated to 3 decimal places for all
%parameters. See further notes in calvet_2006.m file.
%---------------------------------------------------------------------------------------------------------------

m11=2-m01;
m12=2-m02;
st=[1 2 3 4];
s1=[m01 m02];
s2=[m01 m12];
s3=[m11 m02];
s4=[m11 m12];

for n = 1:kbar
    eval(['P', int2str(n), '= st;']);
end
s=[]; 
for d=1:kbar
    if d<kbar
   s = [s 'P' int2str(d) ','];
    else s = [s 'P' int2str(d)];
    end
end

B=zeros(4^kbar,kbar);
b=['B=allcomb(' s ');'];
evalc(b);

for n = 1:kbar
    eval(['M', int2str(n), '= zeros(length(B),2);']);
end

for j=1:kbar
    b=B(:,j);
    S=repmat(s1,length(B),1);
    eval(['M', int2str(j),'(b==1,:)=S(b==1,:);']);
    S=repmat(s2,length(B),1);
    eval(['M', int2str(j),'(b==2,:)=S(b==1,:);']);
    S=repmat(s3,length(B),1);
    eval(['M', int2str(j),'(b==3,:)=S(b==1,:);']);
    S=repmat(s4,length(B),1);
    eval(['M', int2str(j),'(b==4,:)=S(b==1,:);']);
end
clearvars S B b


vol1=M1(:,1); 
vol2=M1(:,2);
if kbar>1
for k=2:kbar
     vol1=eval(['cat(2,vol1,' 'M', int2str(k),'(:,1))']);
     vol2=eval(['cat(2,vol2,' 'M', int2str(k),'(:,2))']);
end
end

g_m11=(prod(vol1,2))';
g_m12=(prod(vol2,2))';
gm=[sqrt(g_m11);sqrt(g_m12)];

end


