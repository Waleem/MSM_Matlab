function [LL, LLs] = get_MSM_likelihood(para, kbar,n, data,A_template)

[LL,LLs] = MSM_likelihood(para,kbar,n,data,A_template);

LL=-LL;

end