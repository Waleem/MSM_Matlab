function [LL, LLs] = get_ARX_MSM_likelihood(para, kbar, Y, xmat, n, model, p)

[LL, LLs] = ARX_MSM_likelihood(para, kbar, Y, xmat, n, model, p);
LL=-LL;

end