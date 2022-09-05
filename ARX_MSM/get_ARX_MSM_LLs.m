function LLs = get_ARX_MSM_LLs(para, kbar, Y, xmat, n, model, p)

[~, LLs] = ARX_MSM_likelihood(para, kbar, Y, xmat, n, model, p);


end