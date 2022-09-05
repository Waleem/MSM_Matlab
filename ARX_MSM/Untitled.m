m=1;
k=8;
[output(k,m).para, LL(m,k), output(k,m).se, ~, ~, summary(m,k)] = ARX_MSM(lnp, p, k, m, X, n, 2);