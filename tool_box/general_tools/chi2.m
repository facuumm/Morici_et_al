n1 = 46; N1 = 278;
n2 = 12; N2 = 66;
n3 = 11; N3 = 142;
x1 = [repmat('a',N1,1); repmat('b',N2,1) ; repmat('c',N3,1)];
x2 = [repmat(1,n1,1); repmat(2,N1-n1,1); repmat(1,n2,1); repmat(2,N2-n2,1) ; repmat(1,n3,1); repmat(2,N3-n3,1)];
[tbl,chi2stat,pval] = crosstab(x1,x2)
