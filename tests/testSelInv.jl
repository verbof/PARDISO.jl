using PARDISO

println("****************   COMPLEX SYMM. SEL. INV.   ********************");

p = ParDiSo(6,1);
initPARDISO(p);

N = 50;
d = rand(Complex128, N-2);
T = triu(sprand(N,N, 0.5));

R = N*speye(Complex128, N) + spdiagm(d, -2, N,N) + spdiagm(d, +2, N,N) + T + T.';

println(issym(R) ? "* R is SYMM." : "* R is NOT SYMM.");
println(issym(SparseMatrixCSC(N,N, R.colptr, R.rowval, ones(length(R.nzval)))) ? "* R is STRUCT. SYMM. \n" : "* R is NOT STRUCT. SYMM. \n");

RR = SparsePardisoCSR(R);

 checkPARDISO(p, RR);
smbfctPARDISO(p, RR);
factorPARDISO(p, RR);

invRR = invertPARDISO(p, RR);

freePARDISO(p);

