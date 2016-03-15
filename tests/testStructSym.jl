using PARDISO

println("********************** REAL STRUCT. SYMM.  **********************");

p = ParDiSo(1,0);
initPARDISO(p);

N = 10;
d = rand(N-2);
T = triu(sprand(N,N, 0.5));

R = N*speye(N) + spdiagm(d-1, -2, N,N) + spdiagm(d+1, +2, N,N) + T - T'

println(issym(R) ? "* R is SYMM." : "* R is NOT SYMM.");
println(issym(SparseMatrixCSC(N,N, R.colptr, R.rowval, ones(length(R.nzval)))) ? "* R is STRUCT. SYMM. \n" : "* R is NOT STRUCT. SYMM. \n");

RR = SparsePardisoCSR(R);

checkPARDISO(p, RR);


smbfctPARDISO(p, RR);
factorPARDISO(p, RR);

x = rand(N);
b = R*x;

X = solvePARDISO(p, RR, b);
println("Relative error on the solution: ", norm(X-x)/norm(x));

freePARDISO(p);

println("********************** COMPLEX STRUCT. SYMM.  **********************");

p = ParDiSo(3,1);
initPARDISO(p);

M = 100;
d = rand(Complex128, M-2);
V = triu(sprand(M,M, 0.5));

C = M*speye(Complex128, M) + spdiagm(d-1, -2, M,M) + spdiagm(d+1, +2, M,M) + V - im*V.';

println(issym(C) ? "* C is SYMM." : "* C is NOT SYMM.");
println(issym(SparseMatrixCSC(M,M, C.colptr, C.rowval, ones(length(C.nzval)))) ? "* C is STRUCT. SYMM. \n" : "* C is NOT STRUCT. SYMM. \n");

CC = SparsePardisoCSR(C);

checkPARDISO(p, CC);


smbfctPARDISO(p, CC);
factorPARDISO(p, CC);


y = rand(M);
d = C*y;

Y = solvePARDISO(p, CC, d);

println("Relative error on the solution: ", norm(Y-y)/norm(y));


freePARDISO(p);

