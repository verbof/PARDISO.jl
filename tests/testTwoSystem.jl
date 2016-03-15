# Factorize two different matrices (one real, one complex) and solve for multiple rhs and free memory

using PARDISO
include("getDivGrad.jl");

A = getDivGrad(8,8,8);
N = size(A,1);

println("********************** Testing REAL drivers **********************");
println("\n\n");

pardisoR = ParDiSo(-2, 1);

initPARDISO(pardisoR);

detA = det(A);

X = rand(N,5);      # Random solution
b = A*X;            # RHS's
X = colwise(X);     # column-wise representation

A = SparsePardisoCSR(A);

println("PARDISO CHECK MATRIX");
@time checkPARDISO(pardisoR, A);


println("SYMBOLIC FACTORIZATION + FACTORIZATION");
#pardisoR.iparm[33] = 1;

@time smbfctPARDISO(pardisoR, A);
@time factorPARDISO(pardisoR, A);

println();
if pardisoR.iparm[33] == 1
    println("Determinant of the matrix:     ", exp(pardisoR.dparm[33]));
    println("Absolute error on determinants: ", abs(exp(pardisoR.dparm[33]) - detA)/abs(detA));
end

println();
println("SOLVE REAL SYSTEM");
@time x = solvePARDISO(pardisoR, A, b);      # Solution computed by PARDISO

println("Residual:     ", norm(X-x)/norm(X));   # Relative error on the solution
println("Total memory: ", memoryPARDISO(pardisoR), " kylobites.\n");

freePARDISO(pardisoR);


println("\n\n");
println("********************** Testing COMPLEX drivers **********************");
println("\n\n");


B = getDivGrad(8,8,8);
M = size(B,1);
f = rand(M-1) + im*rand(M-1);
B = B + spdiagm(f,+1, M,M) + spdiagm(conj(f),-1, M,M);

pardisoC = ParDiSo(-4, 1);

initPARDISO(pardisoC);

Y = rand(Complex128, M, 10);        # Random solution
c = B*Y;                            # RHS's
Y = colwise(Y);                     # column-wise representation

#println(ishermitian(B) ? "HERMITIAN" : "NON-HERMITIAN");

B = SparsePardisoCSR(B);

println("PARDISO CHECK MATRIX");
@time checkPARDISO(pardisoC, B);


println("SYMBOLIC FACTORIZATION + FACTORIZATION");

@time smbfctPARDISO(pardisoC, B);
@time factorPARDISO(pardisoC, B);

println();

println("SOLVE COMPLEX SYSTEM");
@time y = solvePARDISO(pardisoC, B, c);         # Solution computed by PARDISO

println("Residual:     ", norm(Y-y)/norm(Y));   # Relative error on the solution
println("Total memory: ", memoryPARDISO(pardisoC), " kylobites.\n");

freePARDISO(pardisoC);

println("******* TESTS DONE! *******");
println("\n\n");
println("\n\n");
