# Factorize two different matrices (one real, one complex) and solve for multiple rhs and free memory

using PARDISO
using Base.Test
include("getDivGrad.jl");


A  = getDivGrad(2,2,2);
A2 = getDivGrad(2,2,2);


n  = size(A,1);
n2 = size(A2,1);
A  = A + im*spdiagm(rand(n),0);


println("******* DEBUG *******");
println("PARDISO INITIALISATION");

pardiso = ParDiSO(+2, 1);

a = initPARDISO(pardiso);

#d = 0.5 * ones(9);
#A = speye(10) + spdiagm(d,+1, 10,10) + spdiagm(d,-1, 10,10);
A = speye(10) + sprand(10,10, 0.1); A = A + A';
detA = det(A);
#b    = [sum(A,2) 2*sum(A,2) 3*sum(A,2)];
#b = sum(A,2);

X = rand(10,5);     # Random solution

b = A*X;            # RHS's
X = colwise(X);     # column-wise representation

A = SparsePardisoCSR(A);

println("PARDISO CHECK MATRIX");

@time checkPARDISO(pardiso, A);

println("SYMBOLIC FACTORIZATION + FACTORIZATION");

@time smbfctPARDISO(pardiso, A);
@time factorPARDISO(pardiso, A);


if pardiso.iparm[33] == 1
    println("Determinant of the matrix:     ", exp(pardiso.dparm[33]));
    println("Residual between determinants: ", abs(exp(pardiso.dparm[33]) - detA));
end

println("SOLVE REAL SYSTEM");

@time x = solvePARDISO(pardiso, A, b);      # Solution computed by PARDISO

println("Residual: ", norm(X-x)/norm(X));   # Relative error on the solution

@test norm(X-x)/norm(X) < 1e-15;

println("*********************");
exit();

#println("Solve complex system")
#tic();
#x = applyPARDISO(F1,rhs);
#toc();
#err = zeros(nrhs)
#for i=1:nrhs
#        err[i] =  norm(A*x[:,i]-rhs[:,i]) / norm(rhs[:,i]);
#end
#@test maximum(err) < 1e-14

@test maximum(err) < 1e-14

println("Free memory")
destroyPARDISO(F1);
destroyPARDISO(F2);

println("DONE!")


