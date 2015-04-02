# Factorize two different matrices (one real, one complex) and solve for multiple rhs and free memory

using PARDISO
using Base.Test
include("getDivGrad.jl");


A  = getDivGrad(2,2,2);
A2 = getDivGrad(2,2,2);


n  = size(A,1);
n2 = size(A2,1);
A  = A + im*spdiagm(rand(n),0);

#println(A);

nrhs = 10;
rhs = randn(n,nrhs) + im*randn(n,nrhs);
rhs2 = randn(n2,nrhs);

println("******* DEBUG *******");
println("PARDISO INITIALISATION");
pardiso = ParDiSO(-2, 1);

a = initPARDISO(pardiso);

#println(pardiso.pt);

println("PARDISO CHECK MATRIX");

#d = rand(9);
d = 0.5 * ones(9);
B = speye(10) + spdiagm(d,+1, 10,10) + spdiagm(d,-1, 10,10);

detB = det(B);

B = convert(SparseMatrixCSC{Float64,Int32},B);


B = SparsePardisoCSR(B);

#println(B.rowptr);
#println(B.colval);

@time checkPARDISO(pardiso, B);

println("SYMBOLIC FACTORIZATION + FACTORIZATION");

@time smbfctPARDISO(pardiso, B);
@time factorPARDISO(pardiso, B);


if pardiso.iparm[33] == 1
    println("Determinant of the matrix:     ", exp(pardiso.dparm[33]));
    println("Residual between determinants: ", abs(exp(pardiso.dparm[33]) - detB));
end

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

println("Solve real system")
tic();
x2 = applyPARDISO(F2,rhs2);
toc();
err = zeros(nrhs)
for i=1:nrhs
        err[i] =  norm(A2*x2[:,i]-rhs2[:,i]) / norm(rhs2[:,i]);
end
@test maximum(err) < 1e-14

println("Free memory")
destroyPARDISO(F1);
destroyPARDISO(F2);

println("DONE!")


