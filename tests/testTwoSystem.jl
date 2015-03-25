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

#println(rhs);

#println("Factorize complex matrix")
#tic();
#F1 = factorPARDISO(A,1);
#toc();

println("Factorize real matrix")
tic();


println(A2.colptr);
println(A2.rowval);
println(A2.nzval);

pardiso = ParDiSO(-2, 1);

initPARDISO(pardiso);


#F2 = factorPARDISO(A2,2);
sym = 0;
ooc = 0;


println("DEBUGGING: \n\n\n");

toc();
exit()

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


