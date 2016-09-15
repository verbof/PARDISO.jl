# PARDISO
Julia drivers for PARDISO libraries.
This interface is based on the type **ParDiSO**. This type contains the information needed by PARDISO to perform the essential operations.

## Installation
Clone this repository in the julia packages folder:

    cd $JULIA_PKGDIR                                        (default: ~/.julia/v0.4)
    git clone https://github.com/verbof/PARDISO.jl PARDISO
    cd PARDISO/src
                    
Set correctly the `FC` variable to the correct value for your system in `Makefile`; then, modify the correct PARDISO shared library needed for compiling.


## The ParDiSO type
The basic constructor 

        ParDiSO(matrixtype::Integer, msglevel::Integer);

creates a variable that is going to be used for a matrix of the the type `matrixtype` and with a verbose-level `msglevel` (`0`: no messages from PARDISO, `1`: verbose mode).

## The SparsePardisoCSR type
Since PARDISO is based on a **CSR (Compressed Sparse Row)** representation for sparse matrices, the Julia's `SparseMatrixCSC{T<:Integer} <: AbstractSparseMatrix{T,Int32}` have to be converted into `SparsePardisoCSR` format. In order to do so, the constructors

        SparsePardisoCSR{Float64}(A::SparseMatrixCSC{Float64,Int})
        SparsePardisoCSR{Complex128}(A::SparseMatrixCSC{Complex128,Int})

are provided. **NOTE** that this class saves just the upper triangular part of the original matrix if it belongs to one of the following classes:
* Real Symmetric Positive Definite
* Real Symmetric Positive Undefinite
* Complex Hermitian Positive Definite
* Complex Hermitian Undefinite
* Complex Symmetric

## PARDISO initialization
Before performing any operation using PARDISO, it is needed to initialize the **ParDiSO** object using the `initPARDISO` function.

_*Important*: **never** change the values of `pt`_.


## PARDISO usage
The PARDISO solver relies on the matrix factorization to solve the systems. So, PARDISO needs to go through a preliminar *factorization* fase.

### Check Matrix (just for debugging purposes)
Checks if the input matrix is consistent according to the matrix-type.

        checkPARDISO{Tnzval}(pardiso::ParDiSO, A::SparsePardisoCSR{Tnzval})

### Symbolic factorization
Performs a reordering of the rows and a symbolic factorization on the matrix in order to reduce the fill-in in the factorization phase.

        smbfctPARDISO{Tnzval}(pardiso::ParDiSO, A::SparsePardisoCSR{Tnzval})

### Factorization
Computes an **LU** factorization for the matrix and stores internally the factors.

        factorPARDISO{Tnzval}(pardiso::ParDiSO, A::SparsePardisoCSR{Tnzval})

### Solve system
Solves the system **Ax = b** using the internal factors and performs an iterative refinement.

        solvePARDISO{Tnzval}(pardiso::ParDiSO, A::SparsePardisoCSR{Tnzval}, n_rhs::Int64, b::Vector{Tnzval})

### Free memory
Realeases internal PARDISO memory.

        freePARDISO(pardiso::ParDiSO)

# Example code

```
using PARDISO;

M = 1000;
f = rand(M-1) + im*rand(M-1);
A = speye(M) + spdiagm(f,+1, M,M) + spdiagm(f,-1, M,M);      # Creates a COMPLEX SYMMETRIC matrix

pardiso = ParDiSO(+6, 1);

initPARDISO(pardiso);

X = rand(Complex128, M, 5);         # Random solution (5 Right-Hand Sides)
b = A*X;                            # RHS's
X = colwise(X);                     # column-wise representation (see PARDISO.jl)

A = SparsePardisoCSR(A);

println("PARDISO CHECK MATRIX");
checkPARDISO(pardiso, A);


println("SYMBOLIC FACTORIZATION + FACTORIZATION");
smbfctPARDISO(pardiso, A);
factorPARDISO(pardiso, A);

println();
println("SOLVE COMPLEX SYSTEM");
@time x = solvePARDISO(pardiso, A, b);                          # Solution computed by PARDISO

println("Residual:     ", norm(X-x)/norm(X));                   # Relative error on the solution
println("Total memory: ", memoryPARDISO(pardiso), " kylobites.\n");

freePARDISO(pardiso);
```
