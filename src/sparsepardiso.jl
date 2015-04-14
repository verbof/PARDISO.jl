#######
# Compressed Sparse Row (CSR) matrices needed by PARDISO
#######

type SparsePardisoCSR{T<:Number} <: AbstractSparseMatrix{T,Int32}

    upper::Bool;                # Matrix is stored in upper-triangular form
    n::Int32;                   # Number of rows/columns
    rowptr::Vector{Int32};      # i-th row is in rowptr[i]:(rowptr[i+1]-1)
    colval::Vector{Int32};      # Column values of nonzeros
    nzval::Vector{T};           # Nonzero values

end 

    
SparsePardisoCSR{T}(upper::Bool, rowptr::Vector{Int32}, colval::Vector{Int32}, nzval::Vector{T}) =
                    SparsePardisoCSR(upper, Int32(length(rowptr)-1), rowptr, colval, nzval);

SparsePardisoCSR{Float64}(A::SparseMatrixCSC{Float64,Int}) =
begin
    if A.m != A.n
        error("Matrix must be square, but size = ($m, $n).\n");
        #throw(PardisoMatrixNotSquareException(A));
    else
        if issym(A)
            SparsePardisoCSR(true, 
                Int32(A.n), 
                convert(Vector{Int32},tril(A).colptr), 
                convert(Vector{Int32},tril(A).rowval), 
                tril(A).nzval);
        else
            SparsePardisoCSR(false,
                Int32(A.n), 
                convert(Vector{Int32},(A').colptr),
                convert(Vector{Int32},(A').rowval),
                A.nzval);
        end
    end
end 

SparsePardisoCSR{Complex128}(A::SparseMatrixCSC{Complex128,Int}) =
begin
    if A.m != A.n
        error("Matrix must be square, but size = ($m, $n).\n");
        #throw(PardisoMatrixNotSquareException(A));
    else
        if issym(A)
            SparsePardisoCSR(true, 
                Int32(A.n), 
                convert(Vector{Int32},tril(A).colptr), 
                convert(Vector{Int32},tril(A).rowval), 
                tril(A).nzval);
        
        elseif ishermitian(A)
            SparsePardisoCSR(true, 
                Int32(A.n), 
                convert(Vector{Int32},tril(A).colptr), 
                convert(Vector{Int32},tril(A).rowval), 
                tril(A.').nzval);

        else
            SparsePardisoCSR(false,
                Int32(A.n), 
                convert(Vector{Int32},(A.').colptr),
                convert(Vector{Int32},(A.').rowval),
                (A.').nzval);
        end
    end
end 

size(S::SparsePardisoCSR)     = (S.n, S.n)
nnz(S::SparsePardisoCSR)      = int(S.rowptr[end]-1)
countnz(S::SparsePardisoCSR)  = countnz(S.nzval)
density(S::SparsePardisoCSR)  = countnz(S)/prod(size(S));
nonzeros(S::SparsePardisoCSR) = S.nzval


function Base.showarray(io::IO, S::SparsePardisoCSR;
                        header::Bool=true, limit::Bool=Base._limit_output,
                        rows = Base.tty_size()[1], repr=false)

    if header
        print(io, S.n, "x", S.n, " sparse CSR matrix with ", nnz(S), " ", eltype(S), " entries:")
    end

    if limit
        half_screen_rows = div(rows - 8, 2)
    else
        half_screen_rows = typemax(Int)
    end
    pad = ndigits(S.n)
    k = 0
    sep = "\n\t"
    for row = 1:S.n, k = S.rowptr[row] : (S.rowptr[row+1]-1)
        if k < half_screen_rows || k > nnz(S)-half_screen_rows
            print(io, sep, '[', rpad(row, pad), ", ", lpad(S.colval[k], pad), "]  =  ")
            showcompact(io, S.nzval[k])
        elseif k == half_screen_rows
            print(io, sep, '\u22ee')
        end
        k += 1
    end

end

## Constructors

copy(S::SparsePardisoCSR) = 
    SparsePardisoCSR(S.upper, S.n, copy(S.rowptr), copy(S.colval), copy(S.nzval));


function full{T}(S::SparsePardisoCSR{T})
    A = zeros(T, S.n, S.n)
    for row = 1 : S.n, k = S.rowptr[row] : (S.rowptr[row+1]-1)
        A[row, S.colval[k]] = S.nzval[k]
    end
    return A
end


float(S::SparsePardisoCSR)    = SparsePardisoCSR(S.upper, S.n, copy(S.rowptr), copy(S.colval), float(copy(S.nzval)))
complex(S::SparsePardisoCSR)  = SparsePardisoCSR(S.upper, S.n, copy(S.rowptr), copy(S.colval), complex(copy(S.nzval)))
#complex(A::SparsePardisoCSR, 
#        B::SparsePardisoCSR)  = A + im*B



export density;
