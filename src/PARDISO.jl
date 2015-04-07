module PARDISO

    importall Base;

    type ParDiSO 

        pt::Vector{Int64};
        maxfct::Int32;
        mnum::Int32;
        mtype::Int32;
        phase::Int32;
        iparm::Vector{Int32};
        solver::Int32;
        msglvl::Int32;
        error_::Int32;
        dparm::Vector{Float64};


        ParDiSO(matrixtype::Integer, msglevel::Integer) = 
        begin
            new(zeros(Int64,64), 
                1,
                1,
                int32(matrixtype),
                0,
                zeros(Int32,64),
                0,
                int32(msglevel),
                0,
                zeros(Float64,64)
                );

        end
    end
    
    show(io::IO, p::ParDiSO) = 
    begin 
        println(io, "ParDiSO(");
        println(io, "        maxfct: $(p.maxfct)");
        println(io, "        mnum:   $(p.mnum)");
        println(io, "        mtype:  $(p.mtype)");
        println(io, "        phase:  $(p.phase)");
        println(io, "        solver: $(p.solver)");
        println(io, "        msglvl: $(p.msglvl)");
        println(io, "        error_: $(p.error_)");
        println(io, "       );");
    end

    printPARDISO(p::ParDiSO) = 
    begin
        println("ParDiSO(");
        println("        maxfct: $(p.maxfct)");
        println("        mnum:   $(p.mnum)");
        println("        mtype:  $(p.mtype)");
        println("        phase:  $(p.phase)");
        println("        solver: $(p.solver)");
        println("        msglvl: $(p.msglvl)");
        println("        error_: $(p.error_)");
        println("        ********");
        println("        iparm:  $(p.iparm)");
        println("        dparm:  $(p.dparm)");
        println("       );");

        
    end


    #type PardisoMatrixNotSquareException <:Exception 
    #    A::Symbol
    #end
    #Base.showerror(io::IO, e::PardisoMatrixNotSquare) = print(io, e.A " not defined");

    type SparsePardisoCSR{T<:Number} <: AbstractSparseMatrix{T,Int32}

        upper::Bool;                # Matrix is stored in upper-triangular form
        n::Int32;                   # Number of rows/columns
        rowptr::Vector{Int32};      # i-th row is in rowptr[i]:(rowptr[i+1]-1)
        colval::Vector{Int32};      # Column values of nonzeros
        nzval::Vector{T};           # Nonzero values

    end 

        
    SparsePardisoCSR{T}(upper::Bool, rowptr::Vector{Int32}, colval::Vector{Int32}, nzval::Vector{T}) =
                        SparsePardisoCSR(upper, int32(length(rowptr)-1), rowptr, colval, nzval);

    SparsePardisoCSR{Float64}(A::SparseMatrixCSC{Float64,Int}) =
    begin
        if A.m != A.n
            error("Matrix must be square, but size = ($m, $n).\n");
            #throw(PardisoMatrixNotSquareException(A));
        else
            if issym(A)
                SparsePardisoCSR(true, 
                    int32(A.n), 
                    convert(Vector{Int32},tril(A).colptr), 
                    convert(Vector{Int32},tril(A).rowval), 
                    tril(A).nzval);
            else
                SparsePardisoCSR(false,
                    int32(A.n), 
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
                    int32(A.n), 
                    convert(Vector{Int32},tril(A).colptr), 
                    convert(Vector{Int32},tril(A).rowval), 
                    tril(A).nzval);
            
            elseif ishermitian(A)
                SparsePardisoCSR(true, 
                    int32(A.n), 
                    convert(Vector{Int32},tril(A).colptr), 
                    convert(Vector{Int32},tril(A).rowval), 
                    tril(A.').nzval);

            else
                SparsePardisoCSR(false,
                    int32(A.n), 
                    convert(Vector{Int32},(A.').colptr),
                    convert(Vector{Int32},(A.').rowval),
                    (A.').nzval);
            end
        end
    end 

    size(S::SparsePardisoCSR) = (S.n, S.n)
    nnz(S::SparsePardisoCSR)  = int(S.rowptr[end]-1)
    
    colwise(x::Array) = squeeze(reshape(x, prod(size(x)), 1), 2);

    #PARDISOLIB = "$(ENV["HOME"])/.julia/v0.3/PARDISO/lib/PADISO"

    include("pardisobase.jl")

    export  ParDiSO, SparsePardisoCSR, printPARDISO, colwise;
    export  initPARDISO,
            checkPARDISO,
            smbfctPARDISO,
            factorPARDISO,
            solvePARDISO,
            memoryPARDISO,
            freePARDISO;

end
