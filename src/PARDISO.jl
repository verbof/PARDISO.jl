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
    
    type SparsePardisoCSR

        upper::Bool;            # Matrix is stored in upper-triangular form
        n::Int32;               # Number of rows/columns
        rowptr::Vector{Int32};  # i-th row is in rowptr[i]:(rowptr[i+1]-1)
        colval::Vector{Int32};  # Column values of nonzeros
        nzval::Vector{Float64}; # Nonzero values
        
        SparsePardisoCSR(upper::Bool, rowptr::Vector{Int32}, colval::Vector{Int32}, nzval::Vector{Float64}) =
                        new(upper, int32(length(rowptr)-1), rowptr, colval, nzval);

        SparsePardisoCSR(A::SparseMatrixCSC{Float64,Int32}) =
        begin
            if A.m != A.n
                error("Matrix must be square, but size = ($m, $n).\n");
                #throw(PardisoMatrixNotSquareException(A));
            else
                if issym(A)
                    new(true, int32(A.n), tril(A).colptr, tril(A).rowval, tril(A).nzval);
                else
                    new(true, int32(A.n), (A').colptr,    (A').rowval,    A.nzval);
                end
            end
        end 

    end


    size(S::SparsePardisoCSR) = (S.n, S.n)
    nnz(S::SparsePardisoCSR)  = int(S.rowptr[end]-1)
    

    #PARDISOLIB = "$(ENV["HOME"])/.julia/v0.3/PARDISO/lib/PADISO"

    include("pardisobase.jl")

    export  ParDiSO, SparsePardisoCSR, printPARDISO;
    export  initPARDISO,
            checkPARDISO,
            smbfctPARDISO,
            factorPARDISO,
            nothingPARDISO;

end
