module PARDISO

importall Base;

    type ParDiSO 

        pt::Array{Int64};
        maxfct::Int32;
        mnum::Int32;
        mtype::Int32;
        phase::Int32;
        iparm::Array{Int32};
        solver::Int32;
        msglvl::Int32;
        error_::Int32;
        dparm::Array{Float64};


        ParDiSO(matrixtype::Int64, msglevel::Int64) = 
        begin
            new(zeros(Int64,64), 
                1,
                1,
                matrixtype,
                0,
                zeros(Int32,64),
                0,
                msglevel,
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



    type SparseMatrixCSR{Tv,Ti<:Integer} <: AbstractSparseMatrix{Tv,Ti}
        m::Int              # Number of rows
        n::Int              # Number of columns
        colptr::Vector{Ti}  # Column i is in colptr[i]:(colptr[i+1]-1)
        rowval::Vector{Ti}  # Row values of nonzeros
        nzval::Vector{Tv}   # Nonzero values


        SparseMatrixCSR{Tv,Ti}(m::Integer, n::Integer, colptr::Vector{Ti}, rowval::Vector{Ti}, nzval::Vector{Tv}) =
                               SparseMatrixCSR(int(m), int(n), colptr, rowval, nzval)

        SparseMatrixCSR{Tv,Ti}(A::SparseMatrixCSC{Tv,Ti}) =
                               SparseMatrixCSR(A.m, A.n, (A').colptr, (A').rowval, A.nzval)

    end
    



    PARDISOLIB = "$(ENV["HOME"])/.julia/v0.3/PARDISO/lib/PADISO"

    include("pardisobase.jl")

    export  ParDiSO, SparseMatrixCSR, printPARDISO;
    export  initPARDISO,
            checkPARDISO,
            smbfctPARDISO,
            factorPARDISO,
            nothingPARDISO;

end
