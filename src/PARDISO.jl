module PARDISO

    importall Base;

    type ParDiSO 

        pt::Vector{Int64};      # Pointer to PARDISO memory
        maxfct::Int32;          # Maximunm number of factors w/ same non-zero structure
        mnum::Int32;            # Index of the matrix, 1 <= mnum <= maxfct.
        mtype::Int32;           # Matrix Type
        phase::Int32;           # PARDISO Phase
        iparm::Vector{Int32};   # Integer parameters
        solver::Int32;          # Solver. 0: sparse direct; 1: multi-recursive iterative
        msglvl::Int32;          # Message Level. 0: no messages; 1: output statistics
        error_::Int32;          # Error. On output indicates the error detected
        dparm::Vector{Float64}; # Float parameters


        ParDiSO(matrixtype::Integer = 0, msglevel::Integer = 0) = 
        begin
            new(zeros(Int64,64), 
                1,
                1,
                Int32(matrixtype),
                0,
                zeros(Int32,64),
                0,
                Int32(msglevel),
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

    colwise(x::Array) = squeeze(reshape(x, prod(size(x)), 1), 2);
                        # Array{T,2} ---> Array{T, 1}

    #PARDISOLIB = "$(ENV["HOME"])/.julia/v0.3/PARDISO/lib/PADISO"

    include("sparsepardiso.jl")
    include("pardisobase.jl")
    include("pardisoext.jl");

    export  ParDiSO, SparsePardisoCSR, printPARDISO, colwise;
    export  initPARDISO,
            checkPARDISO,
            smbfctPARDISO,
            factorPARDISO,
            solvePARDISO,
            memoryPARDISO,
            freePARDISO,
            pardisoSolveSystem;

end
