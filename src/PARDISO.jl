module PARDISO

importall Base;

    type ParDiSO 

        pt::Array{Int64};
        maxfct::Int64;
        mnum::Int64;
        mtype::Int64;
        phase::Int64;
        iparm::Array{Int32};
        solver::Int64;
        msglvl::Int64;
        error_::Int64;
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

    PARDISOLIB = "$(ENV["HOME"])/.julia/v0.3/PARDISO/lib/PADISO"

    include("pardisobase.jl")

    export  ParDiSO;
    export  initPARDISO,
            checkPARDISO,
            symfctPARDISO,
            factorPARDISO;

end
