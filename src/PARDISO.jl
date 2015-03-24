module PARDISO

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


        ParDiSO(matrixtype::Int64, msglevel::Int64) = begin
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

                                                        if "OMP_NUM_THREADS" in keys(ENV)
                                                            iparm[3] = int32(ENV["OMP_NUM_THREADS"]);
                                                        else
                                                            error("Set environment variable OMP_NUM_THREADS...");
                                                        end
                                                      end 

    end


    PARDISOLIB = "$(ENV["HOME"])/.julia/v0.3/PARDISO/lib/PADISO"

    include("pardisobase.jl")

    export  ParDiSO;
    export  initPARDISO,
            checkPARDISO,
            symfctPARDISO,
            factorPARDISO;

end
