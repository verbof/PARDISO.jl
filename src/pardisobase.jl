
# export errorPARDISO,
#        initPARDISO,
#        checkPARDISO,
#        smbfctPARDISO,
#        factorPARDISO,
#        solvePARDISO



function errorPARDISO(error_)

    if error_ != 0
  
        print("ParDiSO error: $(error_) --- ");

        if error_ == -1
            print("input inconsistent\n");
        elseif error_ == -2
            print("not enough memory\n");
        elseif error_ == -3
            print("reordering problem\n");
        elseif error_ == -4
            print("zero pivot, numerical factorization or iterative refinement problem\n");
        elseif error_ == -5
            print("unclassified (internal) error\n");
        elseif error_ == -6
            print("preordering failed (matrix type 11, 13 only)\n");
        elseif error_ == -7
            print("diagonal matrix problem\n");
        elseif error_ == -8
            print("32 bit integer overflow problem\n");
        elseif error_ == -10
            print("no license file pardiso.lic found\n");
        elseif error_ == -11
            print("license is expired.\n");
        elseif error_ == -12
            print("wrong username or hostname.\n");
        end

    else
        println();
    end

end


function nothingPARDISO()

    ccall((:pardiso_nothing, "/users/verbof/.julia/v0.3/PARDISO/lib/PARDISO"), Void, ());

end



function initPARDISO(pardiso::ParDiSO)

    ccall( (:pardiso_init_, "/users/verbof/.julia/v0.3/PARDISO/lib/PARDISO"), Void,
           (Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int32}, Ptr{Float64},  Ptr{Int64}),
            pardiso.pt, &pardiso.mtype, &pardiso.solver, pardiso.iparm, pardiso.dparm, &pardiso.error_);

    errorPARDISO(pardiso.error_);

    if "OMP_NUM_THREADS" in keys(ENV)
        pardiso.iparm[3] = int32(ENV["OMP_NUM_THREADS"]);
    else
        error("Set environment variable OMP_NUM_THREADS before running this code.");
    end

    
end



function checkPARDISO(pardiso::ParDiSO, A::SparseMatrixCSC{Float64,Int32})

    if abs(pardiso.mtype) == 2

        if issym(A)
    
            M = tril(A);

            ccall( (:pardiso_checkmatrix_, "/users/verbof/.julia/v0.3/PARDISO/lib/PARDISO"), Void,
                    ( Ptr{Int64},   Ptr{Int64}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int64}),
                    &pardiso.mtype, &(M.n),     M.nzval,      M.colptr,   M.rowval,   &pardiso.error_);

        else
            error("pardiso.mtype = ±2, but the matrix is not symmetric.");
        end

    elseif pardiso.mtype == 11 || pardiso.mtype == 13

        ccall( (:pardiso_checkmatrix_, "/users/verbof/.julia/v0.3/PARDISO/lib/PARDISO"), Void,
                    ( Ptr{Int64},    Ptr{Int64}, Ptr{Float64},  Ptr{Int32}, Ptr{Int32}, Ptr{Int64}),
                    pardiso.mtype,   &(A.n),     A.nzval,       A.colptr,   A.rowval,   pardiso.error_);

    end
        

    if pardiso.error_ != 0
        error("Error in consistency of matrix: $(pardiso.error_)");
    else
        println(" -> Matrix checked by PARDISO!");
    end

    errorPARDISO(pardiso.error_);

end


function smbfctPARDISO(pardiso::ParDiSO, A::SparseMatrixCSC{Float64,Int32})
	# Reordering and Symbolic Factorisation + memory allocation for factors

    pardiso.phase = 11;     # just analysis
    ddum   = 0.0;           # double dummy
    idum   = 0;             # integer dummy
    idum32 = int32(0);      # 32-bit integer dummy
    nrhs   = 1;

    pardiso.iparm[33] = 1;   # compute determinant of the matrix

    if A.n != A.m
	error("factorPARDISO: Matrix must be square!");
   	end

    println("Analysis of real matrix");

    M = tril(A);

    printPARDISO(pardiso);
    println();


    ccall( (:pardiso_call_, "/users/verbof/.julia/v0.3/PARDISO/lib/PARDISO"),
        Void,
        (Ptr{Int64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
        Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
        Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64},
        Ptr{Int32}, Ptr{Float64}),
        pardiso.pt, &pardiso.maxfct, &pardiso.mnum, &pardiso.mtype, &pardiso.phase, &(M.n),
        M.nzval, M.colptr, M.rowval, &idum32,
        &idum, pardiso.iparm, &pardiso.msglvl, &ddum, &ddum,
        &pardiso.error_, pardiso.dparm);
    
    if error_ != 0
        error("Error in symbolic factorisation of matrix: $(pardiso.error_)");
    end

    errorPARDISO(pardiso.error_);

    println("Symbolic factorisation is complete.");

    if iparm[33] != 0
        println("Logarithm of the determinant = $(dparm[33])");
    end

end



function factorPARDISO(pardiso::ParDiSO, A::SparseMatrixCSC{Float64,Int32})
	# Generate LU-factorization of a real matrix A

    pardiso.phase = 22;     # just numerical factorisation
    ddum  = 0.0;            # double dummy
    idum  = 0;              # integer dummy
    idum32 = int32(0);      # 32-bit integer dummy

    if A.n != A.m
		error("factorPARDISO: Matrix must be square!")
   	end

    println("Factorize real matrix")
    

    M = tril(A);

    println(M);
    println(full(M));
    println();

    ccall( (:pardiso_call_, "/users/verbof/.julia/v0.3/PARDISO/lib/PARDISO"),
        Void,
        (Ptr{Int64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
        Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
        Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64},
        Ptr{Int32}, Ptr{Float64}),
        pardiso.pt, &pardiso.maxfct, &pardiso.mnum, &pardiso.mtype, &pardiso.phase, &(M.n),
        M.nzval, M.colptr, M.rowval, &idum32,
        &idum, pardiso.iparm, &pardiso.msglvl, &ddum, &ddum,
        &pardiso.error_, pardiso.dparm);
#=
    ccall( (:pardiso_call_, "/users/verbof/.julia/v0.3/PARDISO/lib/PARDISO"),
        Void,
        (Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64},
        Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
        Ptr{Int64}, Ptr{Int32}, Ptr{Int64}, Ptr{Float64}, Ptr{Float64},
        Ptr{Int64}, Ptr{Float64}),
        pardiso.pt, &pardiso.maxfct, &pardiso.mnum, &pardiso.mtype, &pardiso.phase, &(M.n),
        M.nzval, M.rowval, M.colptr, &idum32,
        &idum, pardiso.iparm, &pardiso.msglvl, &ddum, &ddum,
        &pardiso.error_, pardiso.dparm);
=#



    errorPARDISO(pardiso.error_);

end



function solvePARDISO(pardiso::ParDiSO, A::SparseMatrixCSC{Float64,Int32}, b::Vector{Float64})
	# Generate LU-factorization of a real matrix A

    pardiso.iparm[12] = 1;  # solve A^T * x = b;
    pardiso.phase = 33;     # just numerical factorisation
    idum  = 0;              # integer dummy

    if A.n != A.m
		error("factorPARDISO: Matrix must be square!")
   	end

    println("Factorize real matrix")

    ccall( (:pardiso_call, "/users/verbof/.julia/v0.3/PARDISO/lib/PARDISO"),
        Void,
        (Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64},
        Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
        Ptr{Int64}, Ptr{Int32}, Ptr{Int64}, Ptr{Float64}, Ptr{Float64},
        Ptr{Int64}, Ptr{Float64}),
        pardiso.pt, &paridso.maxfct, &pardiso.mnum, &pardiso.mtype, &pardiso.phase, &(A.n),
        A.nzval, A.colptr, A.rowval, &idum,
        &idum, pardiso.iparm, &pardiso.msglvl, b, x,
        &pardiso.error_, pardiso.dparm);


    errorPARDISO(pardiso.error_);

    return x

    #factor = PARDISOfactorization(p,myid(),n,true,facTime);
	# finalizer(factor,destroyPARDISO)
    #return factor

end
