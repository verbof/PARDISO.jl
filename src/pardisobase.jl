
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


function initPARDISO(pardiso::ParDiSO)


    #const PARDISOLIB = "/users/verbof/.julia/v0.3/PARDISO/lib/PARDISO";


    ccall( (:pardiso_init_, "/users/verbof/.julia/v0.3/PARDISO/lib/PARDISO"), Void,
           (Ptr{Int64}, Ptr{Int32},     Ptr{Int32},      Ptr{Int32},    Ptr{Float64},  Ptr{Int32}),
            pardiso.pt, &pardiso.mtype, &pardiso.solver, pardiso.iparm, pardiso.dparm, &pardiso.error_);

    errorPARDISO(pardiso.error_);

    if "OMP_NUM_THREADS" in keys(ENV)
        pardiso.iparm[3] = int32(ENV["OMP_NUM_THREADS"]);
    else
        error("Set environment variable OMP_NUM_THREADS before running this code.");
    end

    
end



function checkPARDISO(pardiso::ParDiSO, A::SparsePardisoCSR)

    if abs(pardiso.mtype) == 2
        if !A.upper
            error("pardiso.mtype = ±2 but the matrix is not in upper triangular form (A.upper = false)");
        end
    end

    ccall( (:pardiso_checkmatrix_, "/users/verbof/.julia/v0.3/PARDISO/lib/PARDISO"), Void,
            ( Ptr{Int32},   Ptr{Int32}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}),
            &pardiso.mtype, &(A.n),     A.nzval,      A.rowptr,   A.colval,   &pardiso.error_);

    if pardiso.error_ != 0
        error("Error in consistency of matrix: $(pardiso.error_)");
    else
        println(" -> Matrix checked by PARDISO!");
    end

    errorPARDISO(pardiso.error_);

end


function smbfctPARDISO(pardiso::ParDiSO, A::SparsePardisoCSR)
	# Reordering and Symbolic Factorisation + memory allocation for factors

    pardiso.phase = 11;     # just analysis
    ddum   = 0.0;           # double dummy
    idum   = 0;             # integer dummy
    idum32 = int32(0);      # 32-bit integer dummy
    nrhs   = 1;

    pardiso.iparm[33] = 1;   # compute determinant of the matrix if iparm[33]=1;

    println("Analysis of real matrix");

    println();

    ccall( (:pardiso_call_, "/users/verbof/.julia/v0.3/PARDISO/lib/PARDISO"),
        Void,
        (Ptr{Int64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
         Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
         Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64},
         Ptr{Int32}, Ptr{Float64}),
        pardiso.pt, &pardiso.maxfct, &pardiso.mnum, &pardiso.mtype, &pardiso.phase, &(A.n),
        A.nzval, A.rowptr, A.colval, &idum32,
        &idum, pardiso.iparm, &pardiso.msglvl, &ddum, &ddum,
        &pardiso.error_, pardiso.dparm);
    
    if pardiso.error_ != 0
        error("Error in symbolic factorisation of matrix: $(pardiso.error_)");
    end

    errorPARDISO(pardiso.error_);

    println("Symbolic factorisation is complete.");

end



function factorPARDISO(pardiso::ParDiSO, A::SparsePardisoCSR)
	# Generate LU-factorization of a real matrix A

    pardiso.phase = 22;     # just numerical factorisation
    ddum   = 0.0;           # double dummy
    idum32 = int32(0);      # 32-bit integer dummy

    println("Factorize real matrix")
    
    ccall( (:pardiso_call_, "/users/verbof/.julia/v0.3/PARDISO/lib/PARDISO"),
        Void,
        (Ptr{Int64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
        Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
        Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64},
        Ptr{Int32}, Ptr{Float64}),
        pardiso.pt, &pardiso.maxfct, &pardiso.mnum, &pardiso.mtype, &pardiso.phase, &(A.n),
        A.nzval, A.rowptr, A.colval, &idum32,
        &idum32, pardiso.iparm, &pardiso.msglvl, &ddum, &ddum,
        &pardiso.error_, pardiso.dparm);

    errorPARDISO(pardiso.error_);

end



function solvePARDISO(pardiso::ParDiSO, A::SparsePardisoCSR, n_rhs::Int64, b::Vector{Float64})
	# Computes the solution for the system Ax = b, where b is 1 RHS.

    pardiso.phase = 33;        # solve + iterative refinement
    ddum   = 0.0;              # double dummy
    idum32 = int32(0);         # 32-bit integer dummy
    nrhs   = int32(n_rhs);     # number of right-hand-sides
    x      = zeros(nrhs*A.n);  # solution of the system

    pardiso.iparm[6] = 0;      # puts the solution in x (b is NOT overwritten)

    println("Factorize real matrix")

    ccall( (:pardiso_call_, "/users/verbof/.julia/v0.3/PARDISO/lib/PARDISO"),
        Void,
        (Ptr{Int64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
        Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
        Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64},
        Ptr{Int32}, Ptr{Float64}),
        pardiso.pt, &pardiso.maxfct, &pardiso.mnum, &pardiso.mtype, &pardiso.phase, &(A.n),
        A.nzval, A.rowptr, A.colval, &idum32,
        &nrhs, pardiso.iparm, &pardiso.msglvl, b, x,
        &pardiso.error_, pardiso.dparm);


    errorPARDISO(pardiso.error_);

    return x

end



function solvePARDISO(pardiso::ParDiSO, A::SparsePardisoCSR, b::Array{Float64})
	# Computes the solution for the system Ax = b, where b is the COLUMNWISE
    # representation of the RHS's.

    s    = size(b);
    nrhs = 0;

    if length(s) == 0
        error("Provide at least one RHS!");
        
    elseif length(s) == 1
        if s[1]%A.n == 0
            B    = b;
            nrhs = int(s[1]/A.n);
        else
            error("length(b) = ", s[1], " is not a multiple of A.n = ", A.n, "\n");
        end
        
    else
        if s[1] == A.n
            B    = colwise(b);
            nrhs = s[2];
        else
            error("The dimension of the ", s[2], "RHS's provided are not consistent with ", A.n, "\n");
        end
        
    end

    return solvePARDISO(pardiso, A, nrhs, B);

end
