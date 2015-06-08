
# export errorPARDISO,
#        initPARDISO,
#        checkPARDISO,
#        smbfctPARDISO,
#        factorPARDISO,
#        solvePARDISO,
#        memoryPARDISO
#        freePARDISO;



function errorPARDISO(error_)

    if error_ != 0
  
        print("ParDiSO error: $(error_) --- ");

        if error_ == -1
            print("Input inconsistent\n");
        elseif error_ == -2
            print("Not enough memory\n");
        elseif error_ == -3
            print("Reordering problem\n");
        elseif error_ == -4
            print("Zero pivot, numerical factorization or iterative refinement problem\n");
        elseif error_ == -5
            print("Unclassified (internal) error\n");
        elseif error_ == -6
            print("Preordering failed (matrix type 11, 13 only)\n");
        elseif error_ == -7
            print("Diagonal matrix problem\n");
        elseif error_ == -8
            print("32-bit integer overflow problem\n");
        elseif error_ == -10
            print("No license file pardiso.lic found\n");
        elseif error_ == -11
            print("License is expired.\n");
        elseif error_ == -12
            print("Wrong username or hostname.\n");
        end

    else
        # Nothing!
    end

end


function initPARDISO(pardiso::ParDiSO)


    ccall( (:pardiso_init_, PARDISOLIB), Void,
           (Ptr{Int64}, Ptr{Int32},     Ptr{Int32},      Ptr{Int32},    Ptr{Float64},  Ptr{Int32}),
            pardiso.pt, &pardiso.mtype, &pardiso.solver, pardiso.iparm, pardiso.dparm, &pardiso.error_);

    errorPARDISO(pardiso.error_);

    if "OMP_NUM_THREADS" in keys(ENV)
        pardiso.iparm[3] = Int32(parse(ENV["OMP_NUM_THREADS"]));
    else
        error("Set environment variable OMP_NUM_THREADS before running this code.");
    end

    
end



function checkPARDISO{Tnzval}(pardiso::ParDiSO, A::SparsePardisoCSR{Tnzval})

    if abs(pardiso.mtype) == 2 ||   #    real Symmetric Pos.Def. / Undef.
       abs(pardiso.mtype) == 4 ||   # complex Hermitian Pos.Def. / Undef.
           pardiso.mtype  == 6      # complex Symmetric

        if !A.upper                 # Checks if the matrix is stored in upper triangular format. If not, error(...).
            error("pardiso.mtype = ±2 or ±4 or +6, but the matrix is not in upper triangular form (A.upper == false)");
        end
        
    end

    if Tnzval == Float64
        ccall( (:pardiso_checkmatrix_, PARDISOLIB), Void,
                ( Ptr{Int32},   Ptr{Int32}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}),
                &pardiso.mtype, &(A.n),     A.nzval,      A.rowptr,   A.colval,   &pardiso.error_);

    elseif Tnzval == Complex128
        ccall( (:pardiso_checkmatrix_z_, PARDISOLIB), Void,
                ( Ptr{Int32},   Ptr{Int32}, Ptr{Complex128}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}),
                &pardiso.mtype, &(A.n),     A.nzval,         A.rowptr,   A.colval,   &pardiso.error_);
    end

    if pardiso.error_ != 0
        error("Error in consistency of matrix: $(pardiso.error_)");
    end

    errorPARDISO(pardiso.error_);

end


function smbfctPARDISO{Tnzval}(pardiso::ParDiSO, A::SparsePardisoCSR{Tnzval})
	# Reordering and Symbolic Factorisation + memory allocation for factors

    pardiso.phase = 11;     # just analysis
    ddum   = 0.0;           # double dummy
    idum   = 0;             # integer dummy
    nrhs   = 1;             # number of RHS's


    if Tnzval == Float64
        ccall( (:pardiso_call_, PARDISOLIB),
            Void,
            (Ptr{Int64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
            Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
            Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64},
            Ptr{Int32}, Ptr{Float64}),
            pardiso.pt, &pardiso.maxfct, &pardiso.mnum, &pardiso.mtype, &pardiso.phase, &(A.n),
            A.nzval, A.rowptr, A.colval, &idum,
            &idum, pardiso.iparm, &pardiso.msglvl, &ddum, &ddum,
            &pardiso.error_, pardiso.dparm);

    elseif Tnzval == Complex128
        ccall( (:pardiso_call_z_, PARDISOLIB),
            Void,
            (Ptr{Int64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
            Ptr{Complex128}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
            Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Complex128}, Ptr{Complex128},
            Ptr{Int32}, Ptr{Float64}),
            pardiso.pt, &pardiso.maxfct, &pardiso.mnum, &pardiso.mtype, &pardiso.phase, &(A.n),
            A.nzval, A.rowptr, A.colval, &idum,
            &idum, pardiso.iparm, &pardiso.msglvl, &ddum, &ddum,
            &pardiso.error_, pardiso.dparm);

    end

    if pardiso.error_ != 0
        error("Error in symbolic factorisation of matrix: $(pardiso.error_)");
    end

    errorPARDISO(pardiso.error_);

    println("Symbolic factorisation completed.");

end



function factorPARDISO{Tnzval}(pardiso::ParDiSO, A::SparsePardisoCSR{Tnzval})
	# Generate LU-factorisation of a real matrix A

    pardiso.phase = 22;     # just numerical factorisation
    ddum = 0.0;             # double dummy
    idum = 0;               # 32-bit integer dummy

    
    if Tnzval == Float64
        println("Factorizing real matrix...")
        
        ccall( (:pardiso_call_, PARDISOLIB),
            Void,
            (Ptr{Int64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
            Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
            Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64},
            Ptr{Int32}, Ptr{Float64}),
            pardiso.pt, &pardiso.maxfct, &pardiso.mnum, &pardiso.mtype, &pardiso.phase, &(A.n),
            A.nzval, A.rowptr, A.colval, &idum,
            &idum, pardiso.iparm, &pardiso.msglvl, &ddum, &ddum,
            &pardiso.error_, pardiso.dparm);

    elseif Tnzval == Complex128
        println("Factorizing complex matrix...")

        ccall( (:pardiso_call_z_, PARDISOLIB),
            Void,
            (Ptr{Int64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
            Ptr{Complex128}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
            Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Complex128}, Ptr{Complex128},
            Ptr{Int32}, Ptr{Float64}),
            pardiso.pt, &pardiso.maxfct, &pardiso.mnum, &pardiso.mtype, &pardiso.phase, &(A.n),
            A.nzval, A.rowptr, A.colval, &idum,
            &idum, pardiso.iparm, &pardiso.msglvl, &ddum, &ddum,
            &pardiso.error_, pardiso.dparm);

    end

    errorPARDISO(pardiso.error_);

    println("Factorization completed.");

end



function solvePARDISO{Tnzval}(pardiso::ParDiSO, A::SparsePardisoCSR{Tnzval}, n_rhs::Int64, b::Array{Tnzval})
	# Computes the solution for the system Ax = b, where b is 1 RHS.

    pardiso.phase = 33;                 # solve + iterative refinement
    ddum   = 0.0;                       # double dummy
    idum   = 0;                         # integer dummy
    nrhs   = n_rhs;                     # number of right-hand-sides
    x      = zeros(Tnzval, nrhs*A.n);   # solution of the system

    pardiso.iparm[6] = 0;               # Overwrites b with the soluzion if iparm[6] = 1;

    if Tnzval == Float64
        println("Solving real system...");

        ccall( (:pardiso_call_, PARDISOLIB),
            Void,
            (Ptr{Int64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
            Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
            Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64},
            Ptr{Int32}, Ptr{Float64}),
            pardiso.pt, &pardiso.maxfct, &pardiso.mnum, &pardiso.mtype, &pardiso.phase, &(A.n),
            A.nzval, A.rowptr, A.colval, &idum,
            &nrhs, pardiso.iparm, &pardiso.msglvl, b, x,
            &pardiso.error_, pardiso.dparm);

    elseif Tnzval == Complex128
        println("Solving complex system...");
    
        ccall( (:pardiso_call_z_, PARDISOLIB),
            Void,
            (Ptr{Int64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
            Ptr{Complex128}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
            Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Complex128}, Ptr{Complex128},
            Ptr{Int32}, Ptr{Float64}),
            pardiso.pt, &pardiso.maxfct, &pardiso.mnum, &pardiso.mtype, &pardiso.phase, &(A.n),
            A.nzval, A.rowptr, A.colval, &idum,
            &nrhs, pardiso.iparm, &pardiso.msglvl, b, x,
            &pardiso.error_, pardiso.dparm);
    end

    errorPARDISO(pardiso.error_);

    return x

end



#=============================================================================
#                       NOT ANYMORE (IT IS USELESS)
# The following function checks if the variable b provided is a column vector
# containing - at least one - some RHS's or the RHS's are provided a a matrix. 
# In the second case, the vector is reshaped.
#
#         1   2  ... nrhs
#       +---+---+---+---+                   
#       | * | * |   | * |                         1           2           3
#       | * | * |   | * |                   +-----------+-----------+------···
# b  =  | * | * |...| * |   --->  RESHAPED: | * * * * * | * * * * * | * * * ···
#       | * | * |   | * |                   +-----------+-----------+------···
#       | * | * |   | * |
#       +---+---+---+---+
#
=============================================================================#

function solvePARDISO(pardiso::ParDiSO, A::SparsePardisoCSR, b::Array)
	# Computes the correct number of RHS's for an arbitrary vector b;

    s    = size(b);
    nrhs = 0;

    if length(s) == 0
        error("Provide at least one RHS!");
        
    elseif length(s) == 1           # if b is column vector
        if s[1]%A.n == 0            # check if b has a correct number of entries
            nrhs = Int(s[1]/A.n);   # compute the number of RHS's stored in b
        else
            error("length(b) = ", s[1], " is not a multiple of A.n = ", A.n, "\n");
        end
        
    else                            # if b is not a column vector
        if s[1] == A.n              # check if it has a correct number of rows
            nrhs = s[2];            # the number of RHS is the number of columns of b
        else
            error("The dimension of the ", s[2], "RHS's provided is not consistent with A.n = ", A.n, "\n");
        end
        
    end


    if eltype(A) == Complex128 && eltype(b) == Float64
        return solvePARDISO(pardiso, A, nrhs, complex(b));
    end
    
    return solvePARDISO(pardiso, A, nrhs, float(b));

end


function memoryPARDISO(pardiso::ParDiSO)
    # Returns the peak memory consumption by PARDISO, in kB (see PARDISO Manual).

    return max(pardiso.iparm[15], pardiso.iparm[16]+pardiso.iparm[17]);

end



function freePARDISO(pardiso::ParDiSO)
    # Releases all the PARDISO memory.

    idum = Int32(0);            # integer dummy
    ddum = 0.0;                 # double dummy

    pardiso.phase = Int32(-1);  # Release internal memory for all matrices

    ccall( (:pardiso_call_, PARDISOLIB),
        Void,
        (Ptr{Int64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
        Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
        Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64},
        Ptr{Int32}, Ptr{Float64}),
        pardiso.pt, &pardiso.maxfct, &pardiso.mnum, &pardiso.mtype, &pardiso.phase, &idum,
        &ddum, &idum, &idum, &idum,
        &idum, pardiso.iparm, &pardiso.msglvl, &ddum, &ddum,
        &pardiso.error_, pardiso.dparm);

    errorPARDISO(pardiso.error_);

end

