
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



function checkPARDISO{Tnzval}(pardiso::ParDiSO, A::SparsePardisoCSR{Tnzval})

    if abs(pardiso.mtype) == 2 || 
       abs(pardiso.mtype) == 4 || 
           pardiso.mtype  == 6
        
        if !A.upper
            error("pardiso.mtype = ±2 or ±4 or +6, but the matrix is not in upper triangular form (A.upper == false)");
        end
        
    end

    if Tnzval == Float64
        ccall( (:pardiso_checkmatrix_, "/users/verbof/.julia/v0.3/PARDISO/lib/PARDISO"), Void,
                ( Ptr{Int32},   Ptr{Int32}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}),
                &pardiso.mtype, &(A.n),     A.nzval,      A.rowptr,   A.colval,   &pardiso.error_);

    elseif Tnzval == Complex128
        ccall( (:pardiso_checkmatrix_z_, "/users/verbof/.julia/v0.3/PARDISO/lib/PARDISO"), Void,
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

    #pardiso.iparm[33] = 1;   # compute determinant of the matrix if iparm[33]=1;

    if Tnzval == Float64
        ccall( (:pardiso_call_, "/users/verbof/.julia/v0.3/PARDISO/lib/PARDISO"),
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
        ccall( (:pardiso_call_z_, "/users/verbof/.julia/v0.3/PARDISO/lib/PARDISO"),
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
	# Generate LU-factorization of a real matrix A

    pardiso.phase = 22;     # just numerical factorisation
    ddum = 0.0;             # double dummy
    idum = 0;               # 32-bit integer dummy

    
    if Tnzval == Float64
        println("Factorizing real matrix...")
        
        ccall( (:pardiso_call_, "/users/verbof/.julia/v0.3/PARDISO/lib/PARDISO"),
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

        ccall( (:pardiso_call_z_, "/users/verbof/.julia/v0.3/PARDISO/lib/PARDISO"),
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



function solvePARDISO{Tnzval}(pardiso::ParDiSO, A::SparsePardisoCSR{Tnzval}, n_rhs::Int64, b::Vector{Tnzval})
	# Computes the solution for the system Ax = b, where b is 1 RHS.

    pardiso.phase = 33;                 # solve + iterative refinement
    ddum   = 0.0;                       # double dummy
    idum   = 0;                         # integer dummy
    nrhs   = n_rhs;                     # number of right-hand-sides
    x      = zeros(Tnzval, nrhs*A.n);   # solution of the system

    pardiso.iparm[6] = 0;      # puts the solution in x (b is NOT overwritten)

    if Tnzval == Float64
        println("Solving real system...");

        ccall( (:pardiso_call_, "/users/verbof/.julia/v0.3/PARDISO/lib/PARDISO"),
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
    
        ccall( (:pardiso_call_z_, "/users/verbof/.julia/v0.3/PARDISO/lib/PARDISO"),
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
#
# The following function checks if the RHS provided is a column vector
# containing - at least one - RHS('s) or the RHS's are provided a a matrix. In
# the second case, the vector is reshaped.
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
	# Computes the solution for the system Ax = b, where b is the COLUMNWISE
    # representation of the RHS's.

    s    = size(b);
    nrhs = 0;

    if length(s) == 0
        error("Provide at least one RHS!");
        
    elseif length(s) == 1           # if b is column vector
        if s[1]%A.n == 0            # check if b has a correct number of entries
            B    = b;
            nrhs = int(s[1]/A.n);   # compute the number of RHS's stored in b
        else
            error("length(b) = ", s[1], " is not a multiple of A.n = ", A.n, "\n");
        end
        
    else                            # if b is not a column vector
        if s[1] == A.n              # check if it has a correct number of rows
            B    = colwise(b);      # create the columnwise version of b and store it in B
            nrhs = s[2];            # the number of RHS is the number of columns of b
        else
            error("The dimension of the ", s[2], "RHS's provided are not consistent with ", A.n, "\n");
        end
        
    end


    if eltype(A) == Complex128 && eltype(B) == Float64
        B = complex(B);
    end


    return solvePARDISO(pardiso, A, nrhs, B);

end


function memoryPARDISO(pardiso::ParDiSO)

    # Returns the peak memory, in kB, consumption by PARDISO.

    return max(pardiso.iparm[15], pardiso.iparm[16]+pardiso.iparm[17]);

end



function freePARDISO(pardiso::ParDiSO)
    
    # Free all the PARDISO memory.

    idum = int32(0);
    ddum = 0.0;

    pardiso.phase = int32(-1);

    ccall( (:pardiso_call_, "/users/verbof/.julia/v0.3/PARDISO/lib/PARDISO"),
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

