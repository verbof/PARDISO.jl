function solvePARDISO(A::SparseMatrixCSC, rhs::Array, x::Array=[], sym=0,ooc=0,tr=0)

error("debugging !!!")
    error("debugging !!!")
    error("debugging !!!")

	n    = size(rhs,1)
	nrhs = size(rhs,2)
	
	if isempty(x);
		x = zeros(eltype(A),n,nrhs)
	else
		if size(x)!=(n,nrhs); 
			error("applyPARDISO: wrong size of x provided")
		end
		# make x complex if necessary
		x = isreal(A) ? x : complex(x)
	end
	
	if norm(rhs)==0
		x = rhs
		return x
	end
	
	# factorization
	factor = factorPARDISO(A,sym,ooc)
	
	# solve system
	x = applyPARDISO(factor,rhs,x,tr)
	
	# free memory
	destroyPARDISO(factor)
   return x
end

function factorPARDISO(A::SparseMatrixCSC{Complex128,Int},sym=0,ooc=0)
	# Generate LU-factorization of complex matrix A, 'sym' can be: 0=unsymmetric, 1=symm. pos def, 2=general symmetric

error("debugging !!!")
    error("debugging !!!")
    error("debugging !!!")

  

    if size(A,1) != size(A,2)
		error("factorPARDISO: Matrix must be square!")
  	end
	n  = size(A,1);
 	PARDISOstat = [0];
	tic()
    p  = ccall( (:factor_pardiso_cmplx_, "/users/verbof/.julia/v0.3/PARDISO/lib/PARDISO"),
    		 Int64, ( Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Complex128}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}),
             &n, &sym, &ooc,  convert(Ptr{Complex128}, pointer(A.nzval)), A.rowval, A.colptr, PARDISOstat);
    checkPARDISOerror(PARDISOstat);
	facTime = toq()
	factor = PARDISOfactorization(p,myid(),n,false,facTime);
	# finalizer(factor,destroyPARDISO)
	return factor
end

function factorPARDISO(A::SparseMatrixCSC{Float64,Int64},sym=0,ooc=0)
	# Generate LU-factorization of a real matrix A, 'sym' can be: 0=unsymmetric, 1=symm. pos def, 2=general symmetric


    rowval = convert(Array{Int32}, A.rowval);
    colptr = convert(Array{Int32}, A.colptr);

    if size(A,1) != size(A,2)
		error("factorPARDISO: Matrix must be square!")
   	end
    println("olaf Factorize real matrix")
	n  = size(A,1);
    PARDISOstat = [0];

    pt = 0;

    tic()
	p  = ccall( (:factor_pardiso_, "/users/verbof/PARDISO/lib/PARDISO"),
                Int64, 
                ( Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}),
                &pt,         &n,           &sym,       &ooc,       A.nzval,      rowval,   colptr,   PARDISOstat);

    checkPARDISOerror(PARDISOstat);
	facTime = toq()
	factor = PARDISOfactorization(p,myid(),n,true,facTime);
	# finalizer(factor,destroyPARDISO)
	return factor
end


function checkPARDISOerror(PARDISOstat)
	# check if PARDISO reported error
	
	if PARDISOstat[1]==-10
	     error("PARDISO: Numerically singular matrix.");
	elseif PARDISOstat[1]==-13
	     error("PARDISO: memory allocation error")
	elseif PARDISOstat[1]==-40
	     error("PARDISO: matrix is not positive definite")
	elseif PARDISOstat[1]==-90
	     error("PARDISO: Error in out-of-core management.Probably there is not enough disk space.")
	elseif PARDISOstat[1]<0
	     error( @sprintf("PARDISO: error --> %d <--. Please refer to Ch. 7 of PARDISO User's guide!",PARDISOstat[1]))
	end
end

function applyPARDISO(factor::PARDISOfactorization,rhs::Array,x::Array=[],tr=0)
	# solve for right hand side(s)

	id1 = myid(); id2 = factor.worker
	if id1 != id2
		warn("Worker $id1 has no access to PARDISO factorization stored on $id2. Trying to remotecall!")
		return remotecall_fetch(factor.worker,applyPARDISO,factor,rhs,x,tr)
	end
	
	 n    = size(rhs,1)
	 nrhs = size(rhs,2)
	 rhs  = reshape(rhs,n,nrhs)
	
	# check size of rhs, allocate space for x
	if n != factor.n;  error("applyPARDISO: wrong size of rhs"); end
	x = isempty(x) ? x = zeros(n,nrhs) : x
	if size(x)!=(n,nrhs); error("applyPARDISO: wrong size of x provided"); end
	 
	if factor.real && !isreal(rhs)
		error("applyPARDISO: rhs must be real.")	
	end
	
	
	# make x complex if necessary
	x   = factor.real ? x   : complex(x)
	rhs = factor.real ? rhs : complex(rhs)
	
	ptr = factor.ptr
	if factor.real
		ccall( (:solve_pardiso_, "/users/verbof/.julia/v0.3/PARDISO/lib/PARDISO"),
				Int64, (Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int64} ),
							&ptr,      &nrhs, 	 	rhs, x, &tr)
	else
		ccall( (:solve_pardiso_cmplx_, "`/users/verbof/.julia/v0.3/PARDISO/lib/PARDISO"),
					Int64, (Ptr{Int64}, Ptr{Int64}, Ptr{Complex128}, Ptr{Complex64}, Ptr{Int64} ),
					&ptr,        &nrhs,   convert(Ptr{Complex128}, pointer(rhs)),   convert(Ptr{Complex128}, pointer(x)),   &tr)
	end
	return x
end

function destroyPARDISO(factor::PARDISOfactorization)
	 #  free memory
	id1 = myid(); id2 = factor.worker;
	if myid() != factor.worker
		warn("Worker $id1 cannot destroy PARDISO factorization stored on worker $id2. Trying remotecall_fetch.")
		remotecall_fetch(id2,destroyPARDISO,factor)
		return
	end
	if factor.real
		ccall( (:destroy_pardiso_, "/users/verbof/.julia/v0.3/PARDISO/lib/PARDISO"),
				Int64, (Ptr{Int64}, ), &factor.ptr )
	else
	 	ccall( (:destroy_pardiso_cmplx_, "/users/verbof/.julia/v0.3/PARDISO/lib/PARDISO"),
	 			Int64, (Ptr{Int64}, ), &factor.ptr )
	end
	factor.ptr = -1
	factor.n    = -1
	factor.time = -1.0

end
