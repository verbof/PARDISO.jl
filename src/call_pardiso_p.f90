!-------------------------------------------------
!-------------------------------------------------

subroutine pardiso_init(pt, mtype, solver, iparm, dparm, ierr)
!DIR$ ATTRIBUTES DLLEXPORT :: pardiso_init
!DIR$ ATTRIBUTES ALIAS: 'pardiso_init_':: pardiso_init

implicit none

TYPE PARDISO_STRUC
   INTEGER :: NN
END TYPE PARDISO_STRUC 


integer(kind=8),intent(out):: pt(64)                ! pardiso pointer
integer(kind=4),intent(in):: mtype
integer(kind=4),intent(in):: solver
integer(kind=4),intent(out):: iparm(64) 
real(kind=8),intent(out):: dparm(64)

integer(kind=4),intent(out):: ierr

integer:: k


!  .. PARDISO license check and initialize solver

      CALL pardisoinit(pt, mtype, solver, iparm, dparm, ierr)

!  .. Numbers of Processors ( value of OMP_NUM_THREADS )

    !DO k = 1,64
    !    iparm(k) = k;
    !    !pt(k)    = k;
    !END DO

    !DO k = 1,64
    !    dparm(k) = 0.0 + k;
    !END DO

        IF (ierr .NE. 0) THEN
        IF (ierr.EQ.-10 ) WRITE(*,*) 'No license file found'
        IF (ierr.EQ.-11 ) WRITE(*,*) 'License is expired'
        IF (ierr.EQ.-12 ) WRITE(*,*) 'Wrong username or hostname'
        STOP
      ELSE
        WRITE(*,*) '[PARDISO]: License check was successful ... '
      END IF

      WRITE(*,*) '[PARDISO]: Initialisation completed!'
     
end subroutine pardiso_init

!-------------------------------------------------
!-------------------------------------------------

subroutine pardiso_checkmatrix(mtype, n, A, iA, jA, ierr)
!DIR$ ATTRIBUTES DLLEXPORT :: pardiso_checkmatrix
!DIR$ ATTRIBUTES ALIAS: 'pardiso_checkmatrix_':: pardiso_checkmatrix

implicit none

TYPE PARDISO_STRUC
   INTEGER :: NN
END TYPE PARDISO_STRUC 


integer(kind=4),intent(in):: mtype
integer(kind=4),intent(in):: n
real(kind=8),intent(out):: A(*)
integer(kind=4),intent(in):: ja(*), iA(n+1)
integer(kind=4),intent(out):: ierr

!
!  .. pardiso_chk_matrix(...)
!     Checks the consistency of the given matrix.
!     Use this functionality only for debugging purposes

    WRITE(*,*) 'PARDISO is checking the matrix...'

      CALL pardiso_chkmatrix(mtype, n, A, iA, jA, ierr);

      IF (ierr .NE. 0) THEN
         WRITE(*,*) 'The following ERROR was detected: ', ierr
         STOP
      ENDIF
     
end subroutine pardiso_checkmatrix

!-------------------------------------------------
!-------------------------------------------------

subroutine pardiso_call( pt, maxfct, mnum, mtype, phase, n, A, jA, iA, perm, nrhs, iparm, msglvl, b, x, ierr, dparm)
!DIR$ ATTRIBUTES DLLEXPORT :: pardiso_call
!DIR$ ATTRIBUTES ALIAS: 'pardiso_call_':: pardiso_call

implicit none

TYPE PARDISO_STRUC
   INTEGER :: NN
END TYPE PARDISO_STRUC 

integer(kind=8)::pt(64)

integer(kind=4),intent(in):: maxfct, mnum, mtype, phase 
integer(kind=4),intent(in):: n
real(kind=8),intent(in):: A(*)
integer(kind=4),intent(in):: jA(*), iA(n+1)

integer(kind=4),intent(in):: perm
integer(kind=4),intent(in):: nrhs, msglvl
integer(kind=4):: iparm(64)
real(kind=8):: dparm
real(kind=8):: b(*), x(*)


integer(kind=4),intent(out):: ierr

      CALL pardiso (pt, maxfct, mnum, mtype, phase, n, A, iA, jA, perm, nrhs, iparm, msglvl, b, x, ierr, dparm)
     
      WRITE(*,*) 'Reordering completed ... '

      IF (ierr .NE. 0) THEN
        WRITE(*,*) 'The following ERROR was detected: ', ierr
        STOP
      END IF

      !WRITE(*,*) 'Number of nonzeros in factors   = ',iparm(18)
      !WRITE(*,*) 'Number of factorization MFLOPS  = ',iparm(19)

end subroutine pardiso_call

!-------------------------------------------------




