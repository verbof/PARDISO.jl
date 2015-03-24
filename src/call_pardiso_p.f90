!-------------------------------------------------

subroutine pardiso_init(pt, mtype, solver, iparm, dparm, ierr)
!DIR$ ATTRIBUTES DLLEXPORT :: pardiso_init
!DIR$ ATTRIBUTES ALIAS: 'pardiso_init_':: pardiso_init

implicit none

TYPE PARDISO_STRUC
   INTEGER :: NN
END TYPE PARDISO_STRUC 


integer(kind=8):: pt                ! pardiso pointer
integer(kind=8),intent(in):: mtype
integer(kind=8),intent(in):: solver
integer(kind=8),intent(out):: iparm(64) 
real(kind=8),intent(out):: dparm(64)

integer(kind=8),intent(out):: ierr


!  .. PARDISO license check and initialize solver

      CALL pardisoinit(pt, mtype, solver, iparm, dparm, ierr)

!  .. Numbers of Processors ( value of OMP_NUM_THREADS )

      iparm(3) = 1

      IF (ierr .NE. 0) THEN
        IF (ierr.EQ.-10 ) WRITE(*,*) 'No license file found'
        IF (ierr.EQ.-11 ) WRITE(*,*) 'License is expired'
        IF (ierr.EQ.-12 ) WRITE(*,*) 'Wrong username or hostname'
        STOP
      ELSE
        WRITE(*,*) '[PARDISO]: License check was successful ... '
      END IF
     
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


integer(kind=8),intent(in):: mtype
integer(kind=8),intent(in):: n
real(kind=8),intent(out):: A(*)
integer(kind=8),intent(in):: ja(*), iA(n+1)
integer(kind=8),intent(out):: ierr

!  .. pardiso_chk_matrix(...)
!     Checks the consistency of the given matrix.
!     Use this functionality only for debugging purposes

      CALL pardiso_chkmatrix  (mtype, n, A, iA, jA, ierr);

      IF (ierr .NE. 0) THEN
         WRITE(*,*) 'The following ERROR was detected: ', ierr
        STOP
      ENDIF
     
end subroutine pardiso_checkmatrix

!-------------------------------------------------
!-------------------------------------------------

subroutine pardiso_call( pt, maxfct, mnum, mtype, phase, n, A, jA, iA, perm, nrhs, iparm, msglvl, b, x, ierr, dparm)
!DIR$ ATTRIBUTES DLLEXPORT :: factor_pardiso
!DIR$ ATTRIBUTES ALIAS: 'factor_pardiso_':: factor_pardiso

implicit none

TYPE PARDISO_STRUC
   INTEGER :: NN
END TYPE PARDISO_STRUC 

integer(kind=8)::pt(64)

integer(kind=8),intent(in):: maxfct, mnum, mtype, phase 
integer(kind=8),intent(in):: n
real(kind=8),intent(in):: A(*)
integer(kind=8),intent(in):: jA(*), iA(n+1)

integer(kind=8),intent(in):: perm, nrhs, msglvl
integer(kind=8):: iparm
real(kind=8):: dparm
real(kind=8):: b, x


integer(kind=8),intent(out):: ierr

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




