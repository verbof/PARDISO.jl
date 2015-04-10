subroutine pardiso_checkmatrix_z(mtype, n, A, iA, jA, ierr)
!DIR$ ATTRIBUTES DLLEXPORT :: pardiso_checkmatrix_z
!DIR$ ATTRIBUTES ALIAS: 'pardiso_checkmatrix_z_':: pardiso_checkmatrix_z

implicit none

integer(kind=4),intent(in):: mtype
integer(kind=4),intent(in):: n
complex(kind=8),intent(out):: A(*)
integer(kind=4),intent(in):: ja(*), iA(n+1)
integer(kind=4),intent(out):: ierr

!integer:: i, j

        !do i =1, n
        !  do j=ia(i), ia(i+1)-1
        !    write(*,*) i, ja(j), a(j)
        !  end do
        !end do

!
!  .. pardiso_chk_matrix(...)
!     Checks the consistency of the given matrix.
!     Use this functionality only for debugging purposes

    WRITE(*,*) 'PARDISO is checking the matrix...'

      CALL pardiso_chkmatrix_z(mtype, n, A, iA, jA, ierr);

      IF (ierr .NE. 0) THEN
         WRITE(*,*) 'The following ERROR was detected: ', ierr
         STOP
      ENDIF
     
end subroutine pardiso_checkmatrix_z

!-------------------------------------------------
!-------------------------------------------------

subroutine pardiso_call_z(pt, maxfct, mnum, mtype, phase, n, A, iA, jA, perm, nrhs, iparm, msglvl, b, x, ierr, dparm)
!DIR$ ATTRIBUTES DLLEXPORT :: pardiso_call_z
!DIR$ ATTRIBUTES ALIAS: 'pardiso_call_z_':: pardiso_call_z

implicit none

integer(kind=8)::pt(64)

integer(kind=4),intent(in):: maxfct, mnum, mtype, phase 
integer(kind=4),intent(in):: n
complex(kind=8),intent(in):: A(*)
integer(kind=4),intent(in):: jA(*), iA(n+1)

integer(kind=4),intent(in):: perm
integer(kind=4),intent(in):: nrhs, msglvl
integer(kind=4):: iparm(64)
real(kind=8):: dparm(64)
complex(kind=8):: b(*), x(*)

integer(kind=4),intent(out):: ierr


      CALL pardiso(pt, maxfct, mnum, mtype, phase, n, A, iA, jA, perm, nrhs, iparm, msglvl, b, x, ierr, dparm)
     
      WRITE(*,*)

      IF (ierr .NE. 0) THEN
        WRITE(*,*) 'The following ERROR was detected: ', ierr
        STOP
      END IF

      !WRITE(*,*) 'Number of nonzeros in factors   = ', iparm(18)
      !WRITE(*,*) 'Number of factorization MFLOPS  = ', iparm(19)
      !IF (iparm(33) .EQ. 1) THEN
      !  WRITE(*,*) 'log(det(A)) = ', dparm(33)
      !ENDIF

end subroutine pardiso_call_z

!-------------------------------------------------

