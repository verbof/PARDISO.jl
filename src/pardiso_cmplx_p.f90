module pardiso_cmplx_mod

   save
  ! integer,parameter:: ui_out = 10  ! unit for pardiso log file
   !integer,parameter:: st_unit=11  ! unit for output statistics file


   INCLUDE 'zpardiso_struc.h'

   !logical,private:: pardiso_initialized = .false.


contains

!---------------------------------------------------------------

subroutine init( pardiso_par, sym )

implicit none

TYPE (Zpardiso_STRUC),intent(inout):: pardiso_par
integer(kind=8),intent(in):: sym  ! =0 unsymmetric, =1 symm. pos def, =2 general symm.

   !pardiso_par%COMM = MPI_COMM_WORLD
   pardiso_par%SYM = sym  ! =0 unsymmetric, =1 symm. pos def, =2 general symm.
   pardiso_par%PAR = 1

   pardiso_par%JOB = -1  ! initialize
   CALL Zpardiso(pardiso_par)


   !if (pardiso_par%MYID == 0) then
    !  open(unit=ui_out, file='pardiso.log', action='write')
   !end if

   !pardiso_initialized = .true.

return
end subroutine init

!---------------------------------------------------------------



!---------------------------------------------------------------

subroutine factor_matrix( pardiso_par, ooc, ierr )

implicit none

TYPE (Zpardiso_STRUC),intent(inout):: pardiso_par
integer(kind=8),intent(in):: ooc ! = 0 in-core factorization, = 1 out-of-core factorization
integer(kind=8),intent(out):: ierr

!pardiso_par%icntl(2) = 6  ! output stream for diagnostics
pardiso_par%icntl(4) = 0 ! 1  ! amount of output


pardiso_par%icntl(2) = 0 ! ui_out  ! output stream for diagnostic printing
pardiso_par%icntl(3) = 0 ! ui_out  ! output stream for global information

!pardiso_par%icntl(14) = 40 ! % increase in working space

pardiso_par%icntl(22) = ooc ! out-of-core factorization

!pardiso_par%icntl(11) = 0 ! return statistics
! pardiso_par%icntl(7) = 5  ! ordering type
!pardiso_par%cntl(1) = 0.d0  ! rel. threshold for numerical pivoting


pardiso_factorization: do

   pardiso_par%JOB = 1  !  analysis
   CALL Zpardiso(pardiso_par)

   pardiso_par%JOB = 2  !  factorization
   CALL Zpardiso(pardiso_par)

   ierr = pardiso_par%INFOG(1)


   if ( pardiso_par%INFOG(1) == -9 .or. pardiso_par%INFOG(1) == -8 ) then
      ! Main internal real/complex workarray S too small.
      pardiso_par%icntl(14) = pardiso_par%icntl(14) + 10
      if ( pardiso_par%MYID == 0 ) then
         write(*,30) pardiso_par%icntl(14)
         30 format(/'pardiso percentage increase in the estimated working space',  &
                   /'increased to',i4)
      end if

   !else if (pardiso_par%INFOG(1) == -13) then
   !   if ( pardiso_par%MYID == 0 ) then
   !      write(*,40)
   !      40 format(/'pardiso memory allocation error.'/)
   !   end if
   !   stop

   !else if (pardiso_par%INFOG(1) == -10) then
   !   if ( pardiso_par%MYID == 0 ) then
   !      write(*,45)
   !      45 format(/'pardiso ERROR: Numerically singular matrix.'/)
   !   end if
   !   stop

   !else if (pardiso_par%INFOG(1) == -40) then
   !   if ( pardiso_par%MYID == 0 ) then
   !      write(*,46)
   !      46 format(/'pardiso ERROR: matrix is not positive definite.'/)
   !   end if
   !   stop

   !else if ( pardiso_par%INFOG(1) < 0 ) then
   !   if ( pardiso_par%MYID == 0 ) then
   !      write(*,20) pardiso_par%INFOG(1), pardiso_par%INFOG(2)
   !      20 format(/'ERROR occured in pardiso!',/,'INFOG(1), INFOG(2) ', 2i6,/)
   !   end if
   !   stop

   else
      exit pardiso_factorization ! factorization successful
   end if

end do pardiso_factorization


if ( pardiso_par%MYID == 0 ) then
   !flush(ui_out)
   ! Turn off pardiso output.
   pardiso_par%icntl(2) = 0  ! output stream for diagnostic printing
   pardiso_par%icntl(3) = 0  ! output stream for global information
end if


return
end subroutine factor_matrix

!--------------------------------------------------------------

subroutine solve( pardiso_par, nrhs, rhs, x, transpose )
! Solve A*x = rhs

implicit none

TYPE (Zpardiso_STRUC),intent(inout):: pardiso_par
integer(kind=8),intent(in):: nrhs   ! # of right-hand-sides
complex(kind=8),intent(in):: rhs(nrhs * pardiso_par%N)
complex(kind=8),intent(out),target:: x(nrhs * pardiso_par%N)  ! solution
logical,intent(in):: transpose  ! if .true. take the transpose

   x = rhs

   pardiso_par%RHS => x

   ! The following is significant only on the host cpu.
   pardiso_par%NRHS = nrhs  ! # of right-hand-sides
   pardiso_par%LRHS = pardiso_par%N  ! size of system

   if (transpose) then
      pardiso_par%icntl(9) = 0  ! for solving A'x = b
   else
      pardiso_par%icntl(9) = 1  ! for solving Ax = b
   end if

   if (transpose) then
      pardiso_par%RHS = conjg(pardiso_par%RHS)
   end if

   pardiso_par%JOB = 3  ! Solve the system.
   CALL Zpardiso(pardiso_par)
   ! At this point pardiso_par%RHS (rhs) contains the solution.

   if (transpose) then
      pardiso_par%RHS = conjg(pardiso_par%RHS)
   end if

return
end subroutine solve

!--------------------------------------------------------------

subroutine destroy(pardiso_par)

implicit none
TYPE (Zpardiso_STRUC),intent(inout):: pardiso_par

!if (.not. pardiso_initialized)  return

!  Destroy the instance (deallocate internal data structures)
pardiso_par%JOB = -2
CALL Zpardiso(pardiso_par)





return
end subroutine destroy


end module pardiso_cmplx_mod
