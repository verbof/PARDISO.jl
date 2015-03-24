


module pardiso_mod

   save
  ! integer,parameter:: ui_out = 10  ! unit for pardiso log file
   !integer,parameter:: st_unit=11  ! unit for output statistics file


   INCLUDE 'dpardiso_struc.h'

   !logical,private:: pardiso_initialized = .false.


contains

!---------------------------------------------------------------

subroutine init( pardiso_par, sym )

implicit none

TYPE (Dpardiso_STRUC),intent(inout):: pardiso_par
integer(kind=8),intent(in):: sym  ! =0 unsymmetric, =1 symm. pos def, =2 general symm.x

   !pardiso_par%COMM = MPI_COMM_WORLD
   pardiso_par%SYM = sym  ! =0 unsymmetric, =1 symm. pos def, =2 general symm.
   pardiso_par%PAR = 1

   pardiso_par%JOB = -1  ! initialize
   CALL Dpardiso(pardiso_par)


   !if (pardiso_par%MYID == 0) then
    !  open(unit=ui_out, file='pardiso.log', action='write')
   !end if

   !pardiso_initialized = .true.

return
end subroutine init

!---------------------------------------------------------------

subroutine convert_to_pardiso_format( pardiso_par, n, A,jA,iA, ierr )
! A is actually a transpose.
implicit none

TYPE (Dpardiso_STRUC),intent(inout):: pardiso_par
integer(kind=8),intent(in):: n  ! # of rows in A
real(kind=8),intent(in):: A(*)
integer(kind=8),intent(in):: jA(*), iA(n+1)
integer(kind=8),intent(out):: ierr

integer nonz, i,j, j1,j2, ind, jcol, istat


nonz = iA(n+1) - 1
!n = size(iA) - 1

pardiso_par%N = n

if (pardiso_par%SYM == 0) then  ! not symmetric matrix
   pardiso_par%NZ = nonz
   allocate( pardiso_par%IRN(nonz), pardiso_par%JCN(nonz), pardiso_par%A(nonz), stat=istat)
   if ( istat /= 0 ) then
      ierr = -13  ! allocation error
      return
   end if


   do i = 1, n
      pardiso_par%JCN( iA(i) : iA(i+1)-1 )  =  i
   end do  ! i

   pardiso_par%IRN = jA(1:nonz)
   pardiso_par%A    = A(1:nonz)

else  ! symmetric matrix
   ! Keep only the lower half of the matrix (row >= column).

   nonz = nonz/2 + n  ! should be +n/2, but I'm using n just in case.
   allocate( pardiso_par%IRN(nonz), pardiso_par%JCN(nonz), pardiso_par%A(nonz), stat=istat)
   if ( istat /= 0 ) then
      ierr = -13  ! allocation error
      return
   end if

   ind = 0
   do i = 1, n

      j1 = iA(i)
      j2 = iA(i+1) - 1
      do j = j1, j2
         jcol = jA(j)

         if (i >= jcol) then
            ind = ind + 1
            pardiso_par%A(ind) = A(j)
            pardiso_par%JCN(ind) = jcol
            pardiso_par%IRN(ind) = i
         end if

      end do  ! j
   end do  ! i

   pardiso_par%NZ = ind
   !if (nonz < ind) call errorstop('nonz < ind')  ! debug
end if

ierr = 0
return
end subroutine convert_to_pardiso_format

!---------------------------------------------------------------

subroutine factor_matrix( pardiso_par, ooc, ierr )

implicit none

TYPE (Dpardiso_STRUC),intent(inout):: pardiso_par
integer(kind=8),intent(in):: ooc ! = 0 in-core factorization, = 1 out-of-core factorization
integer(kind=8),intent(out):: ierr

pardiso_par%icntl(2) = -1  ! output stream for diagnostics
pardiso_par%icntl(4) = -1 ! 1  ! amount of output


pardiso_par%icntl(2) = 0 ! ui_out  ! output stream for diagnostic printing
pardiso_par%icntl(3) = 0 ! ui_out  ! output stream for global information

!pardiso_par%icntl(14) = 40 ! % increase in working space

pardiso_par%icntl(22) = ooc ! out-of-core factorization

!pardiso_par%icntl(11) = 0 ! return statistics
! pardiso_par%icntl(7) = 5  ! ordering type
!pardiso_par%cntl(1) = 0.d0  ! rel. threshold for numerical pivoting


pardiso_factorization: do

   pardiso_par%JOB = 1  !  analysis
   CALL Dpardiso(pardiso_par)

   pardiso_par%JOB = 2  !  factorization
   CALL Dpardiso(pardiso_par)

   ierr = pardiso_par%INFOG(1)


   if ( pardiso_par%INFOG(1) == -9 .or. pardiso_par%INFOG(1) == -8 ) then
      ! Main internal real/complex workarray S too small.
      pardiso_par%icntl(14) = pardiso_par%icntl(14) + 10
      if ( pardiso_par%MYID == 0 ) then
         write(*,30) pardiso_par%icntl(14)
         30 format(/'pardiso percentage increase in the estimated working space',  &
                   /'increased to',i4)
      end if

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

TYPE (Dpardiso_STRUC),intent(inout):: pardiso_par
integer(kind=8),intent(in):: nrhs   ! # of right-hand-sides
real(kind=8),intent(in):: rhs(nrhs * pardiso_par%N)
real(kind=8),intent(out),target:: x(nrhs * pardiso_par%N)  ! solution
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


   pardiso_par%JOB = 3  ! Solve the system.
   CALL Dpardiso(pardiso_par)
   ! At this point pardiso_par%RHS (rhs) contains the solution.


return
end subroutine solve

!--------------------------------------------------------------

subroutine destroy(pardiso_par)

implicit none
TYPE (Dpardiso_STRUC),intent(inout):: pardiso_par

!if (.not. pardiso_initialized)  return

!  Destroy the instance (deallocate internal data structures)
pardiso_par%JOB = -2
CALL Dpardiso(pardiso_par)

if (associated(pardiso_par%A)) then
   deallocate(pardiso_par%IRN, pardiso_par%JCN, pardiso_par%A)
end if

!if (pardiso_par%MYID == 0) then
!   close(ui_out)
!end if

return
end subroutine destroy


end module pardiso_mod
