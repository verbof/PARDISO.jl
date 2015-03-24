

!-------------------------------------------------

function factor_pardiso_cmplx( n, sym, ooc, A, jA, iA, ierr )  result(pm_out)
!DIR$ ATTRIBUTES DLLEXPORT :: factor_pardiso_cmplx
!DIR$ ATTRIBUTES ALIAS: 'factor_pardiso_cmplx_':: factor_pardiso_cmplx

use pardiso_cmplx_mod, only: init, convert_to_pardiso_format, factor_matrix, destroy
implicit none

INCLUDE 'zpardiso_struc.h'

integer(kind=8):: pm_out   ! pardiso pointer
integer(kind=8),intent(in):: n  ! # of rows in A
integer(kind=8),intent(in):: sym  ! =0 unsymmetric, =1 symm. pos def, =2 general symm.
integer(kind=8),intent(in):: ooc  ! = 0 in-core factorization, = 1 out-of-core factorization
complex(kind=8),intent(in):: A(*)
integer(kind=8),intent(in):: jA(*), iA(n+1)

integer(kind=8),intent(out):: ierr
! integer,intent(out):: ierr  ! =0 no error, < 0 error (memory, singular, etc.)

TYPE(Zpardiso_STRUC),pointer:: ppardiso_par
TYPE(Zpardiso_STRUC):: pardiso_par

pointer ( pm, pardiso_par )

allocate(ppardiso_par)
pm = loc(ppardiso_par)

call init( pardiso_par, sym)
call convert_to_pardiso_format(pardiso_par, n, A,jA,iA, ierr )
if (ierr < 0) return  ! allocation error

call factor_matrix(pardiso_par, ooc, ierr)

if ( ierr < 0 ) then
   ! Error in factorization.
   call destroy(pardiso_par)
   pm_out = 0
else
   pm_out = pm
end if

return
end  function factor_pardiso_cmplx

!-------------------------------------------------

subroutine solve_pardiso_cmplx( pm_in, nrhs, rhs, x, transpose )
! Solve A*x = rhs

!DIR$ ATTRIBUTES DLLEXPORT :: solve_pardiso_cmplx
!DIR$ ATTRIBUTES ALIAS: 'solve_pardiso_cmplx_':: solve_pardiso_cmplx
use pardiso_cmplx_mod, only: solve

implicit none

INCLUDE 'zpardiso_struc.h'

integer(kind=8),intent(in):: pm_in  ! pardiso pointer
integer(kind=8),intent(in):: nrhs  ! # of right-hand-sides
complex(kind=8),intent(in):: rhs(*)   ! right-hand-side
complex(kind=8),intent(out):: x(*)    ! solution
integer(kind=8),intent(in):: transpose   ! =1 for transpose

TYPE(Zpardiso_STRUC):: pardiso_par
pointer ( pm, pardiso_par )
pm = pm_in

call solve(pardiso_par, nrhs, rhs, x, (transpose==1) )

return
end subroutine solve_pardiso_cmplx

!-------------------------------------------------

subroutine destroy_pardiso_cmplx( pm_in )
!  Destroy the instance (deallocate internal data structures)

!DIR$ ATTRIBUTES DLLEXPORT :: destroy_pardiso_cmplx
!DIR$ ATTRIBUTES ALIAS: 'destroy_pardiso_cmplx_':: destroy_pardiso_cmplx
use pardiso_cmplx_mod, only: destroy

implicit none
INCLUDE 'zpardiso_struc.h'

integer(kind=8),intent(in):: pm_in  ! pardiso pointer
TYPE(Zpardiso_STRUC):: pardiso_par
pointer ( pm, pardiso_par )
pm = pm_in

call destroy(pardiso_par)

return
end subroutine destroy_pardiso_cmplx

!-------------------------------------------------

