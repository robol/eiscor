#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_rot2_vec2gen
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the generator for a Givens rotation represented by 2
! real numbers: A strictly real cosine and a scrictly real sine.  The first
! column is constructed to be parallel with the real vector [A,B]^T. 
!
! The rot2 refers to a rotation desribed by two double and the vec2 to a vector
! of lenght two described by two double.
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  A               REAL(8) 
!                    the first component of the real vector [A,B]^T
!
!  B               REAL(8) 
!                    the second component of the real vector [A,B]^T
!
! OUTPUT VARIABLES:
!
!  C               REAL(8)
!                    on exit contains the generator for the cosine
!                    component of the Givens rotation
!
!  S               REAL(8)
!                    on exit contains the generator for the sine
!                    component of the Givens rotation
!
!  NRM             REAL(8)
!                    on exit contains the norm of the vector [A,B]^T
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! EXCEPTIONAL CASES
!
!    A    |    B    |    C    |    S    |   NRM 
! ------- | ------- | ------- | ------- | -------
!     0d0 !     0d0 !     1d0 !     0d0 |   0d0
! +/- INF | +/- XdX | +/- 1d0 |     0d0 |   INF
! +/- XdX | +/- INF |     0d0 | +/- 1d0 |   INF
! +/- INF | +/- INF |     NAN |     NAN |   NAN
!     NAN | +/- XdX |     NAN |     NAN |   NAN
! +/- XdX |     NAN |     NAN |     NAN |   NAN
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_rot2_vec2gen(A,B,C,S,NRM)

  implicit none
  
  ! input variables
  real(8), intent(in) :: A,B
  real(8), intent(inout) :: C,S,NRM
  
  ! construct rotation
  if (B.EQ.0d0) then
     if (A.NE.A) then
        C = A
        S = A
        NRM =A
     else if (A.LT.0d0) then
        C = -1d0
        S = 0d0
        NRM = -A
     else
        C = 1d0
        S = 0d0
        NRM = A
     end if
  else if (abs(A**2 + B**2 - 1d0) .LT. 3d0*EISCOR_DBL_EPS) then
     C = A
     S = B
     NRM = 1d0
  else if (abs(A).GE.abs(B)) then
     S = B/A
     if (A.LT.0) then
        NRM = -sqrt(1.d0 + S**2)
     else
        NRM = sqrt(1.d0 + S**2)
     end if
     C =  1.d0/NRM
     S =  S*C
     NRM =  A*NRM
  else
     C = A/B;
     if (B.LT.0) then
        NRM = -sqrt(1.d0 + C**2)
     else
        NRM = sqrt(1.d0 + C**2)
     end if
     S =  1.d0/NRM
     C =  C*S
     NRM =  B*NRM
  end if
           
end subroutine d_rot2_vec2gen
