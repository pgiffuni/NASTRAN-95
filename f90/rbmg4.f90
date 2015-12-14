SUBROUTINE rbmg4
!*****
! RBMG4 COMPUTES MR FROM THE MATRIX EQUATION
!      MR = MRR + DM(T) * MLR + MLR(T) * DM + DM(T) * MLL * DM
!*****
 INTEGER :: scr1,scr2,dm
 INTEGER :: scr3
!*****
!     INPUT DATA FILES
!*****
 DATA dm,mll,mlr,mrr/101,102,103,104/
!*****
!     OUTPUT DATA FILES
!*****
 DATA mr/201/
!*****
!     SCRATCH DATA FILES
!*****
 DATA scr1,scr2,scr3/301,302,303/
!*****
!     COMPUTE MR
!*****
 CALL elim(mrr,mlr,mll,dm,mr,scr1,scr2,scr3)
 RETURN
END SUBROUTINE rbmg4
