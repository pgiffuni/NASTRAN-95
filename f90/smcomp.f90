SUBROUTINE smcomp (*,zi,zr,zd)
     
! DRIVER PROGRAM FOR SYMMETRIC DECOMPOSITION.  SUBROUTINE SMCPH1 READS
! THE INPUT MATRIX AND STORES THE DATA EITHER IN MEMORY OR ON THE
! SPILL FILE.  SUBROUTINE SMCPH2 IS THEN CALLED TO PERFORM THE
! MATRIX DECOMPOSITION.
 
INTEGER, INTENT(IN OUT)                  :: zi(4)
REAL, INTENT(IN OUT)                     :: zr(4)
DOUBLE PRECISION, INTENT(IN OUT)         :: zd(4)
 
INTEGER :: module(5), begn, END
INCLUDE          'SMCOMX.COM'

!  mcb   - matrix control block for input matrix
!  lll   - matrix control block for lower triangular matrix
!  dbc   - dbc(1) = available scratch file, dbc(2-7) are not used
!  scr1, scr2, scr3 - three available scratch files
!  lcore - amount of open core available for use
!  ddr   - d.p. values of (real, imaginary) for scaled value of determinant
!  power - scale factor to apply to determinant, determinant=det * 10**power
!  mindd - d.p. value for minimum value of diagonal elements
!  chlsky - cholesky option when =1, i.e., form c matrix

DATA  module / 4HSMCO, 4HMP  , 3*4H    /
DATA  begn   / 4HBEGN /
DATA  END    / 4HEND  /

ierror = 0
ncol   = mcb(2)
module( 3 ) = begn
sturm  = 0
CALL conmsg ( module, 5, 0 )
CALL smcph1 ( zi, zr, zd )
IF ( ierror == 1 ) GO TO 701
IF ( ierror /= 0 ) GO TO 700
CALL smcph2 ( zi, zr, zd )
IF ( ierror == 1 ) GO TO 701

! print roots information if this is an eigenvalue problem, and keep
! two largest shift point data if several shift point movings are involved.

IF ( shftpt > 0. ) WRITE ( nout, 901 ) sturm, shftpt
901   FORMAT( 20X, i5, ' ROOTS BELOW ', 1P,e14.6 )
IF ( sturm /= 0 ) GO TO 100
IF ( KEEP  <= 0 ) GO TO 700
sturm  = KEEP
shftpt = ptshft
GO TO 700
100   IF ( KEEP > sturm ) GO TO 700
jj     = KEEP
rs     = ptshft
KEEP   = sturm
ptshft = jj
shftpt = rs
700   module( 3 ) = END
CALL conmsg ( module, 5, 0 )
IF ( ierror /= 0 ) RETURN 1
GO TO 777
701   CONTINUE
module( 3 ) = END
CALL conmsg ( module, 5, 0 )
777   CONTINUE

RETURN
END SUBROUTINE smcomp
