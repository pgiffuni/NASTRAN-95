SUBROUTINE flunam (lu,filnam)
     
!     THIS ROUTINE FORMULATES A FORTRAN LOGICAL UNIT NAME FROM A
!     LOGICAL UNIT NUMBER       =       =       =    ===
 
!     INPUT  LU          e.g. LU = 8
!     OUTPUT FILNAM           FILNAM = 'fort.08'  NOTE - IS .08 NOT .8
 
 
 INTEGER, INTENT(IN)                      :: lu
 CHARACTER (LEN=7), INTENT(OUT)           :: filnam
 CHARACTER (LEN=5) :: fort
 CHARACTER (LEN=5) :: for5
 CHARACTER (LEN=7) :: for7
 
 EQUIVALENCE (for5,for7)
 DATA  fort/ 'fort.' /
 
 j = lu + 100
 WRITE  (for7,10) j
 10 FORMAT (4X,i3)
 for5   = fort
 filnam = for7
 RETURN
END SUBROUTINE flunam
