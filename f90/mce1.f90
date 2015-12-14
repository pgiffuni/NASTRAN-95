SUBROUTINE mce1
     
!     MCE1 PARTITIONS RG INTO RM AND RN
!     THEN SOLVES THE MATRIX EQUATION RM * GM = -RN.
 
 
 INTEGER :: uset   ,rg    ,gm    ,scr1  ,scr2  ,scr3  ,rm   , rn     ,l     ,u
 COMMON /BLANK/ uset   ,rg    ,gm    ,scr1  ,scr2  ,scr3  ,rm   ,  &
     rn     ,l     ,u     ,mcb(7)
 
!     SET INPUT, OUTPUT AND SCRATCH FILES
 
 uset = 101
 rg   = 102
 gm   = 201
 scr1 = 304
 scr2 = 305
 scr3 = 301
 rm   = 302
 rn   = 303
 l    = 306
 u    = 307
 
!     PARTITION RG INTO RM AND RN
 
 CALL mce1a
 
!     TEST FOR RM DIAGONAL
 
 mcb(1) = rm
 CALL rdtrl (mcb)
 IF (mcb(5) == 1 .AND. mcb(6) == 1) GO TO 50
 IF (mcb(5) == 2 .AND. mcb(6) == 2) GO TO 50
 
!     RM IS NOT DIAGONAL, DECOMPOSE RM THEN SOLVE FOR GM
!     BY FORWARD-BACKWARD SUBSTITUTION.
 
 CALL mce1b
 CALL mce1c
 RETURN
 
!     RM IS DIAGONAL, COMPUTE GM = -RM(-1) * RN
 
 50 CALL mce1d
 RETURN
END SUBROUTINE mce1
