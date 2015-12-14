SUBROUTINE fbsinv (x,y,iobuff)
     
!     SINGLE PRECISION VERSION
 
!     FBSINV IS A SPECIAL FORWARD-BACKWARD SUBSTITUTION ROUTINE FOR
!     INVPWR. IT OPERATES ON CONJUNCTION WITH SDCOMP.
!     THE ARITHMETIC PRECISION IS THAT OF THE INPUT FILE
 
!     FILEL  = MATRIX CONTROL BLOCK FOR THE LOWER TRIANGLE
!     X      = THE LOAD VECTOR
!     Y      = THE SOLUTION VECTOR
!     IOBUFF = NOT USED
 
 
 REAL, INTENT(IN)                         :: x(1)
 REAL, INTENT(OUT)                        :: y(1)
 INTEGER, INTENT(IN OUT)                  :: iobuff
 INTEGER :: filel   ,parm(3)  ,iblk(15)
 
 COMMON /fbsx  / filel(7)
 EQUIVALENCE    (filel(3),nrow), (filel(5),ltype)
 DATA    parm  / 4H      ,4HFBSI, 4HNV   /
 
!     FORWARD PASS
 
 parm(1) = filel(1)
 iblk(1) = filel(1)
 IF (ltype == 2) GO TO 20
 IF (ltype /= 1) GO TO 50
 
!     TRANSFER THE SINGLE PRECISION LOAD VECTOR TO THE SOLUTION VECTOR
 
 DO  i = 1,nrow
   y(i) = x(i)
 END DO
 CALL fbs1 (iblk,y,y,nrow)
 GO TO 40
 
!     TRANSFER THE DOUBLE PRECISION LOAD VECTOR TO THE SOLUTION VECTOR
 
 20 nrow2 = 2*nrow
 DO  i = 1,nrow2
   y(i) = x(i)
 END DO
 CALL fbs2 (iblk,y,y,nrow2)
 
 40 CALL REWIND (filel)
 CALL skprec (filel,1)
 RETURN
 
!     FATAL ERRORS
 
 50 CALL mesage (-7,parm(1),parm(2))
 RETURN
END SUBROUTINE fbsinv
