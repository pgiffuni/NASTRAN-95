SUBROUTINE pcoord (pen)
     
!     PLOTS A COORDINATE TRIAD AT THE LOWER RIGHT CORNER OF A STRUCTURAL
!     PLOT. THIS ROUTINE IS CALLED ONLY BY DRAW
 
!     WRITTEN BY G.CHAN/UNISYS      10/1990
 
 
 INTEGER, INTENT(IN OUT)                  :: pen
 INTEGER :: fpltit,sym(2)
 COMMON /xxparm/ idumm(215),fpltit
 COMMON /pltdat/ idum20(20),size,idum2(2),cntchr(2)
 COMMON /drwaxs/ g(3,4)
 
 
!     COMPUTE THE ORIGIN, WHICH IS A FUNCTION OF FRAME SIZE, CHARACTER
!     VERTICAL AND HORIZONTAL SCALES, PRESENCE OF PTITLE LINE, AND THE
!     OVERALL SIZE OF THE TRIAD
 
!     ALL THE NUMERIC MULTIPLIERS USED BELOW WERE WORKED OUT WITH FRAME
!     SIZE OF 1023.0. THEY SHOULD BE APPLICABLE TO FRAME SIZE OF 3000.
 
 x2 = 0.0
 y2 = 0.0
 y1 = 0.0
 DO  i = 1,3
   IF (g(2,i) > x2) x2 = g(2,i)
   IF (g(3,i) < y2) y2 = g(3,i)
   IF (g(3,i) > y1) y1 = g(3,i)
 END DO
 de = 1.8*cntchr(1)
 sf = 2.7*cntchr(2)
 IF (fpltit == 0) sf = 1.3*sf
 sf = sf/(y1-y2)
 x1 = size - x2*sf - de
 y1 = -y2*sf
 IF (fpltit /= 0) y1 = y1 + 0.8*cntchr(2)
 ep = 0.0001
 of = -1.
 IF (g(2,1) <= ep .AND. g(2,2) <= ep .AND. g(2,3) <= ep) of = +1.
 IF (of == +1.) x1 = x1 - de
 
!     DRAW THE X-Y-Z COORDINATE TRIAD
!     DRAW A CIRCLE AT THE ORIGIN IF ANY AXIS IS IN LINE WITH VIEWER
 
 sym(1) = 6
 sym(2) = 0
 de = 0.8*cntchr(1)
 of = 1.3*of*cntchr(1)
 DO  i = 1,3
   x2 = g(2,i)*sf + x1
   y2 = g(3,i)*sf + y1
   CALL line (x1,y1,x2,y2,pen,0)
   IF (ABS(g(2,i))+ABS(g(3,i)) >= ep) GO TO 20
   CALL symbol (x1,y1,sym,0)
   CALL tipe (x2+of,y2,1,g(i,4),1,0)
   CYCLE
   20 IF (g(2,i) > 0.0) CALL tipe (x2+de,y2,1,g(i,4),1,0)
   IF (g(2,i) <= 0.0) CALL tipe (x2-de,y2,1,g(i,4),1,0)
 END DO
 
 RETURN
END SUBROUTINE pcoord
