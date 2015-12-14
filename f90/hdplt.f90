SUBROUTINE hdplt (x1,y1,ij,im)
     
!     PLOTS POINTS GOVERNED BY THE VALUE OF IM.
 
!     NOTE THAT CALL PLOT(X,Y,2) MEANS MOVE PEN FROM THE CURRENT
!     POSITION TO THE POINT,(X,Y),WITH THE PEN DOWN.
 
!     CALL PLOT(X,Y,3) MEANS MOVE THE PEN FROM THE CURRENT POSITION
!     TO THE POINT,(X,Y), WITH THE PEN UP.
 
 
 REAL, INTENT(IN)                         :: x1(4)
 REAL, INTENT(IN)                         :: y1(4)
 INTEGER, INTENT(OUT)                     :: ij
 INTEGER, INTENT(OUT)                     :: im
 LOGICAL :: debug
 INTEGER :: ppen
 
 COMMON /drwdat/ dum(3),ppen
 COMMON /system/ ibuf,nout
 DATA    debug / .false./
 
 IF (debug) WRITE (nout,1000)ij,im,(x1(i),i=1,4),(y1(j),j=1,4)
 1000 FORMAT (7H hdplt ,2I3,8F12.5)
 IF (im == 1) GO TO 20
 xvalue = (x1(2))/x1(4)
 yvalue = (y1(2))/y1(4)
 IF (ij == 0) GO TO 10
 CALL line (xold,yold,xvalue,yvalue,ppen,0)
 xold = xvalue
 yold = yvalue
 GO TO 30
 10 CONTINUE
 xold = xvalue
 yold = yvalue
 ij = 1
 GO TO 30
 20 CONTINUE
 ij = 0
 30 CONTINUE
 RETURN
END SUBROUTINE hdplt
