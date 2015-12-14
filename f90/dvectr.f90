SUBROUTINE dvectr (gpt,x,u,pen)
     
 
 INTEGER, INTENT(IN)                      :: gpt(1)
 REAL, INTENT(IN)                         :: x(3,1)
 REAL, INTENT(IN)                         :: u(2,1)
 INTEGER, INTENT(IN OUT)                  :: pen
 
 
 COMMON /BLANK/ ngp
 
 CALL line (0,0,0,0,0,-1)
 
!     DO NOT DRAW A VECTOR AT ANY GRID POINT WHOSE INDEX .LE. 0.
 
 DO  i = 1,ngp
   j  = gpt(i)
   IF (j <= 0) CYCLE
   x1 = x(2,j)
   y1 = x(3,j)
   x2 = u(1,j)
   y2 = u(2,j)
   CALL line (x1,y1,x2,y2,pen,0)
 END DO
 
 CALL line (0,0,0,0,0,1)
 RETURN
END SUBROUTINE dvectr
