SUBROUTINE gptsym (gplst,x,u,sym,deform)
     
 
 INTEGER, INTENT(IN)                      :: gplst(1)
 REAL, INTENT(IN)                         :: x(3,1)
 REAL, INTENT(IN)                         :: u(2,1)
 INTEGER, INTENT(IN OUT)                  :: sym(2)
 INTEGER, INTENT(IN OUT)                  :: deform
 
 
 COMMON /BLANK/ ngp
 
 CALL symbol (0,0,0,-1)
 
!     IF THE GRID POINT INDEX IS 0 (NOT IN SET) OR NEGATIVE (EXCLUDED),
!     NEVER PUT A SYMBOL AT THAT GRID POINT.
 
 DO  i = 1,ngp
   j  = gplst(i)
   IF (j <= 0) CYCLE
   IF (deform /= 0) GO TO 105
   xx = x(2,j)
   yy = x(3,j)
   GO TO 106
   105 xx = u(1,j)
   yy = u(2,j)
   106 CALL symbol (xx,yy,sym,0)
 END DO
 
 CALL symbol (0,0,0,1)
 RETURN
END SUBROUTINE gptsym
