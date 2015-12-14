SUBROUTINE mpyl (a,b,ncola,nrowa,ncolb,c)
     
!     SINCE CDC FORTRAN 5 IMPOSES NO LONGER EXACT NO. OF DUMMY ARGUMENT
!     LIST FOR SUBROUTINE AND ENTRY POINTS, THIS ROUTINE IS NOW MACHINE
!     INDEPENDENT.
 
 
 REAL, INTENT(IN)                         :: a(ncola,nrowa)
 REAL, INTENT(IN)                         :: b(ncolb,ncola)
 INTEGER, INTENT(IN)                      :: ncola
 INTEGER, INTENT(IN)                      :: nrowa
 INTEGER, INTENT(IN)                      :: ncolb
 REAL, INTENT(OUT)                        :: c(ncolb,nrowa)
 
 DIMENSION d(nrowa,ncola),x(3),y(3),vect(3)
 
!     SIMPLE MATRIX MULTIPLICATION
 
 DO  n= 1,ncolb
   DO  l= 1,nrowa
     c(n,l) = 0.0
     DO  m= 1,ncola
       c(n,l) = c(n,l)+b(n,m)*a(m,l)
     END DO
   END DO
 END DO
 RETURN
 
 ENTRY norm (x,y)
!     ================
 
!     NORMALIZE X VECTOR
 
 y(1) = x(1)*x(1)+x(2)*x(2)+x(3)*x(3)
 IF (y(1) == 0.0)  GO TO 15
 w    = 1./SQRT(y(1))
 x(1) = x(1)*w
 x(2) = x(2)*w
 x(3) = x(3)*w
 15 RETURN
 
 ENTRY cross (x,y,vect)
!     ======================
 
!     CROSS PRODUCT
 
 vect(1) = x(2)*y(3)-x(3)*y(2)
 vect(2) = y(1)*x(3)-x(1)*y(3)
 vect(3) = x(1)*y(2)-y(1)*x(2)
 RETURN
 
 ENTRY mpylt (d,b,ncola,nrowa,ncolb,c)
!     =====================================
 
!     TRANSPOSE MULTIPLY
 
 DO  n= 1,ncolb
   DO  l= 1,nrowa
     c(n,l) = 0.0
     DO  m= 1,ncola
       c(n,l) = c(n,l)+d(l,m)*b(n,m)
     END DO
   END DO
 END DO
 RETURN
END SUBROUTINE mpyl
