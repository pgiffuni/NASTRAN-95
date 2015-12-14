SUBROUTINE x trn  y1(x,y,alpha1)
!     SUBROUTINE X TRNS Y (X,Y,ALPHA)
!*******
!     X TRNS Y  FORMS THE DOT PRODUCT X TRANSPOSE * Y = ALPHA
!*******
!     DOUBLE PRECISION   X(1)      ,Y(1)     ,ALPHA
 
 REAL, INTENT(IN)                         :: x(1)
 REAL, INTENT(IN)                         :: y(1)
 DOUBLE PRECISION, INTENT(OUT)            :: alpha1
 
 
 COMMON   /invpwx/  aaa       ,ncol
 
 alpha = 0.0
 DO  i=1,ncol
   alpha = alpha + x(i)*y(i)
 END DO
 alpha1 = alpha
 RETURN
END SUBROUTINE x trn  y1
