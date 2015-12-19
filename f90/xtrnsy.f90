SUBROUTINE xtrnsy (x,y,alpha)
!*******
!     X TRNS Y  FORMS THE DOT PRODUCT X TRANSPOSE * Y = ALPHA
!*******
 
 DOUBLE PRECISION, INTENT(IN)             :: x(1)
 DOUBLE PRECISION, INTENT(IN)             :: y(1)
 DOUBLE PRECISION, INTENT(OUT)            :: alpha
 
 COMMON   /invpwx/  aaa       ,ncol
 
 alpha = 0.d0

 DO  i=1,ncol
   alpha = alpha + x(i)*y(i)
 END DO
 
 RETURN
END SUBROUTINE xtrnsy
