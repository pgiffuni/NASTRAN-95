SUBROUTINE cx_trn_y (x,y,alpha)
!*******
!     CX TRN Y FORMS THE DOT PRODUCT X TRANSPOSE * Y = ALPHA WHERE
!     X AND Y ARE COMPLEX
!*******
 
 DOUBLE PRECISION, INTENT(IN)             :: x(1)
 DOUBLE PRECISION, INTENT(IN)             :: y(1)
 DOUBLE PRECISION, INTENT(OUT)            :: alpha(2)
 COMMON   /cinvpx/  aaa       ,ncol
 
 
 ncol2 = ncol+ncol
 alpha(1) = 0.d0
 alpha(2) = 0.d0
 
 DO  i = 1,ncol2,2
   alpha(1) = alpha(1)+x(i)*y(i)-x(i+1)*y(i+1)
   alpha(2) = alpha(2)+x(i)*y(i+1)+x(i+1)*y(i)
 END DO
 
 RETURN
END SUBROUTINE cx_trn_y
