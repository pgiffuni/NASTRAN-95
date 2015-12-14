SUBROUTINE dloop (x,y,mpy,END)
!*******
!     DLOOP IMPROVES THE EFFICIENCY OF AN INNER DCOMP LOOP
!*******
 
 DOUBLE PRECISION, INTENT(OUT)            :: x(1)
 DOUBLE PRECISION, INTENT(IN)             :: y(1)
 DOUBLE PRECISION, INTENT(IN)             :: mpy
INTEGER, INTENT(IN)                      :: END

!*******
!     DDLOOP IMPROVES THE EFFICIENCY OF THE ACTIVE ROW LOOP
!*******
DOUBLE PRECISION :: a         ,b(1)     ,c(1)
DOUBLE PRECISION :: xx(1),yy(1)


DO  i = 1,END
x(i) = x(i)+mpy*y(i)
END DO
RETURN
!*******************************
ENTRY ddloop (a,b,c,endd)
DO  i = 1,endd
  a = a-b(i)*c(i)
END DO
RETURN
!************
!     ENTRY FOR ANOTHER LOOP
!***********
ENTRY xloop(xx,yy,nn)
DO  i = 1,nn
  xx(i) = yy(i)
END DO
RETURN
END SUBROUTINE dloop
