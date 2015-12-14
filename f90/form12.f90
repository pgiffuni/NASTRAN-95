SUBROUTINE form12 (u0, udot0, u1, p0, p1, deltt, ibuf)
!*******
!     FORM12 GENERATES THE STARTING VECTORS FOR THE INTEGRATION MODULE
 
!     THIS ROUTINE IS SUITABLE FOR DOUBLE PRECISION OPERATION
!*******
 
 DOUBLE PRECISION, INTENT(IN)             :: u0(1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: udot0(1)
 DOUBLE PRECISION, INTENT(OUT)            :: u1(1)
 DOUBLE PRECISION, INTENT(OUT)            :: p0(1)
 DOUBLE PRECISION, INTENT(OUT)            :: p1(1)
 REAL, INTENT(IN)                         :: deltt
 INTEGER, INTENT(IN OUT)                  :: ibuf(1)
 
 
 
 
 COMMON /BLANK/  dummy(5)  ,istart
 COMMON   /trdxx /  ifilk(7)  ,ifilm(7) ,ifilb(7)
 
 nrow = ifilk(2)
!*******
!     FORM U(-1)
!*******
 DO  i = 1,nrow
   p1(i) = 0.0D0
   u1(i) = u0(i)-deltt*udot0(i)
 END DO
 IF (istart >= 0) GO TO 30
 DO  i = 1, nrow
   p0(i) = 0.0D0
 END DO
!*******
!     FORM P0
!*******
 CALL matvc2(u0(1),p0(1),ifilk(1),ibuf)
 CALL matvc2(udot0(1),p0(1),ifilb(1),ibuf)
!*******
!     FORM P(-1)
!*******
 CALL matvc2(udot0(1),p1(1),ifilk(1),ibuf)
 DO  i = 1,nrow
   p1(i) = p0(i)-deltt*p1(i)
 END DO
 RETURN
 
!     ALTERNATE STARTING METHOD
 
 30 CALL matvc2 (u0(1), p1(1), ifilk(1), ibuf)
 CALL matvc2 (udot0(1), p1(1), ifilb(1), ibuf)
 DO  i = 1, nrow
   p0(i) = 0.5D0*(p0(i) + p1(i))
   udot0(i) = - udot0(i)*deltt
 END DO
 
!     ADD UDOT CONTRIBUTION
 
 CALL matvc2 (udot0(1), p1(1), ifilk(1), ibuf)
 
!     RESTORE UDOT
 
 DO  i = 1, nrow
   udot0(i) = - udot0(i)/deltt
 END DO
 RETURN
END SUBROUTINE form12
