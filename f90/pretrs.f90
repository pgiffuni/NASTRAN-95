SUBROUTINE pretrs (cstmx,ncstmx)
     
!     PRETRS SETS UP EVENTUAL CALLS TO TRANSS.  FOR A MODULE TO USE
!     TRANSS A CALL TO PRETRS MUST BE INITIATED BY THE MODULE DRIVER
!     ONCE AND ONLY ONCE.  CSTMX IS ARRAY OF COORDINATE SYSTEM
!     TRANSFORMATION MATRICES AND MCSTMX IS THE LENGTH OF THIS ARRAY.
 
!     THE CSTMX ARRAY MUST BE WITHIN OPEN CORE BOUND, AND THERE IS NO
!     CHECK ON THIS CONDITION
 
!     GIVEN THE ARRAY ECPT OF LENGTH 4, THE FIRST WORD BEING AN INTEGER
!     COORDINATE SYSTEM IDENTIFICATION NUMBER AND THE NEXT WORDS BEING
!     THE REAL COORDINATES OF A POINT IN BASIC COORDINATES, THIS ROUTINE
!     COMPUTES THE TRANSFORMATION (DIRECTION COSINE) MATRIX TA WHICH
!     WILL MAP A VECTOR FROM THE LOCAL SYSTEM LABELED ECPT(1) TO BASIC
!     COORDINATES.
 
!     REVISED  7/92 BY G.CHAN/UNISYS. NEW REFERENCE TO CSTM ARRAY SUCH
!     THAT THE SOURCE CODE IS UP TO ANSI FORTRAN 77 STANDARD.
 
 
 REAL, INTENT(IN OUT)                     :: cstmx(1)
 INTEGER, INTENT(IN)                      :: ncstmx
 INTEGER :: offset
 REAL :: ke
 DIMENSION  ecpt(4),ta(9),tl(9),ke(9),xn(3)
 COMMON /zzzzzz/ cstm(1)
 EQUIVALENCE     (fl1,int1),(fl2,int2)
 
 ncstm  = ncstmx
 offset = locfx(cstmx(1)) - locfx(cstm(1))
 IF (offset < 0) CALL errtrc ('pretrs  ',1)
 icheck = 123456789
 RETURN
 
 
 ENTRY transs (ecpt,ta)
!     ======================
 
 fl1 = ecpt(1)
 IF (int1 == 0) GO TO 90
 IF (icheck /= 123456789) CALL errtrc ('pretrs  ',10)
 DO  j = 1,ncstm,14
   i   = j + offset
   fl2 = cstm(i)
   IF (int1 /= int2) CYCLE
   kk = i
   fl2 = cstm(i+1)
   SELECT CASE ( int2 )
     CASE (    1)
       GO TO 20
     CASE (    2)
       GO TO 40
     CASE (    3)
       GO TO 40
   END SELECT
 END DO
 
!     THE COORDINATE SYSTEM ID. COULD NOT BE FOUND IN THE CSTM.
 
 CALL mesage (-30,25,int1)
 
!     THE COORDINATE SYSTEM IS RECTANGULAR.
 
 20 DO  j = 1,9
   k = kk + 4 + j
   ta(j) = cstm(k)
 END DO
 RETURN
 
 40 xn(1) = ecpt(2) - cstm(kk+2)
 xn(2) = ecpt(3) - cstm(kk+3)
 xn(3) = ecpt(4) - cstm(kk+4)
 x = cstm(kk+5)*xn(1) + cstm(kk+ 8)*xn(2) + cstm(kk+11)*xn(3)
 y = cstm(kk+6)*xn(1) + cstm(kk+ 9)*xn(2) + cstm(kk+12)*xn(3)
 z = cstm(kk+7)*xn(1) + cstm(kk+10)*xn(2) + cstm(kk+13)*xn(3)
 r = SQRT(x**2 + y**2)
 IF (r == 0.0) GO TO 20
 DO  j = 1,9
   k = kk + 4 + j
   ke(j) = cstm(k)
 END DO
 SELECT CASE ( int2 )
   CASE (    1)
     GO TO 60
   CASE (    2)
     GO TO 60
   CASE (    3)
     GO TO 70
 END SELECT
 
!     THE COORDINATE SYSTEM IS CYLINDRICAL.
 
 60 tl(1) = x/r
 tl(2) =-y/r
 tl(3) = 0.0
 tl(4) =-tl(2)
 tl(5) = tl(1)
 tl(6) = 0.0
 tl(7) = 0.0
 tl(8) = 0.0
 tl(9) = 1.0
 GO TO 80
 
!     THE COORDINATE SYSTEM IS SPHERICAL.
 
 70 xl = SQRT(x**2 + y**2 + z**2)
 tl(1) = x/xl
 tl(2) = (x*z)/(r*xl)
 tl(3) =-y/r
 tl(4) = y/xl
 tl(5) = (y*z)/(r*xl)
 tl(6) = x/r
 tl(7) = z/xl
 tl(8) =-r/xl
 tl(9) = 0.0
 80 CALL gmmats (ke(1),3,3,0, tl(1),3,3,0, ta(1))
 RETURN
 
!     THE LOCAL SYSTEM IS BASIC.
 
 90 DO  i = 1,9
   ta(i) = 0.0
 END DO
 ta(1) = 1.0
 ta(5) = 1.0
 ta(9) = 1.0
 RETURN
END SUBROUTINE pretrs
