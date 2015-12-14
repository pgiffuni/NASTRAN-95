SUBROUTINE bmgtns(cstm,ncstm,ecpt,ta)
!/// THIS ROUTINE WAS LIFTED FROM PRETRD AND TRANSD AND CONVERTED
!///// TO HAVE ONE ENTRY POINT
 
! PRETRD SETS UP EVENTUAL CALLS TO TRANSD.  FOR A MODULE TO USE TRANSD
! A CALL TO PRETRD MUST BE INITIATED BY THE MODULE DRIVER ONCE AND ONLY
! ONCE.  CSTM IS ARRAY OF COORDINATE SYSTEM TRANSFORMATION MATRICES
! AND NCSTM IS THE LENGTH OF THIS ARRAY.
 
 
 REAL, INTENT(IN)                         :: cstm(1)
 INTEGER, INTENT(IN)                      :: ncstm
 REAL, INTENT(IN)                         :: ecpt(4)
 DOUBLE PRECISION, INTENT(OUT)            :: ta(9)
 
!*****
! GIVEN THE ECPT ARRAY OF LENGTH 4, THE FIRST WORD BEING AN INTEGER
! COORDINATE SYSTEM IDENTIFICATION NUMBER AND THE NEXT WORDS BEING THE
! REAL COORDINATES OF A POINT IN BASIC COORDINATES, THIS ROUTINE
! COMPUTES THE TRANSFORMATION (DIRECTION COSINE) MATRIX TA WHICH WILL
! MAP A VECTOR FROM THE LOCAL SYSTEM LABELED ECPT(1) TO BASIC COORDI-
! NATES.  TA IS A DOUBLE PRECISION MATRIX.
!*****
 
 
 DOUBLE PRECISION ::  tl(9)  &
     ,                  xn(3)              ,x  &
     ,                  y                  ,z  &
     ,                  r                  ,ke(9) ,                  xl
 
 EQUIVALENCE        (fl1,int1)         ,(fl2,int2)
 
 fl1 = ecpt(1)
 IF(int1 == 0) GO TO 13
 DO  i = 1,ncstm,14
   fl2 = cstm(i)
   IF(int1 /= int2) CYCLE
   kk = i
   fl2 = cstm(i + 1)
   SELECT CASE ( int2 )
     CASE (    1)
       GO TO 7
     CASE (    2)
       GO TO 10
     CASE (    3)
       GO TO 10
   END SELECT
 END DO
 
! THE COORDINATE SYSTEM ID. COULD NOT BE FOUND IN THE CSTM.
 
 CALL mesage (-30,25,int1)
 
! THE COORDINATE SYSTEM IS RECTANGULAR.
 
 7 DO  j = 1,9
   k = kk + 4 + j
   ta(j) = cstm(k)
 END DO
 RETURN
 10 xn(1) = ecpt(2) - cstm(kk+2)
 xn(2) = ecpt(3) - cstm(kk+3)
 xn(3) = ecpt(4) - cstm(kk + 4)
 x = cstm(kk+5)*xn(1)+cstm(kk+8)*xn(2)+cstm(kk+11)*xn(3)
 y = cstm(kk+6)*xn(1)+cstm(kk+9)*xn(2)+cstm(kk+12)*xn(3)
 z = cstm(kk+7)*xn(1)+cstm(kk+10)*xn(2)+cstm(kk+13)*xn(3)
 r = DSQRT(x**2+y**2)
 IF (r == 0.0D0) GO TO 7
 DO  j=1,9
   k=kk+4+j
   ke(j)=cstm(k)
 END DO
 SELECT CASE ( int2 )
   CASE (    1)
     GO TO 11
   CASE (    2)
     GO TO 11
   CASE (    3)
     GO TO 12
 END SELECT
 
! THE COORDINATE SYSTEM IS CYLINDRICAL.
 
 11 tl(1)=x/r
 tl(2)=-y/r
 tl(3)=0.0D0
 tl(4)=-tl(2)
 tl(5)=tl(1)
 tl(6)=0.0D0
 tl(7)=0.0D0
 tl(8)=0.0D0
 tl(9)=1.0D0
 GO TO 125
 
! THE COORDINATE SYSTEM IS SPHERICAL.
 
 12 xl=DSQRT(x*x+y*y+z*z)
 tl(1)=x/xl
 tl(2)=(x*z)/(r*xl)
 tl(3)=-y/r
 tl(4)=y/xl
 tl(5)=(y*z)/(r*xl)
 tl(6)=x/r
 tl(7)=z/xl
 tl(8)=-r/xl
 tl(9)=0.0D0
 125 CALL gmmatd (ke(1),3,3,0, tl(1),3,3,0, ta(1))
 RETURN
 
! THE LOCAL SYSTEM IS BASIC.
 
 13 DO  i=1,9
   ta(i)=0.0D0
 END DO
 ta(1)=1.0D0
 ta(5)=1.0D0
 ta(9)=1.0D0
 RETURN
END SUBROUTINE bmgtns
