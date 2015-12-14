SUBROUTINE presax (iharm)
     
!     THIS ROUTINE APPLIES PRESSURE LOADS TO AXISYMMETRIC SHELL
 
 
 INTEGER, INTENT(IN OUT)                  :: iharm
 LOGICAL :: piez
 INTEGER :: FILE,slt,icard(6),iord(2),NAME(2)
 REAL :: card(6),gpco(4,2)
 COMMON /condas/ pi,twopi,radeg,degrad,s4pisq
 COMMON /zzzzzz/ z(1)
 COMMON /loadx / lc,slt,bgpdt,OLD
 COMMON /system/ ksystm(80)
 EQUIVALENCE     (icard(1),card(1))
 DATA    NAME  / 4HPRES,4HAX  /
 
!     DEFINITION OF VARIABLES
 
!     N        NUMBER OF CURRENT HARMONIC
!     FILE     FILE NAME FOR ERROR MESAGES
!     SLT      STATIC LOADS TABLE
!     CARD     CARD IMAGE OF PRESAX CARD
!     DEGRAD   CONVERSION FACTOR FOR DEGREES TO RADIANS
!     IORD     ARRAY GIVING OPTIMUM ORDER FOR LOOKING UP POINTS IN BGPDT
!     OLD      CURRENT POSITION OF BGPDT
!     GPCO     ARRAY HOLDING BGPDT DATA FOR EACH RING
!     XL       DISTANCE  BETWEEN  RINGS
!     SINSI    SIN  ANGLE BETWEEN RINGS
!     COSSI    COS  ANGLE BETWEEN RINGS
!     ISILA    SIL VALUE  OF CURRENT HARMONIC - RING A
!     ISILB    SIL VALUE  OF CURRENT HARMONIC - RING B
!     IHARM    SUBCASE  INDICATOR  1 = SINE  2 = COSINE
 
 
!     BRING IN PRESAX CARD
 
 FILE = slt
 CALL READ (*910,*920,slt,card(1),6,0,iflag)
 n  = icard(6) + 1
 xi = n - 1
 
!     CONVERT PHI1,PHI2 TO RADIANS
 
 card(4) = card(4)*degrad
 card(5) = card(5)*degrad
 
!     PICK UP BGPDT DATA FOR RINGS
 
!     IF 1ST. RING IS NEGATIVE, THIS IS A SURFACE CHARGE LOAD IN A
!     PIEZOELECTRIC PROBLEM
 
 piez = .false.
 IF (ksystm(78) /= 1 .OR. icard(2) > 0) GO TO 5
 piez = .true.
 icard(2) = -icard(2)
 5 CONTINUE
 CALL permut (icard(2),iord(1),2,OLD)
 DO  i = 1,2
   j = iord(i) + 1
   CALL fndpnt (gpco(1,j-1),icard(j) )
 END DO
 xl = SQRT((gpco(2,2) - gpco(2,1))**2 + (gpco(3,2) - gpco(3,1))**2)
 IF (xl == 0.0) CALL mesage (-30,26,-1)
 sinsi = (gpco(2,2) - gpco(2,1))/xl
 cossi = (gpco(3,2) - gpco(3,1))/xl
 CALL fndsil (icard(2))
 isila = icard(2)
 CALL fndsil (icard(3))
 isilb = icard(3)
 
!     APPLY LOADS TO ALL HARMONICS
 
 IF (n /= 1) GO TO 20
 
!     APPLY LOADS TO ZERO HARMONIC - COSINE SUBCASE ONLY
 
 IF (iharm /= 2) GO TO 90
 pr = (card(5) - card(4))
 GO TO 30
 
!     I .GT. 1  APPLY  SINE AND COSINE FACTORS
 
 20 IF (iharm == 1) GO TO 40
 
!     COSINE CASE
 
 pr = (SIN(xi*card(5)) - SIN(xi*card(4)))/xi
 GO TO 30
 
!     SINE CASE
 
 40 pr = -(COS(xi*card(5)) - COS(xi*card(4)))/xi
 
!     APPLY LOADS
 
 30 pr = pr*card(1)*xl
 prpiez = pr
 prc = pr*cossi
 prs =-pr*sinsi
 pr  = gpco(2,1)/3.0 + gpco(2,2)/6.0
 IF (.NOT.piez) GO TO 35
 
!     PIEZOELECTRIC
 
 prc = 0.
 prs = 0.
 35 CONTINUE
 z(isila  ) = z(isila  ) + prc*pr
 z(isila+2) = z(isila+2) + prs*pr
 IF (piez)  z(isila+3) = z(isila+3) + prpiez*pr
 pr = gpco(2,2)/3.0 + gpco(2,1)/6.0
 z(isilb  ) = z(isilb  ) + prc*pr
 z(isilb+2) = z(isilb+2) + prs*pr
 IF (piez)  z(isilb+3) = z(isilb+3) + prpiez*pr
 90 RETURN
 
!     FILE ERRORS
 
 910 ip1 = -2
 911 CALL mesage (ip1,FILE,NAME(1))
 920 ip1 = -3
 GO TO 911
END SUBROUTINE presax
