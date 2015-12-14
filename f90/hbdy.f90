SUBROUTINE hbdy (ecpt,necpt,iopt,rvect,ivect)
     
!     THIS SUBROUTINE CALCULATES THE GEOMETRIC PROPERTIES OF THE VARIOUS
!     TYPES OF HBDY ELEMENTS. IOPT IS DESCRIBED BELOW
 
!     THE ECPT INPUT DATA IS
 
!     POSITION    DATA
!        1       EL ID
!        2       FLAG
!        3       SIL-1
!        4       SIL-2
!        5       SIL-3
!        6       SIL-4
!        7       SIL-5
!        8       SIL-6
!        9       SIL-7
!       10       SIL-8
!       11       VECTOR V1
!       12       VECTOR V2
!       13       VECTOR V3
!       14       ECPT14
!       15       MAT ID
!       16       A-FACTOR
!       17       EMISSIVITY
!       18       ABSORBTIVIY
!       19       R1
!       20       R2
!       21       CS-1
!       22       X1
!       23       Y1
!       24       Z1
!       25       CS-2
!       26       X2
!       27       Y2
!       28       Z2
!       29       CS-3
!       30       X3
!       31       Y3
!       32       Z3
!       33       CS-4
!       34       X4
!       35       Y4
!       36       Z4
!       37-52    NOT USED
!       53       AVG. EL. TEMP.
 
!     THE VALUE OF FLAG INDICATES THE TYPE OF ELEMENT
 
!       FLAG     TYPE
!       ****     ****
!        1       POINT
!        2       LINE
!        3       REV
!        4       TRIANGLE
!        5       QUADRILATERAL
!        6       ELLIPTIC CYLINDER
!        7       FTUBE
 
 
!     THE OUTPUT DATA IS PLACED IN  VECT AND IVECT
!         THE FORMATS ARE
 
!     POSITION
!          IOPT=  1             2
!      1        EL ID         EL ID
!      2        AREA          AREA
!      3        EMIS          SIL-1
!      4        ---           SIL-2
!      5        SIL-1         SIL-3
!      6        SIL-2         SIL-4
!      7        SIL-3         AREA-1
!      8        SIL-4         AREA-2
!      9        GFACT-1       AREA-3
!     10        GFACT-2       AREA-4
!     11        GFACT-3       N1X
!     12        GFACT-4       N1Y
!     13                      N1Z
!     14                      N2X  -  FOR FLAG = 6 ONLY
!     15                      N2Y  -
!     16                      N2Z  -
 
 
 
 REAL, INTENT(OUT)                        :: ecpt(36)
 INTEGER, INTENT(IN)                      :: necpt(5)
 INTEGER, INTENT(IN OUT)                  :: iopt
 REAL, INTENT(OUT)                        :: rvect(16)
 INTEGER, INTENT(OUT)                     :: ivect(5)
 INTEGER :: flag
 REAL :: dxyz(3), v(3)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /condas/ consts(5)
 EQUIVALENCE     (consts(1),pi), (dxyz(1),dx), (dxyz(2),dy), (dxyz(3),dz)
 
 
 DO  i = 1,16
   rvect(i) = 0.0
   ivect(i) = 0
 END DO
 ivect(1) = necpt(1)
 flag = necpt(2)
 IF (flag <= 0 .OR. flag > 7) GO TO 210
 IF (flag == 7) ecpt(16) = pi*(ecpt(19) + ecpt(20))
 
 SELECT CASE ( flag )
   CASE (    1)
     GO TO 20
   CASE (    2)
     GO TO 30
   CASE (    3)
     GO TO 40
   CASE (    4)
     GO TO 50
   CASE (    5)
     GO TO 60
   CASE (    6)
     GO TO 90
   CASE (    7)
     GO TO 30
 END SELECT
 
!     FLAG = POINT
 
 20 ivect(3) = necpt(3)
 rvect(7) = ecpt(16)
 rvect(2) = ecpt(16)
 CALL sanorm (*110,ecpt(11))
 npts = 1
 GO TO 110
 
!     FLAG = LINE
 
 30 ivect(3) = necpt(3)
 ivect(4) = necpt(4)
 npts = 2
 dx = ecpt(26) - ecpt(22)
 dy = ecpt(27) - ecpt(23)
 dz = ecpt(28) - ecpt(24)
 
 temp = dx**2 + dy**2 + dz**2
 IF (temp <= 1.0E-20) GO TO 210
 
!     AREA CALCULATIONS
 
 rvect(2) = ecpt(16)*SQRT(temp)
 rvect(7) = rvect(2)*0.5
 rvect(8) = rvect(7)
 
!     NORMAL VECTOR CALCULATIONS
 
 temp  =(dx*ecpt(11) + dy*ecpt(12) + dz*ecpt(13))/temp
 rvect(11) = ecpt(11) - temp*dx
 rvect(12) = ecpt(12) - temp*dy
 rvect(13) = ecpt(13) - temp*dz
 
!     NORMALIZE
 
 CALL sanorm (*110,rvect(11))
 GO TO 110
 
!     TYPE= REV
 
 40 ivect(3) = necpt(3)
 ivect(4) = necpt(4)
 npts = 2
 dx = ecpt(26) - ecpt(22)
 dz = ecpt(28) - ecpt(24)
 temp = SQRT(dx**2 +dz**2)*pi
 IF (temp <= 1.0E-20) GO TO 210
 rvect(7) = (2.0*ecpt(22) + ecpt(26))*temp/3.0
 rvect(8) = (2.0*ecpt(26) + ecpt(22))*temp/3.0
 rvect(2) =  rvect(7) + rvect(8)
 
 temp = temp/pi
 rvect(11) = dz/temp
 rvect(13) =-dx/temp
 GO TO 110
 
!     FLAG = AREA3
 
 50 ivect(3) = necpt(3)
 ivect(4) = necpt(4)
 ivect(5) = necpt(5)
 npts = 3
 dx = ecpt(26) - ecpt(22)
 dy = ecpt(27) - ecpt(23)
 dz = ecpt(28) - ecpt(24)
 rvect(7) = ecpt(30) - ecpt(26)
 rvect(8) = ecpt(31) - ecpt(27)
 rvect(9) = ecpt(32) - ecpt(28)
 
!     CALC. NORMAL VECTOR
 
 CALL saxb (dxyz,rvect(7),rvect(11))
 
 CALL sanorm (*210,rvect(11))
 
 rvect(2) = temp/2.0
 rvect(7) = temp/6.0
 rvect(8) = rvect(7)
 rvect(9) = rvect(7)
 
 GO TO 110
 
!     FLAG = AREA4
 
 60 DO  i = 3,6
   ivect(i) = necpt(i)
 END DO
 npts = 4
 DO  i = 1,3
   
!     CALCULATE  DIFFERENCE VECTORS
   
!        R2 - R1
   
   rvect(i+6) = ecpt(i+25) - ecpt(i+21)
   
!        R3 - R1
   
   rvect(i+13) = ecpt(i+29) - ecpt(i+21)
   
!        R4 - R2
   
   v(i) = ecpt(i+33) - ecpt(i+25)
 END DO
 
!        (R3 - R1) X (R4 - R2)
 
 CALL saxb (rvect(14),v,rvect(11))
 
!     2*AREA
 
 temp  = SQRT(rvect(11)**2 + rvect(12)**2 + rvect(13)**2)
 rvect(2) = temp/2.0
 
!     NORMALIZE
 
 CALL sanorm (*210,rvect(11))
 
 CALL saxb (rvect(7),rvect(14),dxyz)
 
!     AREA OF TRIANGLE 123
 
 temp = SQRT(dx**2 + dy**2 + dz**2)/2.0
 
 CALL saxb (rvect(7),v,dxyz)
 
!     AREA OF TRIANGLE 412
 
 dx =  SQRT(dx**2 + dy**2 + dz**2)/2.0
 
!     AREA FOR POINTS
 
 rvect( 7) = (rvect(2)+dx   )/6.0
 rvect( 8) = (rvect(2)+temp )/6.0
 rvect( 9) = (rvect(2)*2.-dx)/6.0
 rvect(10) = (rvect(2)*2.-temp)/6.0
 rvect(14) = 0.0
 rvect(15) = 0.0
 rvect(16) = 0.0
 npts = 4
 GO TO 110
 
!     FLAG = ELCYL
 
 90 ivect(3) = necpt(3)
 ivect(4) = necpt(4)
 npts = 2
 dx = ecpt(26) - ecpt(22)
 dy = ecpt(27) - ecpt(23)
 dz = ecpt(28) - ecpt(24)
 temp = SQRT(dx**2 + dy**2 + dz**2)
 rvect(2) = temp*ecpt(16)
 IF (iopt == 3) rvect(2) = temp
 IF (temp <= 0) GO TO 210
 CALL saxb (ecpt(11),dxyz,rvect(14))
 CALL saxb (dxyz,rvect(14),rvect(11))
 
 CALL sanorm (*210,rvect(11))
 CALL sanorm (*210,rvect(14))
 DO  i = 1,3
   rvect(i+10) = rvect(i+10)*ecpt(20)
   rvect(i+13) = rvect(i+13)*ecpt(19)
 END DO
 rvect(7) = rvect(2)/2.0
 rvect(8) = rvect(7)
 
!     IOPT EQUALS 1
!     CALCULATE G FACTORS. STORE IN NEW LOCATIONS.
!     WORK FROM LAST TO FIRST
 
!     CHECK FOR ZERO AREA
 
 110 area = rvect(2)
 IF (area < 1.0E-20) GO TO 210
 120 IF (iopt >       1) GO TO 170
 DO  i = 1,npts
   j = npts - i + 1
   rvect(j+8) =  rvect(j+6)/area
 END DO
 
 DO  i = 1,4
   j =  5-i
   IF (j -npts > 0) THEN
     GO TO   140
   ELSE
     GO TO   150
   END IF
   140 ivect(j+4) = 0
   CYCLE
   150 ivect(j+4) = ivect(j+2)
 END DO
 
!     STORE EMISSIVITY VALUE
 
 rvect(3) = ecpt(17)
 RETURN
 
!     IOPT EQUALS 2
 
 170 IF (iopt == 2) RETURN
 DO  i = 1,npts
   rvect(i+6) = rvect(i+6)*ecpt(18)
 END DO
 RETURN
 
 210 WRITE  (6,220) uwm,necpt(1)
 220 FORMAT (a25,' 2154, ZERO AREA OR ILLEGAL CONNECTION FOR HBDY ',  &
     'ELEMENT NUMBER',i9)
 area = 1.0
 GO TO 120
END SUBROUTINE hbdy
