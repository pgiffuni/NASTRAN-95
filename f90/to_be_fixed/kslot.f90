SUBROUTINE kslot (itype)
     
!     THIS ROUTINE CALCULATES THE STIFFNESS MATRIX TERMS FOR THE
!     CSLOT3 AND CSLOT4 TWO DIMENSIONAL LAPLACE ELEMENTS
 
!     IOPT-  CSLOT3 = 0,  CSLOT4 = 1
 
!     THE ECPT DATA FOR THESE ELEMENTS ARE
 
!     FIELD   CSLOT3                CSLOT4
!       1       ID                  ID
!       2       SIL1                SIL1
!       3       SIL2                SIL2
!       4       SIL3                SIL3
!       5       RHO                 SIL4
!       6       BULK                RHO
!       7       M                   BULK
!       8       N                   M
!       9       CID1                N
!       10      R1                  CID1
!      11       Z1                  R1
!      12       W1                  Z1
!      13       CID2                W1
!      14       R2                  CID2
!      15       Z2                  R2
!      16       W2                  Z2
!      17       CID3                W2
!      18       R3                  CID3
!      19       Z3                  R3
!      20       W3                  Z3
!      21       TEMP                W3
!      22                           CID4
!      23                           R4
!      24                           Z4
!      25                           W4
!      26                           TEMP
 
 
 INTEGER, INTENT(IN OUT)                  :: itype
 LOGICAL :: nogo
 INTEGER :: necpt(100) ,out
 DOUBLE PRECISION :: coef       ,fir        ,fiz        ,  &
     r          ,z          ,rki        , a2         ,kij
 CHARACTER (LEN=23) :: ufm
 COMMON  /xmssg /  ufm
 COMMON  /system/  sysbuf     ,out        ,nogo
 COMMON  /sma1cl/  iopt4      ,k4ggsw     ,npvt
 COMMON  /sma1et/  ecpt(100)
 COMMON  /sma1io/  dum1(10)   ,ifile
 COMMON  /sma1dp/  coef       ,fir(3)     ,fiz(3)      ,  &
     r(3)       ,z(3)       ,rki         ,  &
     a2         ,kij        ,nneg        ,  &
     ip         ,nptj       ,iret        ,  &
     lri        ,lrj        ,lrk         , ipvt
 EQUIVALENCE       (ecpt(1),necpt(1))
 
 IF (itype > 0) GO TO 50
 IF (ecpt(5) == 0.0 .OR. necpt(7) == 0) RETURN
 k = -1
 10 k = k + 1
 IF (2*necpt(8)-k*necpt(7) < 0) THEN
   GO TO    30
 ELSE IF (2*necpt(8)-k*necpt(7) == 0) THEN
   GO TO    20
 ELSE
   GO TO    10
 END IF
 20 necpt(7) = necpt(7)*2
 30 ecpt(7)  = FLOAT(necpt(7))/2.0
 DO  i = 1,20
   ecpt(i+50) =  ecpt(i)
 END DO
 iret = 4
 GO TO 170
 
!     THE CSLOT4 ELEMENT IS CHECKED FOR VALIDITY AND THE DATA ARE
!     REARRANGED TO CONFORM TO THE CSLOT3 FORMAT
 
 50 IF (ecpt(6) == 0.0 .OR. necpt(8) == 0) RETURN
 k = -1
 60 k = k + 1
 IF (2*necpt(9)-k*necpt(8) < 0) THEN
   GO TO    80
 ELSE IF (2*necpt(9)-k*necpt(8) == 0) THEN
   GO TO    70
 ELSE
   GO TO    60
 END IF
 70 necpt(8) = necpt(8)*2
 80 ecpt(8)  = FLOAT(necpt(8))/2.0
 
 nneg = 0
 ip   = 0
 DO  i = 1,4
   IF (npvt == necpt(i+1)) ip = ip + 1
   DO  j = 1,3
     nj   = i + j - 1
     IF (nj > 4) nj = nj - 4
     nptj = 4*(nj-1) + 11
     r(j) = ecpt(nptj  )
     z(j) = ecpt(nptj+1)
   END DO
   coef = (r(2)-r(1))*(z(3)-z(1)) - (r(3)-r(1))*(z(2)-z(1))
   IF (coef < 0.0) THEN
     GO TO   100
   ELSE IF (coef == 0.0) THEN
     GO TO   220
   ELSE
     GO TO   110
   END IF
   100 nneg = nneg + 1
 END DO
 IF (nneg == 1 .OR. nneg == 3) GO TO 220
 IF (ip /= 1) GO TO 220
 
 DO  i = 1,4
   ecpt(i+50) = ecpt(i)
 END DO
 DO  i = 7,21
   ecpt(i+49) = ecpt(i)
 END DO
 ecpt(55) = ecpt(6)*2.0
 iret = 1
 GO TO 170
 140 ecpt(54) = ecpt( 5)
 ecpt(68) = ecpt(23)
 ecpt(69) = ecpt(24)
 ecpt(70) = ecpt(25)
 iret = 2
 GO TO 170
 150 ecpt(53) = ecpt( 4)
 ecpt(64) = ecpt(19)
 ecpt(65) = ecpt(20)
 ecpt(66) = ecpt(21)
 iret = 3
 GO TO 170
 160 ecpt(52) = ecpt( 3)
 ecpt(60) = ecpt(15)
 ecpt(61) = ecpt(16)
 ecpt(62) = ecpt(17)
 iret = 4
 
!     EACH CSLOT3 ELEMENT OR SUBELEMENT IS FORMULATED AS FOLLOWS
 
 170 IF (necpt(52) /= npvt .AND. necpt(53) /= npvt .AND.  &
     necpt(54) /= npvt) GO TO 200
 coef = 0.0
 a2   = 0.0
 DO  i = 1,3
   j    = i + 1
   IF (j > 3) j = j - 3
   k    = j + 1
   IF (k > 3) k = k - 3
   lri  = 4*i + 56
   lrj  = 4*j + 56
   lrk  = 4*k + 56
   coef = coef + ecpt(lri+2)
   fir(i) = ecpt(lrk  ) - ecpt(lrj  )
   fiz(i) = ecpt(lrj+1) - ecpt(lrk+1)
   a2   = a2 + ecpt(lri)*fiz(i)
   IF (necpt(i+51) == npvt) ipvt = i
 END DO
 IF (a2 == 0.0D0) GO TO 220
 coef = coef*ecpt(57)/(6.0D0*ecpt(55)*DABS(a2))
 i    = npvt
 DO  j = 1,3
   k    = necpt(j+51)
   kij  = coef*(fir(ipvt)*fir(j) + fiz(ipvt)*fiz(j))
   CALL sma1b( kij,k,i,ifile,0.0D0)
 END DO
 200 SELECT CASE ( iret )
   CASE (    1)
     GO TO 140
   CASE (    2)
     GO TO 150
   CASE (    3)
     GO TO 160
   CASE (    4)
     GO TO 210
 END SELECT
 210 RETURN
 
 220 WRITE  (out,230) ufm,necpt(1)
 230 FORMAT (a23,' 2160, BAD GEOMETRY OR ZERO COEFFICIENT FOR SLOT ',  &
     'ELEMENT NUMBER',i18)
 nogo =.true.
 RETURN
END SUBROUTINE kslot
