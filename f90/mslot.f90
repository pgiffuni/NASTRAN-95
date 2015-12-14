SUBROUTINE mslot(itype)
!*****
!     THIS ROUTINE CALCULATES THE MASS MATRIX TERMS FOR THE
!         CSLOT3 AND CSLOT4 TWO DIMENSIONAL LAPLACE ELEMENTS
!                  IOPT-  CSLOT3 = 0,  CSLOT4 = 1
!*****
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
!***** 25                           W4
!***** 26                           TEMP
 
 
 INTEGER, INTENT(IN OUT)                  :: itype
 INTEGER :: necpt(100)
 DOUBLE PRECISION :: coef          ,a2            ,wb  &
     ,r             ,z             ,w ,mij
!*****
 COMMON  /sma2cl/  iopt4,k4ggsw,npvt
 COMMON  /sma2et/  ecpt(100)
 COMMON  /sma2io/  dum1(10),ifile
 COMMON/sma2dp/     coef          ,a2            ,wb  &
     ,r(3)          ,z(3)          ,w(3) ,mij           ,iret          ,ip  &
     ,k
!*****
 EQUIVALENCE  (ecpt(1),necpt(1))
!*****
 IF(itype > 0) GO TO 50
 IF(ecpt(6) == 0.0.OR.necpt(7) == 0 ) RETURN
 k=-1
 10 k=k+1
 IF(2*necpt(8) - k*necpt(7)   < 0) THEN
   GO TO    30
 ELSE IF (2*necpt(8) - k*necpt(7)   == 0) THEN
   GO TO    20
 ELSE
   GO TO    10
 END IF
 20 necpt(7) = necpt(7)*2
 30 ecpt(7) = FLOAT(necpt(7))/2.0
 DO  i=1,20
   ecpt(i+50)= ecpt(i)
 END DO
 iret =4
 GO TO 140
!*****
!     THE CSLOT4 ELEMENT IS CHECKED FOR VALIDITY AND THE DATA ARE
!     REARRANGED TO CONFORM TO THE CSLOT3 FORMAT
!*****
 50 IF(ecpt(7) == 0.0 .OR.necpt(8) == 0 ) RETURN
 k =-1
 60 k =k+1
 IF( 2*necpt(9) - k*necpt(8)  < 0) THEN
   GO TO    80
 ELSE IF ( 2*necpt(9) - k*necpt(8)  == 0) THEN
   GO TO    70
 ELSE
   GO TO    60
 END IF
 70 necpt(8) =necpt(8)*2
 80 ecpt(8) = FLOAT(necpt(8))/2.0
 DO  i=1,4
   ecpt(i+50) = ecpt(i)
 END DO
 DO  i=6,21
   ecpt(i+49) = ecpt(i)
 END DO
 ecpt(56) = ecpt(7)*2.0
 iret =1
 GO TO 140
 110 ecpt(54)= ecpt(5)
 ecpt(68) = ecpt(23)
 ecpt(69) = ecpt(24)
 ecpt(70) = ecpt(25)
 iret =2
 GO TO 140
 120 ecpt(53)= ecpt(4)
 ecpt(64)= ecpt(19)
 ecpt(65)= ecpt(20)
 ecpt(66)= ecpt(21)
 iret =3
 GO TO 140
 130 ecpt(52)= ecpt(3)
 ecpt(60)= ecpt(15)
 ecpt(61)= ecpt(16)
 ecpt(62)= ecpt(17)
 iret =4
!*****
!     EACH CSLOT3 ELEMENT OR SUBELEMENT IS FORMULATED AS FOLLOWS
!*****
 140 IF((necpt(52) /= npvt).AND.(necpt(53) /= npvt).AND.  &
     (necpt(54) /= npvt)) GO TO 170
 DO  i=1,3
   ip = 4*(i-1)+60
   r(i) =ecpt(ip)
   z(i) =ecpt(ip+1)
   w(i) =ecpt(ip+2)
   IF(npvt == necpt(i+51)) ipvt=i
 END DO
 a2 = (r(2)-r(1))*(z(3)-z(1))  -(r(3)-r(1))*(z(2)-z(1))
 wb = w(1) +w(2) +w(3)+w(ipvt)
 coef = DABS(a2)*ecpt(57) /(120.0D0 *ecpt(56))
 i=npvt
 DO  j=1,3
   k = necpt(j+51)
   mij = coef *( wb + w(j) )
   IF (ipvt == j) mij =mij*2.0D0
   CALL sma2b(mij,k,i,ifile,0.0D0)
 END DO
 170 SELECT CASE ( iret )
   CASE (    1)
     GO TO 110
   CASE (    2)
     GO TO 120
   CASE (    3)
     GO TO 130
   CASE (    4)
     GO TO 180
 END SELECT
 180 RETURN
END SUBROUTINE mslot
