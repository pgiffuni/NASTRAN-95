SUBROUTINE sslot1(iopt)
!                  IOPT-  CSLOT3 = 0,  CSLOT4 = 1
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
 
 INTEGER, INTENT(IN)                      :: iopt
 INTEGER :: necpt(100),nout(100)
 DIMENSION sv(24)
 COMMON/sdr2x5/ ecpt(100),out(100)
 COMMON/sdr2x6/ r(4),z(4),rho,fact,a,nc1,nc2,nc3,iret,nr,  &
     col,nr1,nr2,nr3,ii,ij
 EQUIVALENCE (ecpt(1),necpt(1)),(out(1),nout(1)),(out(6),sv(1))
 
 iot = 6
 DO  i=1,30
   nout(i) = 0
 END DO
 nc1 =1
 nc2 =2
 nc3 =3
 IF(iopt /= 0) GO TO 30
!     SET UP FOR THE SLOT3 ELEMENT
!****
 rho = ecpt(5)
 iret=4
 DO  i=1,3
   nr= 4*(i-1)+10
   r(i)=ecpt(nr)
   z(i)=ecpt(nr+1)
 END DO
 GO TO 80
!****
!     THE CSLOT4 ELEMENT IS CALCULATED AS FOLLOWS
 30 CONTINUE
 rho = ecpt(6)*4.0
 DO  i =1,4
   nr = 4*(i-1) +11
   r(i) = ecpt(nr)
   z(i) = ecpt(nr+1)
 END DO
 ncol=6
 iret =1
 GO TO 80
 50 nc3 =4
 iret=2
 GO TO 80
 60 nc2 =3
 iret=3
 GO TO 80
 70 nc1= 2
 iret=4
 80 a =    (r(nc1)*(z(nc2)-z(nc3)) +r(nc2)*(z(nc3)-z(nc1))  &
     + r(nc3)*(z(nc1)-z(nc2)) )
 fact = - rho *a
 sv(nc1)=(z(nc2) -z(nc3))/fact+sv(nc1)
 sv(nc2)=(z(nc3) -z(nc1))/fact+sv(nc2)
 sv(nc3)=(z(nc1) -z(nc2))/fact+sv(nc3)
 
 nr1= 3+iopt +nc1
 nr2= 3+iopt +nc2
 nr3= 3+iopt +nc3
 
 sv(nr1)=(r(nc3)-r(nc2))/fact +sv(nr1)
 sv(nr2)=(r(nc1)-r(nc3))/fact +sv(nr2)
 sv(nr3)=(r(nc2)-r(nc1))/fact +sv(nr3)
 
 SELECT CASE ( iret )
   CASE (    1)
     GO TO 50
   CASE (    2)
     GO TO 60
   CASE (    3)
     GO TO 70
   CASE (    4)
     GO TO 90
 END SELECT
 
 90 CONTINUE
 nr = iopt+3
 IF(iopt == 1) rho =rho/4.0
 DO  i =1,nr
   j=i+1
   IF(j > iopt+3) j =j-iopt-3
   fact = 1.0/(SQRT((r(j)-r(i))**2+(z(j)-z(i))**2)*rho)
   ii=  iopt*(i+1) +4*i+3
   fact = 1.0/(SQRT((r(j)-r(i))**2+(z(j)-z(i))**2)*rho)
   ii=  iopt*(i+1) +4*i+3
   ij = ii +j -i
   sv(ii) = fact
   sv(ij) =-fact
 END DO
 
!*****
!     WRAP UP OUTPUT
!*****
 nout(1)= necpt(1)
 nout(2)= necpt(2)
 nout(3)= necpt(3)
 nout(4)= necpt(4)
 IF(iopt > 0) nout(5)= necpt(5)
 RETURN
END SUBROUTINE sslot1
