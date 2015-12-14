SUBROUTINE gpwg1b (mo,ogpwg,wtmass,ipoint)
     
!     DOUBLE PRECISION VERSION, BY G.CHAN/UNISYS  8/86
 
!     THIS ROUTINE WRITES OGPWG--
!         HEADER
!         MO =  36 D.P.WORDS
!         S  =  9  D.P.WORDS
!         MX,XX,YX,ZX,MY,XY,YY,ZY,MZ,XZ,YZ,ZZ  = 12 D.P.WORDS
!         I  =  9  D.P.WORDS
!         I1P, I2P, I3P = 3 D.P.WORDS
!         Q  =  9  D.P.WORDS
!               78 D.P.WORDS (156 S.P.WORDS) TOTAL
 
 
 INTEGER, INTENT(IN OUT)                  :: mo
 INTEGER, INTENT(IN OUT)                  :: ogpwg
 REAL, INTENT(IN)                         :: wtmass
 INTEGER, INTENT(IN)                      :: ipoint
 DOUBLE PRECISION :: s(3,3),mt(3,3),mtr(3,3),mr(3,3),temp(3,3),  &
                     dz(36),delta,epsi
 INTEGER :: sysbuf, NAME(2),z(150)
 EQUIVALENCE       (dz(1),z(1),iz(1))
 
 COMMON /zzzzzz/ iz(1)
 COMMON /unpakx/ it1,ii,jj,incr
 COMMON /system/ sysbuf
 COMMON /output/ head(1)
 
 DATA    NAME  / 4HGPWG,4H1B   /
 
!     ASSIGN BUFFER
!     OPEN OGPWG, PUT ON OFP HEADER
 
 ibuf = korsz(z)- sysbuf+1
 CALL gopen (mo,z(ibuf),0)
 
!     UNPACK MO  + MOVE TO PARTITIONS
 
 it1  = 2
 incr = 1
 jj   = 6
 ii   = 1
 k    = 1
 DO  i=1,6
   CALL unpack (*10,mo,dz(k))
   GO TO 30
   10 DO  l=1,6
     m = l+k-1
     dz(m) =0.0D0
   END DO
   30 k = k+6
 END DO
 CALL CLOSE (mo,1)
 delta=1.d0/wtmass
 DO  i=1,36
   dz(i) = dz(i)*delta
 END DO
 
!     OPEN OGPWG FOR OUTPUT
 
 CALL gopen (ogpwg,z(ibuf),1)
 DO  i = 104,150
   z(i) = 0
 END DO
 z(101) = 1
 z(102) = 13
 z(103) = ipoint
 z(110) = 78*2
 CALL WRITE (ogpwg,z(101),50,0)
 CALL WRITE (ogpwg,head,96,1)
 
!     PUT MO  ON OGPWG
 
 CALL WRITE (ogpwg,z(1),72,0)
 
!     PARTITION MO INTO MT, MTR, AND MR
!     AND CREATE DIAGONAL S MATRIX
 
 mt(1,1) = dz(1)
 mt(1,2) = dz(2)
 mt(1,3) = dz(3)
 mt(2,1) = dz(7)
 mt(2,2) = dz(8)
 mt(2,3) = dz(9)
 mt(3,1) = dz(13)
 mt(3,2) = dz(14)
 mt(3,3) = dz(15)
 mtr(1,1)= dz(4)
 mtr(2,1)= dz(5)
 mtr(3,1)= dz(6)
 mtr(1,2)= dz(10)
 mtr(2,2)= dz(11)
 mtr(3,2)= dz(12)
 mtr(1,3)= dz(16)
 mtr(2,3)= dz(17)
 mtr(3,3)= dz(18)
 mr(1,1) = dz(22)
 mr(1,2) = dz(23)
 mr(1,3) = dz(24)
 mr(2,1) = dz(28)
 mr(2,2) = dz(29)
 mr(2,3) = dz(30)
 mr(3,1) = dz(34)
 mr(3,2) = dz(35)
 mr(3,3) = dz(36)
 s(1,1)  = 1.0D0
 s(1,2)  = 0.0D0
 s(1,3)  = 0.0D0
 s(2,1)  = 0.0D0
 s(2,2)  = 1.0D0
 s(2,3)  = 0.0D0
 s(3,1)  = 0.0D0
 s(3,2)  = 0.0D0
 s(3,3)  = 1.0D0
 
!     COMPUTE  DETERMINATE OF  MT
 
 delta = DSQRT(mt(1,1)**2 + mt(2,2)**2 + mt(3,3)**2)
 epsi  = DSQRT(mt(2,1)**2 + mt(3,1)**2 + mt(3,2)**2)
 IF (epsi  == 0.0D0) GO TO 60
 epsi = epsi/delta
 IF (delta == 0.0D0) GO TO 45
 IF (epsi < 1.0D-6) GO TO 60
 
!     ROTATE COORDINATES
 
 45 r = epsi
 CALL mesage (42,r,NAME)
 DO  i=1,3
   DO  j=1,3
     temp(i,j)= mt(i,j)
   END DO
 END DO
 
!     COMPUTE EIGENVECTORS OF  MT  BY JACOBY  METHOD
 
 CALL gpwg1c (temp,s,dz(1),iflag)
 IF (iflag > 0) CALL mesage(-7,0,NAME)
 
!     ORDER EIGENVECTORS  SUCH THAT
 
!     TRANSFORM  MT
 
 CALL gmmatd (mt,3,3,0,s,3,3,0,temp)
 CALL gmmatd (s,3,3,1,temp,3,3,0,mt)
 
!     TRANSFORM  MTR
 
 CALL gmmatd (mtr,3,3,0,s,3,3,0,temp)
 CALL gmmatd (s,3,3,1,temp,3,3,0,mtr)
 
!     TRANSFORM  MR
 
 CALL gmmatd (mr,3,3,0,s,3,3,0,temp)
 CALL gmmatd (s,3,3,1,temp,3,3,0,mr)
 
!     OUTPUT S
 
 60 CALL WRITE (ogpwg,s,18,0)
 
!     COMPUTE   MX,XX,YX,ZX
 
 dz(1) = mt(1,1)
 dz(2) = 0.0D0
 dz(3) = 0.0D0
 dz(4) = 0.0D0
 IF (dz(1) == 0.0D0) GO TO 70
 dz(2) = mtr(1,1)/dz(1)
 dz(3) =-mtr(3,1)/dz(1)
 dz(4) = mtr(2,1)/dz(1)
 70 CALL WRITE (ogpwg,dz(1),8,0)
 dz(5) = mt(2,2)
 dz(6) = 0.0D0
 dz(7) = 0.0D0
 dz(8) = 0.0D0
 IF (dz(5) == 0.0D0) GO TO 80
 dz(6) = mtr(3,2)/dz(5)
 dz(7) = mtr(2,2)/dz(5)
 dz(8) =-mtr(1,2)/dz(5)
 80 CALL WRITE (ogpwg,dz(5),8,0)
 dz( 9) = mt(3,3)
 dz(10) = 0.0D0
 dz(11) = 0.0D0
 dz(12) = 0.0D0
 IF (dz(9) == 0.0D0) GO TO 90
 dz(10) =-mtr(2,3)/dz(9)
 dz(11) = mtr(1,3)/dz(9)
 dz(12) = mtr(3,3)/dz(9)
 90 CALL WRITE (ogpwg,dz(9),8,0)
 
!     COMPUTE INERTIAS
 
 temp(1,1) = mr(1,1) - dz(5)*dz(8)*dz(8) - dz(9)*dz(11)*dz(11)
 temp(2,1) =-mr(1,2) - dz(9)*dz(10)*dz(11)
 temp(1,2) = temp(2,1)
 temp(1,3) =-mr(1,3) - dz(5)*dz(6)*dz(8)
 temp(3,1) = temp(1,3)
 temp(2,2) = mr(2,2) - dz(9)*dz(10)*dz(10) - dz(1)*dz(4)*dz(4)
 temp(2,3) =-mr(2,3) - dz(1)*dz(3)*dz(4)
 temp(3,2) = temp(2,3)
 temp(3,3) = mr(3,3) - dz(1)*dz(3)*dz(3) - dz(5)*dz(6)*dz(6)
 
 CALL WRITE (ogpwg,temp,18,0)
 CALL gpwg1c (temp,s,dz(1),iflag)
 IF (iflag > 0) CALL mesage(-7,0,NAME)
 
!     PUT OUT  PRINCIPLE INERTIA-S
 
 CALL WRITE (ogpwg,dz(1),6,0)
 
!     PUT  OUT  Q
 
 CALL WRITE (ogpwg,s,18,0)
 CALL clstab (ogpwg,1)
 
 RETURN
END SUBROUTINE gpwg1b
