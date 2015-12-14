SUBROUTINE fa1pki (fsave,qhhl)
     
!     FA1PKI BUILDS AN INTERPOLATION MATRIX IN CORE FOR PK METHOD
 
!     LAST REVISED  2/91, BY J.PETKAS/LOOKHEED
!     TO ALLOW CALCULATION OF INTERPOLATION MATRIX IN D.P.
 
 
 INTEGER, INTENT(IN)                      :: fsave
 INTEGER, INTENT(IN)                      :: qhhl
 INTEGER :: sysbuf,NAME(2),trl(7),buf1,floop
 REAL :: newm
 DOUBLE PRECISION :: dx1,dx2,det,dz(1)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /BLANK / floop
 COMMON /system/ sysbuf,nout
 COMMON /unpakx/ iout,inn,nnn,incr1
 COMMON /zzzzzz/ z(1)
 COMMON /fa1pkc/ ncore,nk,imvr,ik,ia,iq,icp,iflag
 EQUIVALENCE     (dz(1),z(1))
 DATA    oldm  / -1.0/
 DATA    NAME  / 4HFA1P,4HKI  /
 
 iflag = 0
 IF (oldm /= -1.0) GO TO 20
 ncore= korsz(z)
 buf1 = ncore - sysbuf
 imvr = 1
 
!     PUT M V   IN CORE ON SECOND LOOP RETURN IF SAME MACH
 
 ifle = fsave
 CALL gopen (fsave,z(buf1),0)
 CALL READ (*180,*10,fsave,z(imvr),buf1,1,nwr)
 10 ik = imvr + nwr
 CALL CLOSE (fsave,1)
 20 i = (floop-1)*3 + imvr
 newm = z(i)
 IF (oldm == newm) GO TO 200
 oldm = newm
 iflag= 1
 
!     PUT LIST OF M K'S IN CORE FOR THIS MACH
 
 trl(1) = qhhl
 CALL rdtrl (trl)
 nrow = trl(3)
 ni   = (trl(2)/trl(3))*2
 iout = 3
 inn  = 1
 incr1= 1
 nnn  = nrow
 n2   = nrow*2
 nn   = nrow*nrow
 ifle = qhhl
 CALL OPEN (*180,qhhl,z(buf1),0)
 CALL READ (*180,*180,qhhl,z,-3,0,nwr)
 CALL READ (*180,*180,qhhl,n, 1,0,nwr)
 n  = n + n
 ni = MIN0(ni,n)
 CALL READ (*180,*180,qhhl,z(ik),ni,1,nwr)
 
!     FIND M'S CLOSEST TO NEWM
 
 ia  = ik + ni
 IF (MOD(ia,2) == 0) ia = ia + 1
 rmi = 1.e20
 rms = 0.0
 DO  i = 1,ni,2
   rmx = ABS(z(ik+i-1)-newm)
   rmi = AMIN1(rmi,rmx)
   IF (rmx > rmi) CYCLE
   rms = z(ik+i-1)
 END DO
 rmi = rms
 
!     COUNT K"S
 
 nk = 0
 DO  i = 1,ni,2
   IF (z(ik+i-1) == rmi) GO TO 40
   CYCLE
   40 nk = nk + 1
 END DO
 
!     ALLOCATE CORE FOR A-1 AND Q.  THEN BUILD THEM.
 
 i  = 2*(nk+1)**2
 iq = ia + i
 icp= iq + nn*2*nk
 IF (MOD(icp,2) == 0) icp = icp + 1
 IF (icp+sysbuf+n2 > ncore) CALL mesage (-8,0,NAME)
 
!     BUILD A
 
 j = 0
 DO  i = 1,ni,2
   IF (z(ik+i-1) == rmi) GO TO 60
   CYCLE
   60 z(iq+j) = z(ik+i)
   j = j + 1
 END DO
 nk1 = nk + 1
 n   = 0
 m   = iq - 1
 iad = ia/2  + 1
 icpd= icp/2 + 1
 DO  i = 1,nk1
   dx2 = z(m+i)
   DO  j = 1,nk1
     IF (i == nk1 .AND. j == nk1) EXIT
     IF (j == nk1 .OR.  i == nk1) GO TO 75
     dx1 = z(m+j)
     dz(iad+n) = DABS((dx1-dx2)**3) + (dx1+dx2)**3
     GO TO 80
     75 dz(iad+n) = 1.0D+0
     80 n = n + 1
   END DO
 END DO
 100 dz(iad+n) = 0.0D+0
 
!     MODIFICATION FOR LEVEL 17.7 UPDATE
!     REPLACE ALL CALLS TO INVAER WITH CALLS TO INVERS.
!     NO NEED TO COMPUTE DETERMINANT SINCE IT IS NOT USED SUBSEQUENTLY.
 
 ising = -1
 CALL inverd (nk1,dz(iad),nk1,0,0,det,ising,dz(icpd))
 IF (ising == 2) GO TO 150
 
!     BUILD Q
 
 n  = 0
 in = nn
 l  = 0
 DO  i = 1,ni,2
   IF (z(ik+i-1) == rmi) GO TO 110
   CALL skprec (qhhl,nrow)
   CYCLE
   110 DO  j = 1,nrow
     CALL unpack (*115,qhhl,z(icp))
     GO TO 120
     115 CALL zeroc (z(icp),nrow*2)
     
!     SPLIT REAL AND IMAGINARY DIVIDE IMAGINARY BY K
     
     120 DO  k = 1,n2,2
       z(iq+n) = z(icp+k-1)
       n = n + 1
       z(iq+in) = z(icp+k)/z(ik+i)
       in = in + 1
     END DO
   END DO
   z(ik+l) = z(ik+i)
   l = l  + 1
   n = n  + nn
   in= in + nn
 END DO
 CALL CLOSE (qhhl,1)
 GO TO 200
 
 150 WRITE  (nout,160) ufm,NAME
 160 FORMAT (a23,' 2427, SINGULAR MATRIX FOR INTERPOLATION IN ',2A4)
 CALL mesage (-61,0,0)
 180 CALL mesage (-2,ifle,NAME)
 200 RETURN
END SUBROUTINE fa1pki
