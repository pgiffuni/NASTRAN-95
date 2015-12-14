SUBROUTINE fa1ke (scr1,kfreq,bref,rho,rref,floop,nloop)
     
 
 INTEGER, INTENT(IN OUT)                  :: scr1
 REAL, INTENT(IN OUT)                     :: kfreq
 REAL, INTENT(IN OUT)                     :: bref
 REAL, INTENT(IN OUT)                     :: rho
 REAL, INTENT(IN OUT)                     :: rref
 INTEGER, INTENT(IN OUT)                  :: floop
 INTEGER, INTENT(IN)                      :: nloop
 INTEGER :: mhh,khh, sysbuf,mout,NAME(2)
 INTEGER :: buf1,trl(7),fsave
 
 REAL :: k2b2
 
 COMPLEX :: cz(1)
 
 COMMON /system/ sysbuf
 COMMON /zzzzzz/ z(1)
 COMMON /unpakx/ iout,inn,nnn,incr1
 
 EQUIVALENCE     (z(1),cz(1))
 
 DATA NAME /4HFA1K,4HE   /
 DATA khh  /101/, mhh /103/ , mout /203/ , fsave /201/
 
!     INITILIZE ON FIRST LOOP
 
 IF(floop > 1) GO TO 100
 ncore = korsz(z)
 trl(1)= khh
 CALL rdtrl(trl)
 n  = trl(3)
 nn = n*n
 n2 = n*2
 im = nn * 2
 
!       LOC     SIZE     USE
 
!     IAM1K     N*N*2    A-1 K
!     IKC       N*N*2    K   SCRATCH FOR ALLMAT
!     IMS       N*N*2    M + Q   LAMBDA FOR ALLMAT
!     IPM       N*N*2    M   HELD IN CORE
!     IPK       N*N*2    K   BETWEEN LOOPS
 
 iam1k = 1
 ikc = iam1k + im
 ims = ikc + im
 icp = ims + im
 IF(im*5+sysbuf > ncore) CALL mesage(-8,0,NAME)
 ipm = ncore - im
 ipk = ipm - im
 ncore= ipk -1
 buf1 = ncore - sysbuf
 iout = 3
 inn  = 1
 nnn  = n
 incr1= 1
 
!     PUT MHH AND KHH IN CORE
 
 ifl = khh
 ji  = ipk
 10 CALL gopen(ifl,z(buf1),0)
 DO  i=1,n
   CALL unpack(*15,ifl,z(ji))
   GO TO 16
   15 CALL zeroc(z(ji),n2)
   16 ji = ji + n2
 END DO
 CALL CLOSE(ifl,1)
 IF(ifl == mhh) GO TO 40
 ifl = mhh
 ji  = ipm
 GO TO 10
 
!     WRITE A HEADER ON MOUT
 
 40 CALL gopen(mout,z(buf1),1)
 CALL CLOSE(mout,2)
 
!              2  2
!     SOLVE   K /B  MHH + (RHO*RREF)/2.0 QHH     KHH
 
 100 k2b2 = (kfreq*kfreq) /(bref*bref)
 rr2  = (rho*rref) / 2.0
 iout = 3
 inn  = 1
 nnn  = n
 incr1= 1
 DO  i=1,ikc
   z(i) = 0.0
 END DO
 ji = ims
 CALL gopen(scr1,z(buf1),0)
 DO  i=1,n
   CALL unpack(*115,scr1,z(ji))
   115 ji = ji+n2
 END DO
 CALL CLOSE(scr1,1)
 ick = ikc -1
 ikp = ipk -1
 imp = ipm -1
 ism = ims -1
 j   = nn*2
 DO  i=1,j
   z(i+ism) = z(i+ism) * rr2 + z(i+imp) * k2b2
   z(i+ick) = - z(i+ikp)
 END DO
 CALL incore(z(ims),n,z(ikc),z(iam1k),n)
 
!     GET EIGENVALUES FROM ALLMAT
 
 im = ims + n2
 in = im + n2
 l  = 0
 CALL allmat(z(iam1k),z(ims),z(ikc),0,0,z(im),0,z(in),n,l,0)
 
!     WRITE OUT EIGENVALUES ON MOUT
 
 im = ims/2
 nl = 2*l
 DO  i = 1,l
   IF(cz(i+im) /= (0.0,0.0))cz(i+im) = CSQRT(cz(i+im))
   IF(AIMAG(cz(i+im)) < 0.0) cz(i+im) = - cz(i+im)
 END DO
 CALL gopen(mout,z(buf1),3)
 CALL WRITE(mout,z(ims),nl,1)
 IF(floop >= nloop) GO TO 200
 CALL CLOSE(mout,3)
 RETURN
 
!     LAST LOOP BUILD FSAVE
 
 200 CALL CLOSE(mout,1)
 ibuf2 = buf1 - sysbuf
 CALL gopen(mout,z(buf1),0)
 CALL gopen(fsave,z(ibuf2),0)
 CALL skprec(fsave,3)
 CALL CLOSE(fsave,2)
 CALL gopen(fsave,z(ibuf2),3)
 210 CALL READ(*230,*220,mout,z(1),ibuf2,1,nwr)
 220 CALL WRITE(fsave,z(1),nwr,1)
 GO TO 210
 230 CALL CLOSE(mout,1)
 CALL CLOSE(fsave,1)
 trl(1) = fsave
 trl(2) = nloop
 trl(7) = l
 CALL wrttrl(trl)
 RETURN
END SUBROUTINE fa1ke
