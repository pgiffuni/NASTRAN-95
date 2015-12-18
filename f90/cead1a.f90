SUBROUTINE cead1a (lami,phidi,phidli,lamd,phid,phidl,nfound,nvect,  &
        capp)
     
!     ROUTINE SORTS LAMI, PHIDI AND PHIDLI (INV. POWER), BASED ON LAMI,
!     AND CREATES LAMD, PHID AND PHIDL
 
 
 INTEGER, INTENT(IN)                      :: lami
 INTEGER, INTENT(IN)                      :: phidi
 INTEGER, INTENT(IN)                      :: phidli
 INTEGER, INTENT(IN)                      :: lamd
 INTEGER, INTENT(IN)                      :: phid
 INTEGER, INTENT(IN)                      :: phidl
 INTEGER, INTENT(IN)                      :: nfound
 INTEGER, INTENT(IN)                      :: nvect
 INTEGER, INTENT(IN)                      :: capp
 DOUBLE PRECISION :: zd(1),d1,d2
 INTEGER :: sysbuf,iz(1),ih(7),FILE
 INTEGER :: NAME(2)
 INTEGER :: det,hes
 INTEGER :: filek,filem,fileb
 
 
 COMMON /system/ ksystm(65)
 COMMON /output/head(1)
 COMMON /zzzzzz/ z(1)
 COMMON /cinvpx/ filek(7),filem(7),fileb(7)
 COMMON /condas/    pi       ,twopi    ,radeg    ,degra    , s4pisq
 COMMON /packx/ it1,it2,ii,jj,incur
 
 EQUIVALENCE ( ksystm( 1) , sysbuf )
 EQUIVALENCE (iz(1),z(1)),(z(1),zd(1))
 
 DATA NAME/4HCEAD,4H1A  /
 DATA ih / 7*0 /
 DATA det,hes / 4HDET ,4HHESS /
 
!     INITIALIZE POINTER ARRAY
 
 DO  i=1,nfound
   iz(i)= i
 END DO
 
!     BRING IN  EIGENVALUES
 
 ilama =(nfound+1)/2  +1
 ibuf = korsz(iz)-sysbuf+1
 FILE = lami
 CALL OPEN(*170,lami,iz(ibuf),0)
 k =ilama
 DO  i=1,nfound
   CALL READ(*190,*200,lami,zd(k),4,1,iflag)
   k= k +2
 END DO
 CALL CLOSE(lami,1)
 IF(nfound == 1) GO TO 70
 
 
!     SORT ON SIGN IMAGINARY THEN ON MAG IMAG
 
 jj = nfound-1
 DO  i=1,jj
   ii = i+1
   m  = ilama+ 2*i -2
   DO  j=ii,nfound
     l =  ilama +2*j-2
     
!     SIGN IMAG
     
     d1 = DSIGN(1.0D0,zd(l+1))
     d2 = DSIGN(1.0D0,zd(m+1))
     IF(d1 == d2) GO TO 40
     IF( d1 == 1.0D0)CYCLE
     
!     SWITCH
     
     30 d1 = zd(l)
     zd(l) =zd(m)
     zd(m)=d1
     d1 =zd(l+1)
     zd(l+1)=zd(m+1)
     zd(m+1)=d1
     it1 = iz(j)
     iz(j)= iz(i)
     iz(i)= it1
     CYCLE
     
!     TEST MAGNITIDE IMAG
     
     40 IF(DABS(zd(l+1)) -DABS(zd(m+1)) < 0.0) THEN
       GO TO    30
     ELSE
       GO TO    50
     END IF
     50 CONTINUE
   END DO
 END DO
 
!     PUT OUT LAMA-S IN ORDER GIVEN BY LIST
 
 70 CALL gopen(lamd,iz(ibuf),1)
 ih(2) =1006
 ih(1) = 90
 CALL WRITE(lamd,  ih, 4,0)
 ih(6) = 6
 CALL WRITE(lamd,  ih, 6,0)
 CALL WRITE(lamd,  iz,40,0)
 CALL WRITE(lamd,head,96,1)
 l = 5*nfound +2
 DO  i=1,nfound
   iz(l)=i
   iz(l+1)=iz(i)
   k = 2*i-2+ilama
   z(l+2) =zd(k)
   z(l+3) = zd(k+1)
   z(l+4) = 0.0
   z(l+5) = 0.0
   IF(ABS(z(l+3)) <= 1.0E-3*ABS(z(l+2))) GO TO 80
   z(l+4) = ABS(z(l+3))/twopi
   z(l+5) = -2.0*z(l+2)/ABS(z(l+3))
   80 CALL WRITE(lamd,iz(l),6,0)
 END DO
 CALL CLOSE(lamd,1)
 ih(1) =lamd
 CALL wrttrl(ih)
 
!     BRING IN PHIDI IN ORDER NEEDED AND OUTPUT
 
 ibuf1 = ibuf -sysbuf
 CALL gopen(phid,iz(ibuf1),1)
 it1 = 4
 it2 = 3
 incur =1
 ii =1
 ih(1)=phid
 ih(2)= 0
 ih(4) =2
 ih(5) =3
 ih(6) = 0
 k = 1
 101 IF(iz(k) <= nvect) GO TO 111
 k = k+1
 GO TO 101
 111 FILE = phidi
 ipos =1
 CALL OPEN(*170,phidi,iz(ibuf),0)
 DO  i=1,nvect
   IF (nvect == 1) GO TO 130
   100 l= iz(i)-ipos
   IF(l < 0) THEN
     GO TO   150
   ELSE IF (l == 0) THEN
     GO TO   130
   END IF
   110 CALL skprec(phidi,l)
   
!     BRING IN EIGENVECTORS
   
   130 CALL READ(*190,*140,phidi,zd(ilama),ibuf1-1,0,m)
   GO TO 210
   140 jj= m/4
   ipos = iz(k) +1
   CALL pack(zd(ilama),phid,ih)
   GO TO 159
   
!     PAST VECTOR NEEDED
   
   150 CALL REWIND(phidi)
   ipos =1
   GO TO 100
   159 k = k+1
 END DO
 CALL CLOSE(phid,1)
 CALL CLOSE(phidi,1)
 ih(3) =jj
 CALL wrttrl(ih)
 
!     OUTPUT PHIDL IF NOT PURGED AND IF AT LEAST ONE INPUT MATRIX IS
!     UNSYMMETRIC
 
 ih(1) = phidl
 CALL rdtrl (ih)
 IF (ih(1) < 0) RETURN
 IF (capp /= det .AND. capp /= hes) GO TO 301
 filek(1) = 101
 CALL rdtrl (filek)
 filem(1) = 103
 CALL rdtrl (filem)
 fileb(1) = 102
 CALL rdtrl (fileb)
 301 IF (filek(1) > 0 .AND. filek(4) /= 6) GO TO 302
 IF (filem(1) > 0 .AND. filem(4) /= 6) GO TO 302
 IF (fileb(1) > 0 .AND. fileb(4) /= 6) GO TO 302
 RETURN
 302 CALL gopen (phidl,iz(ibuf1),1)
 CALL makmcb (ih,phidl,0,2,3)
 IF (capp /= det .AND. capp /= hes) GO TO 305
 CALL clvec (lamd,nvect,phidl,ih,ibuf,ibuf1)
 GO TO 395
 305 k = 1
 310 IF (iz(k) <= nvect) GO TO 320
 k = k + 1
 GO TO 310
 320 FILE = phidli
 ipos = 1
 CALL OPEN(*170,phidli,iz(ibuf),0)
 DO  i=1,nvect
   IF (nvect == 1) GO TO 350
   330 l = iz(i) - ipos
   IF (l < 0) THEN
     GO TO   370
   ELSE IF (l == 0) THEN
     GO TO   350
   END IF
   340 CALL skprec (phidli,l)
   
!     BRING IN LEFT EIGENVECTORS
   
   350 CALL READ(*190,*360,phidli,zd(ilama),ibuf1-1,0,m)
   GO TO 210
   360 jj = m/4
   ipos = iz(k) + 1
   CALL pack (zd(ilama),phidl,ih)
   GO TO 380
   
!     PAST VECTOR NEEDED
   
   370 CALL REWIND (phidli)
   ipos = 1
   GO TO 330
   380 k = k + 1
 END DO
 CALL CLOSE (phidli,1)
 395 CALL CLOSE (phidl,1)
 ih(3) = jj
 CALL wrttrl (ih)
 RETURN
 
!     ERROR MESAGES
 
 170 ip1 =-1
 180 CALL mesage(ip1,FILE,NAME)
 190 ip1 =-2
 GO TO 180
 200 ip1 = -3
 GO TO 180
 210 ip1 = -8
 GO TO 180
END SUBROUTINE cead1a
