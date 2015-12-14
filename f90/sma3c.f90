SUBROUTINE sma3c(iflag,k)
     
!     THIS ROUTINE WILL MERGE ZINVS,ZS,STZ,AND STZS INTO KE AND
!       BUILD KE UP TO G SIZE.  IF INFLAG .LT. 0 THERE ARE NO
!       UD-S
 
 
 INTEGER, INTENT(IN OUT)                  :: iflag
 INTEGER, INTENT(OUT)                     :: k(7)
 DOUBLE PRECISION :: a11,b11,d11
 INTEGER :: zinvs,zs,stz,stzs,gei,sysbuf,NAME(2),iz(1)
 DIMENSION block1(20),block2(20)
 
 COMMON /BLANK/luset
 COMMON /zzzzzz/ z(1)
 COMMON /zblpkx/ d11(2),id
 COMMON /system/ sysbuf, dummy(53), iprec
 COMMON /genely/gei,dum(4),stzs(1),zinvs(1),zs(1),stz(1),dum1(62), m,n
 
 EQUIVALENCE (z(1),iz(1))
 DATA NAME / 4HSMA3,4HC    /
 
!     IUI IS POINTER TO UI SET, IUD IS POINTER TO UD SET
 
 iui =1
 iud =m+1
 nz = korsz(z)
 
!     OPEN GEI(WITHOUT REWIND)
 
 nz = nz -sysbuf
 CALL gopen(gei,z(nz+1),2)
 
!     READ IN UI SET
 
 CALL fread(gei,z,-3,0)
 CALL fread(gei,z, m,0)
 
!     READ IN UD
 
 IF (iflag < 0) GO TO 10
 CALL fread(gei,z(iud),n,1)
 
!     OPEN BUFFERS FOR MATRICES
 
 10 llen = m+n+2*sysbuf
 IF(iflag >= 0) llen = llen+3*sysbuf
 IF (llen > nz) GO TO 220
 nz = nz-sysbuf
 CALL gopen(k,z(nz+1),1)
 nz = nz -sysbuf
 CALL gopen(zinvs,z(nz+1),0)
 IF (iflag < 0) GO TO 20
 nz =nz -sysbuf
 CALL gopen(zs,z(nz+1),0)
 nz =nz -sysbuf
 CALL gopen(stz,z(nz+1),0)
 nz =nz -sysbuf
 CALL gopen(stzs,z(nz+1),0)
 
!     LOOP ON LUSET MAKING COLUMNS OF KGG
 
 20 k(2) = 0
 k(3) = luset
 k(4) = 6
 k(5) = 2
 k(6) = 0
 k(7) = 0
 iip = 0
 idp = 0
 DO  i=1,luset
   CALL bldpk (2, iprec, k(1), 0, 0)
   IF( iip >= m )  GO TO 25
   l = iui + iip
   IF (i == iz(l)) GO TO 30
   25 CONTINUE
   IF (iflag < 0) GO TO 160
   IF( idp >= n )  GO TO 160
   l = iud + idp
   IF (i == iz(l)) GO TO 40
   GO TO 160
   
!     USING UI -- ZINVS AND STZ
   
   30 iip = iip +1
   nam1 = zinvs(1)
   nam2 = stz(1)
   GO TO 50
   
!     USING UD ZS AND STZS
   
   40 idp = idp +1
   nam1 = zs(1)
   nam2 = stzs(1)
   
!     MERGE ROUTINE FOR COLUMN
   
   50 iad = 0
   ibd = 0
   ihop = 0
   CALL intpk(*140,nam1,block1(1),2,1)
   60 IF(iflag < 0) GO TO 150
   CALL intpk(*150,nam2,block2(1),2,1)
   70 CALL intpki(a11,ia,nam1,block1(1),iaeol)
   l=  iui +ia -1
   ii = iz(l)
   IF (ihop == 1) GO TO 90
   ihop = 1
   80 CALL intpki(b11,ib,nam2,block2(1),ibeol)
   l = iud +ib -1
   jj = iz(l)
   90 IF (ii-jj < 0) THEN
     GO TO   100
   ELSE IF (ii-jj == 0) THEN
     GO TO   320
   ELSE
     GO TO   120
   END IF
   
!     PUT IN A11
   
   100 d11(1) =a11
   id = ii
   CALL zblpki
   IF (iaeol == 0) THEN
     GO TO    70
   END IF
   110 iad = 1
   ii = 99999
   IF(ibd == 0) THEN
     GO TO   120
   ELSE
     GO TO   160
   END IF
   
!     PUT IN BUU
   
   120 d11(1) = b11
   id = jj
   CALL zblpki
   IF (ibeol == 0) THEN
     GO TO    80
   END IF
   130 ibd = 1
   jj = 99999
   IF(iad == 0) THEN
     GO TO   100
   ELSE
     GO TO   160
   END IF
   
!     NULL NAM1
   
   140 iad =1
   ii = 99999
   GO TO 60
   
!     NO NAM2
   
   150 ibd =1
   jj = 99999
   ihop =1
   GO TO 70
   
!     END OF COLUMN
   
   160 CALL bldpkn(k(1),0,k)
   
!     END LOOP
   
 END DO
 CALL wrttrl (k)
 CALL CLOSE (k(1),1)
 CALL CLOSE (zinvs(1),1)
 IF (iflag < 0) GO TO 180
 CALL CLOSE (stz(1),1)
 CALL CLOSE (stzs(1),1)
 CALL CLOSE (zs(1),1)
 180 RETURN
 
!     ERROR MESAGES
 
 220 CALL mesage(-8,gei,NAME)
 320 CALL mesage(-7,0,NAME)
 RETURN
END SUBROUTINE sma3c
