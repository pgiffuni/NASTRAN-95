SUBROUTINE modac1(casecc,tol,tol1,casezz,caseyy)
     
!     MODAC1 REDUCES THE NUMBER OF ENTRIES ON TOL TO THE TIMES
!         SPECIFIED BY THE OFREQ SET IN CASECC
 
!     CORE IS        OUT AS FOLLOWS ON RETURN
 
!         CONTENTS            LENGTH  TYPE   POINTER
!         --------            ------  ----   -------
!         NEW TIMES           NFN      R     IFN
!         KEEP REMOVE         NFO      I     IKR
 
 INTEGER, INTENT(IN)                      :: casecc
 INTEGER, INTENT(IN)                      :: tol
 INTEGER, INTENT(IN)                      :: tol1
 INTEGER, INTENT(IN)                      :: casezz
 INTEGER, INTENT(IN OUT)                  :: caseyy
 INTEGER :: sysbuf, NAME(2), FILE,ihd(2),mcb(7),ibuf(6)
 INTEGER :: flag
 REAL :: z(1)
 COMMON /system/ sysbuf
 COMMON /modac3/ nfo,nfn,nz, id
 COMMON /zzzzzz/ iz(1)
 EQUIVALENCE   ( z(1), iz(1) )
 DATA NAME /4HMODA,4HC1  /
 
!     BRING  IN  CASECC
 
 lw = 6
 IF(id == 4) lw = 7
 ibuf1 = nz -sysbuf+1
 ibuf2 = ibuf1-sysbuf
 ibuf3 = ibuf2 -sysbuf
 CALL gopen(casecc,iz(ibuf1),0)
 FILE = casecc
 CALL READ(*900,*10,casecc,iz,ibuf2-1,0,ivec)
 CALL mesage(-8,0,NAME)
 10 CONTINUE
 icc = 0
 CALL CLOSE(casecc, 1)
 ifrout =145
 ilsym = 200
 ivec  = ivec+1
 ilist = ivec
 IF(id == 5) GO TO 600
 
!     BRING IN OLD TIME/FREQ  LIST
 
 FILE = tol
 CALL OPEN(*900,tol,iz(ibuf1),0)
 i = ilist
 m = 3
 ix= 2
 nfo = nfo + i
 IF(id == 2 .OR. id == 4) GO TO 25
 20 CALL READ(*910,*30,tol,ibuf,m,0,flag)
 iz(i) =ibuf(m)
 iz(i+1)= 0
 i =  i + ix
 m =1
 GO TO 20
 25 CALL fwdrec(*910,tol)
 CALL fwdrec(*910,tol)
 26 CALL READ(*910,*30,tol,ibuf,lw,0,flag)
 iz(i) = ibuf(4)
!     REIG SHOULD BE ON CYCLES
 IF(id == 4) iz(i) = ibuf(5)
 iz(i+1) = 0
 i = i+2
 IF(i == nfo) GO TO 30
 GO TO 26
 30 CALL CLOSE(tol,1)
 nlist = i -ix
 
!     MATCH LIST OF  SELECTED VALUES WITH TIME LIST IN CORE
 
 35 CONTINUE
 ix = icc + ifrout
 ifset = iz(ix)
 IF ( ifset  <= 0)  GO TO  70
 ix = icc +ilsym
 isetnf = ix + iz(ix)+1
 40 isetf  = isetnf +2
 nsetf  =iz(isetnf+1) + isetf-1
 IF( iz(isetnf) == ifset) GO TO 80
 isetnf = nsetf +1
 IF ( isetnf < ivec) GO TO 40
 ifset = -1
 70 DO   j = ilist,nlist,2
   iz(j+1) = 1
 END DO
 GO TO 200
 80 DO  i = isetf,nsetf
   k = 0
   diff = 1.e25
   real = z(i)
   DO   j = ilist,nlist,2
     IF (iz(j+1) /= 0) CYCLE
     diff1 =  ABS(z(j) - REAL)
     IF( diff1 >= diff) CYCLE
     diff = diff1
     k = j
   END DO
   IF ( k /= 0)  iz(k+1) = 1
 END DO
 
!     SELECTED FREQUENCIES MARKED FOR OUTPUT
 
 200 nfo =(nlist - ilist +2)/2
 
!     MOVE NEW FREQ  TO UPPER
 
 k=1
 DO  i= ilist,nlist,2
   IF( iz(i+1) == 0) CYCLE
   z(k) = z(i)
   k = k +1
 END DO
 nfn = k-1
 DO   i = ilist,nlist,2
   iz(k) = iz(i+1)
   k = k+1
 END DO
 IF(id == 5) RETURN
 FILE =tol1
 CALL OPEN(*800,tol1,iz(ibuf1),1)
 CALL fname(tol1,ihd)
 CALL WRITE(tol1,ihd,2,0)
 IF(id == 2 .OR. id == 4) GO TO 402
 CALL WRITE(tol1,z,nfn,1)
 401 CONTINUE
 CALL CLOSE(tol1,1)
 mcb(1)= tol1
 mcb(2)= nfn
 CALL wrttrl(mcb )
 IF(id == 2) GO TO 500
 800 RETURN
 
!     COPY OVER CLAMA STUFF
 
 402 CALL WRITE(tol1,0,0,1)
 k = nfn + nfo + 1
 nzx  = ibuf3 -k
 FILE = tol
 CALL gopen(tol,iz(ibuf2),0)
 CALL READ(*910,*920,tol,iz(k),146,1,flag)
 CALL WRITE(tol1,iz(k),146,1)
 m = nfn+1
 n = m+nfo -1
 DO  i = m,n
   CALL READ(*910,*920,tol,iz(k),lw,0,flag)
   IF(iz(i) == 0) CYCLE
   CALL WRITE(tol1,iz(k),lw,0)
 END DO
 CALL CLOSE(tol,1)
 CALL WRITE(tol1,0,0,1)
 GO TO 401
 
!      COPY OVER CASECC
 
 500 CALL gopen(casecc,iz(ibuf1),0)
 CALL gopen(casezz,iz(ibuf2),1)
 m = nfn +1
 n = m+nfo-1
 DO  i = m,n
   CALL READ(*511,*520,casecc,iz(k),nzx,0,flag)
   520 IF(iz(i) == 0) CYCLE
   CALL WRITE(casezz,iz(k),flag,1)
 END DO
 511 CALL CLOSE(casecc,1)
 CALL CLOSE(casezz,1)
 mcb(1) = casecc
 CALL rdtrl(mcb)
 mcb(1) = casezz
 CALL wrttrl(mcb)
 RETURN
 
!     STATIC ANALYSIS
 
 600 CONTINUE
 r = 1.0
 nfo = nfo+ilist
 nlist = nfo-2
 DO  i = ilist,nfo,2
   z(i) = r
   iz(i+1) = 0
   r = r+1.
 END DO
 
!     COPY EDT
 
 CALL OPEN(*670,tol,iz(ibuf1),0)
 CALL OPEN(*670,tol1,iz(ibuf2),1)
 FILE = tol
 CALL fname(tol1,ihd)
 CALL WRITE(tol1,ihd,2,0)
 620 CALL READ(*630,*920,tol,iz(nfo+2),nz,0,flag)
 CALL WRITE(tol1,iz(nfo+2),flag,1)
 GO TO 620
 630 CALL CLOSE(tol,1)
 CALL CLOSE(tol1,1)
 mcb(1) = tol
 CALL rdtrl(mcb)
 mcb(1) = tol1
 CALL wrttrl(mcb)
 670 GO TO 35
 
!     ERROR MESSAGES
 
 900 ip1=-1
 901 CALL mesage(ip1,FILE,NAME)
 910 ip1=-2
 GO TO 901
 920 ip1=-3
 GO TO 901
END SUBROUTINE modac1
