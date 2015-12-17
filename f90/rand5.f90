SUBROUTINE rand5 (nfreq,npsdl,ntau,xycb,ltab,ifile,psdf,auto,  &
        nfile)
     
!     THIS ROUTINE COMPUTES RANDOM RESPONSE FOR UNCOUPLED POWER SPECTRAL
!     DENSITY COEFFICIENTS
 
 
 INTEGER, INTENT(IN)                      :: nfreq
 INTEGER, INTENT(IN)                      :: npsdl
 INTEGER, INTENT(IN)                      :: ntau
 INTEGER, INTENT(IN OUT)                  :: xycb
 INTEGER, INTENT(IN)                      :: ltab
 INTEGER, INTENT(IN)                      :: ifile(1)
 INTEGER, INTENT(IN)                      :: psdf
 INTEGER, INTENT(IN)                      :: auto
 INTEGER, INTENT(IN)                      :: nfile
 INTEGER :: iz(1),sysbuf,FILE, NAME(2), mcb1(7),mcb2(7),oldld
 REAL :: q(2),DATA(100)
 
 COMMON /system/ sysbuf
 COMMON /zzzzzz/ z(1)
 
 EQUIVALENCE  (z(1),iz(1))
 
 DATA    NAME,mcb1,mcb2 /4HRAND,4H5   ,14*0/
 DATA    ipsdf,iauto    /  4001,4002       /
! *****
!     DEFINITION OF VARIABLES
! *****
!     NFREQ    NUMBER OF FREQUENCIES
!     NPSDL    NUMBER OF SUBCASES ON PSDL CARDS
!     NTAU     NUMBER OF TIMES
!     XYCB     DATA BLOCK CONTAINING XY USER REQUESTS
!     LTAB     LENGTH OF CORE  USED FOR TABLES BY PRETAB
!     IFILE    ARRY CONTAINING FILE NAMES FOR SORT 2 INPUT FILES
!     PSDF     OUTPUT FILE FOR POWER SPECTRAL DENSITY FUNCTIONS
!     AUTO     OUTPUT FILE FOR AUTOCORRELATION FUNCTIONS
!     NFILE    LENGTH OF IFILE ARRAY
!     MCB1     TRAILER FOR PSDF
!     MCB2     TRAILER FOR AUTO
!     IPSDF    OFP ID FOR  PSDF
!     IAUTO    OFP ID FOR  AUTO
!     LCORE    AVAIABLE CORE FOR ANY LIST
!     IBUF1    BUFFER POINTERS
!     IBUF2
!     IBUF3
!     ITAU     POINTER TO FIRST TAU -1
!     ISAA     POINTER TO FIRST S(AA) -1
!     TAU      TIMES FOR AUTTOCORRELATION
!     SAA      POWER SPECTRAL DENSITY FACTORS
!     ICORE    POINTER  TO FIRST REQUEST -1
!     SYSBUF   LENGTH OF ONE BUFFER
!     NPOINT   NUMBER OF REQUESTS
!     NZ       CORE AVAILABLE FOR STORING PSDF-S
!     IP       POINTER TO FIRST POINT OF CURRENT CORE LOAD
!     NDONE    NUMBER OF REQUESTS PROCESSED
!     OLDLD    LOAD ID OF OLD LOAD SET
!     NDO      NUMBER POSSIBLE TO DO IN CORE
!     ICS      POINTER TO FIRST PSDF VECTOR
!     NLOAD    NUMBER OF PSDL CARDS PROCESSED
!     ICDONE   NUMBER CURRENTLY DONE-- SEVERAL COMP FROM EACH VALUE
!     LOAD     SUBCASE ID FROM INPUT RECORD
!     IF       FORMAT FLAG  IF =0  DATA IS REA/IMAG IF.NE.0 MAG/PHASE
!     LEN      LENGTH OF DATA RECORD
!     Q        MEAN RESPONSE
!     R        AUTO CORRALATION FUNCTION AT TIME TAU
!     IP1      LOCAL POINT POINTER
 
 
! *****
!     CORE LAYOUT DURING EXECUTION
! *****
!     FREQUENCIES   NFREQ OF THEM
!     RANDPS DATA   NPSDL OF THEM  5 WORDS PER CARD
!                   LOAD ID   LOAD ID   X    Y=0. TABLE
!     TAUS          NTAU OF THEM
!     TABLE DATA    LTAB OF IT
!     S(AA)         NFREQ OF THEM  THESE ARE  REEVALUATED WHEN LOAD CHAG
!     REQUESTS      NPOINT OF THEM 5 WORDS PER REQUEST
!                   D.B. ID   COMP O.P. P/P
!     S(J,A)        NO  DO OF THEM      LENGTH = NFREQ
 
 
!     BUFFERS       3 NEEDED
 
 
!     INITALIZE GENERAL VARIABLES, ASSIGN BUFFERS  ETC
 
 mcb1(1) = psdf
 mcb2(1) = auto
 lcore = korsz(z)
 ibuf1 = lcore -sysbuf+1
 ibuf2 = ibuf1 -sysbuf
 ibuf3 = ibuf2 -sysbuf
 itau  = nfreq +5*npsdl
 isaa  = ntau  +ltab+itau
 icore = isaa  +nfreq
 lcore = lcore -icore-3*sysbuf
 icrq  =-lcore
 IF (lcore <= 0) GO TO 980
 
!     OPEN OUTPUT FILES
 
 CALL gopen (psdf,z(ibuf2),1)
 CALL gopen (auto,z(ibuf3),1)
 
!     BEGIN LOOP ON EACH FILE
 
 DO  i = 1,nfile
   
!     BUILD POINT LIST FOR FILE(I)
   
   CALL rand6 (xycb,z(ibuf1),npoint,iz(icore+1),ifile(i),lcore)
   IF (npoint == 0) CYCLE
   nz   = lcore -5*npoint
   icrq =-nz
   IF (nz <= 0) GO TO 980
   
!     OPEN INPUT FILE
   
   FILE  = ifile(i)
   CALL OPEN (*1000,FILE,z(ibuf1),0)
   ip    = icore +1
   ndone = 0
   oldld = 0
   ics   = icore +5*npoint +1
   llist = 5*npoint
   ips   = ip
   llists= llist
   13 ndo   = MIN0(npoint-ndone,nz/nfreq)
   icrq  = MAX0(npoint-ndone,nfreq)
   IF (ndo == 0) GO TO 980
   nload = 0
   
!     ZERO CORE
   
   jj = ics + ndo*nfreq-1
   DO  k = ics,jj
     z(k) = 0.0
   END DO
   icdone = 0
   
!     GET READY TO OBTAIN FIRST VALUE
   
   15 CALL rand2 (ifile(i),iz(ip),load,IF,LEN,llist)
   
!     CHECK FOR NEW LOAD
   
   IF (load == 0) IF (nload-npsdl) 111,100,111
   IF (load == oldld) GO TO 50
   
!     NEW LOAD --EVALUATE S(AA) FUNCTIONS FOR THIS LOAD
   
   j  = nfreq +1
   jj = itau
   DO  k = j,jj,5
     IF (iz(k) == load) GO TO 20
   END DO
   
!     LOAD NOT NEEDED --REJECT
   
   GO TO 15
   
!     GOOD LOAD --EVALUATE
   
   20 oldld = load
   nload = nload +1
   DO  j = 1,nfreq
     jj = isaa +j
     
!                TAB      X    F(X)
     CALL tab (iz(k+4),z(j),z(jj))
     IF (iz(k+4) == 0) z(jj) = 1.0
     z(jj) = z(jj)*z(k+2)
   END DO
   
!     BRING IN DATA
   
   50 IF (LEN > 100) GO TO 970
   DO  j = 1,nfreq
     
!     ACCESS DATA FROM FILE INTO DATA ARRAY
     
     CALL rand2a (DATA(1))
     ip1 = ip
     ii  = icdone
     ll  = isaa +j
     
!     COMPUTE  MAGNITUDE         OF CURRENT COMPONENT
     
     52 IF ((LEN-2)/2 >= iz(ip1+2)) GO TO 53
     
!     REQUEST OUT OF RANGE
     
     CALL mesage (52,iz(ip1),iz(ip1+1))
     iz(ip1+2) = (LEN-2)/2
     53 jj = iz(ip1+2) +2
     IF (IF /= 0) GO TO 51
     
!     REAL + IMAGINARY
     
     k  = jj + LEN/2 -1
     DATA(jj) = SQRT(DATA(jj)*DATA(jj) + DATA(k)*DATA(k))
     
!     COMPUTE POWER SPECTRAL DENSITY FUNCTION
     
     51 k = ics + ii*nfreq-1 +j
     z(k) = z(k) + DATA(jj)*z(ll)*DATA(jj)
     IF (ii == ndo-1) CYCLE
     
!     IS NEXT REQUEST FROM SAME POINT
     
     IF (iz(ip1) /= iz(ip1+5) .OR. iz(ip1+1) /= iz(ip1+6)) CYCLE
     ii = ii +1
     ip1= ip1+ 5
     GO TO 52
   END DO
   llist  = llist - 5*(ii-icdone+1)
   icdone = ii +1
   ip = ip1 +5
!     HAVE I DONE ALL REQUEST(IN CURRENT CORE)
   
   IF (icdone /= ndo) GO TO 15
   
!     HAVE I ADDED IN ALL LOADS
   
   ip = ips
   IF (nload == npsdl) GO TO 100
   
!     START AGAIN ON NEXT LOAD
   
   llist  = llists
   icdone =  0
   GO TO 15
   
!     ALL LOADS FOR CURRENT BUNCH DONE
   
   100 jj = ip
   j  = ndo* 5 +jj-1
   l  = ics - nfreq
   DO  k = jj,j,5
     l = l+ nfreq
     
!     COMPUTE MEAN RESPONSE   Q
     
     CALL rand3 (z(1),z(l),q,nfreq)
     IF (iz(k+3) == 2) GO TO 105
     
!     PSDF REQUESTED
     
!     PUT OUT ID RECORD
     
     mcb1(7) = mcb1(7) +1
     CALL rand1 (psdf,ipsdf,iz(k),iz(k+1),iz(k+4),q)
     
!     PUT OUT DATA RECORD
     
     DO  ll = 1,nfreq
       kk = l +ll -1
       CALL WRITE (psdf,z(ll),1,0)
       CALL WRITE (psdf,z(kk),1,0)
     END DO
     CALL WRITE (psdf,0,0,1)
     105 IF (iz(k+3) == 1) CYCLE
     
!     AUTOCORRELATION REQUESTED
     
     IF (ntau == 0) CYCLE
     CALL rand1 (auto,iauto,iz(k),iz(k+1),iz(k+4),q)
     mcb2(7) = mcb2(7)+1
     
!     PUT OUT DATA RECORD
     
     DO  ll = 1,ntau
       kk = itau + ll
       CALL WRITE (auto,z(kk),1,0)
       
!     COMPUTE AUTO
       
       CALL rand4 (z(1),z(l),z(kk),r,nfreq)
       CALL WRITE (auto,r,1,0)
     END DO
     CALL WRITE (auto,0,0,1)
   END DO
   CALL REWIND (ifile(i))
   ndone = ndone +ndo
   IF (ndone /= npoint) GO TO 200
   111 CALL CLOSE (ifile(i),1)
   CYCLE
   
!     SPILL ON POINT LISTS --GO AGAIN
   
   200 oldld = 0
   ip  = ip + 5*ndo
   ips = ip
   llists = llists-5*ndo
   GO TO 13
   1000 CONTINUE
 END DO
 
!     ALL STUFF DONE --GET OUT
 
 CALL CLOSE (psdf,1)
 CALL CLOSE (auto,1)
 CALL wrttrl (mcb1)
 CALL wrttrl (mcb2)
 RETURN
 
!     FILE + MISC ERRORS
 
 901 CALL mesage (ip1,FILE,NAME)
 RETURN
 970 ip1 = -7
 GO TO 901
 980 ip1 = -8
 FILE= icrq
 GO TO 901
END SUBROUTINE rand5
