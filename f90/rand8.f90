SUBROUTINE rand8(nfreq,npsdl,ntau,xycb,ltab,ifile,psdf,auto,nfile)
     
!     THIS ROUTINE COMPUTES RANDOM RESPONSE FOR COUPLED POWER SPECTRAL
!       DENSITY COEFICIENTS
 
 
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
 REAL :: q(2)
 REAL :: DATA(100)
 
 COMMON /condas/    pi       ,twopi    ,radeg    ,degrad   , s4pisq
 COMMON /system/ sysbuf
 COMMON /zzzzzz/  z(1)
 
 EQUIVALENCE (z(1),iz(1))
 
 DATA  NAME,mcb1,mcb2/4HRAND,4H8   ,14*0/
 DATA ipsdf,iauto /4001,4002/
! *****
!     DEFINITION OF VARIABLES
! *****
!     NFREQ    NUMBER OF FREQUENCIES
!     NPSDL    NUMBER OF PSDL  CARDS
!     NTAU     NUMBER OF TIMES
!     XYCB     DATA BLOCK CONTAINING XY USER REQUESTS
!     LTAB     LENGTH OF CORE USED FOR TABLES BY PRETAB
!     IFILE    ARRAY CONTAINING FILE NAMES FOR SORT 2 INPUT FILES
!     PSDF     OUTPUT FILE FOR POWER SPECTRAL DENSITY FUNCTIONS
!     AUTO     OUTPUT FILE FOR AUTOCORRELATION FUNCTIONS
!     NFILE    LENGTH OF IFILE ARRAY
!     MCB1     TRAILER FOR PSDF
!     MCB2     TRAILER FOR AUTO
!     IPSDF    OFP ID FOR PSDF
!     IAUTO    OFP ID FOR AUTO
!     LCORE    AVAILABLE CORE FOR  LISTS
!     IBUF1    BUFFER POINTERS
!     IBUF2
!     IBUF3
!     ITAU     POINTER TO FIRST TAU-1
!     ISAA     POINTER TO SAB TABLE -1
!     TAU      TIMES FOR AUTOCORRELATION
!     SAB      POWER SPECTRAL DENSITY FACTORS
!     ICORE    POINTER TO FIRST REQUEST-1
!     SYSBUF   LENGTH OF ONE BUFFER
!     NPOINT   TOTAL NUMBER OF REQUESTS
!     NZ       CORE AVAIABLE FOR STORING H VALUES
!     IP       POINTER TO FIRST POINT OF CURRENT CORE LOAD
!     NDONE    NUMBER OF REQUESTS PROCESSED
!     OLDLD    LOAD ID OF OLD LOAD SET
!     NDO      NUMBER POSSIBLE TO DO IN CORE
!     ICS      POINTER TO FIRST H ARRAY
!     NLOAD    NUMBER OF LOADS      PROCESSED ON CURRENT CORE LOAD
!     ICDONE   NUMBER CURRENTLY DONE -- SEVERAL COMP FROM EACH VALUE
!     LOAD     SUBCASE ID FROM INPUT RECORD
!     IF       FORMAT FLAG IF=0  DATA IS REAL/IMAG  IF .NE. 0 MAG/PHASE
!     LEN      LENGTH OF DATA RECORD
!     Q        MEAN  RESPONSE
!     R        AUTOCORRALATION FUNCTION AT TIME TAU
!     IP1      LOCAL POINT POINTER
!     NUNQ     NUMBER OF UNIQUE LOAD ID-S
!     ILOAD    POINTER TO LOAD LIST-1
!     ISJ      POINTER TO SJ ADD AREA-1
!     ICS      H STORAGE -1
 
 
 
! *****
!     CORE LAYOUTDURING EXECUTION
! *****
!     FREQUENCIES   NFREQ OF THEM
!     RANDPS DATA   NPSDL OF THEM  5 WORDS PER CARD
!                   LOAD ID  LOAD ID   X   Y   TABLE
!     TAUS          NTAU OF THEM
!     TABLE DATA    LTAB OF IT
!     S(AB)         NFREQ OF THEM-- THESE ARE REEVALUATED WHEN LOAD CHAN
!     UNIQUE ID-S   NUNQ  OF THEM
!     REQUESTS      NPOINT OF THEM 5 WORDS PER REQUEST
!                   DB   ID   COMP O.P. P/P
!     H-S           LENGTH = 2*NFREQ --REAL+IMAGINARY
!                   NUNQ  H-S PER SET-- NDO SETS
!     SJ COMPUTE    NFREQ OF IT
 
 
!     BUFFERS       3 NEEDED
 
 
 
 
!     INITIALIZE GENERAL VARIABLES--ASSIGN BUFFERS,ETC
 
 mcb1(1)=psdf
 mcb2(1)=auto
 lcore=korsz(z)
 ibuf1=lcore-sysbuf+1
 ibuf2=ibuf1-sysbuf
 ibuf3=ibuf2-sysbuf
 itau=nfreq+5*npsdl
 isaa=ntau+ltab+itau
 lcore=lcore-(isaa+nfreq+3*sysbuf)
 icrq = -lcore
 IF(lcore <= 0) GO TO 980
 
!     BUILD LIST OF UNIQUE LOAD ID-S
!         REPLACE LOAD ID OF PSDL WITH POINTER TO LIST
 
 nunq=0
 iload=isaa+nfreq
 m=iload+1
 k=m-1
 i=nfreq+1
 jj=itau+1
 j=1
 GO TO 4
 
!         SEARCH LIST OF UNIQUE ID-S
 
 5 DO   l=m,k
   IF(iz(i) == iz(l)) GO TO 9
 END DO
 GO TO 4
 
!         NEXT PSDL CARD
 
 2 IF(j == 0) GO TO 7
 i=i+1
 j=0
 GO TO 5
 
!         SAVE LOAD ID
 
 4 k=k+1
 nunq=nunq+1
 iz(k)=iz(i)
 l=k
 
!         REPLACE ID WITH POINTER INTO LIST
 
 9 iz(i)=l-m+1
 GO TO 2
 
!         NEXT PSDL CARD
 
 7 i=i+4
 j=1
 IF(i /= jj) GO TO 5
 
!         COMPUTE MINIMUM CORE
 
 mincr=nunq*nfreq*2+nfreq
 icore=iload+nunq
 lcore=lcore-nunq
 icrq = mincr - lcore
 IF(lcore-mincr <= 0) GO TO 980
 
!         OPEN OUTPUT FILES
 
 CALL gopen(psdf,z(ibuf2),1)
 CALL gopen(auto,z(ibuf3),1)
 
!         BEGIN LOOP ON EACH FILE
 
 DO   i=1,nfile
   
!         BUILD POINT LIST FOR FILE(I)
   
   CALL rand6 (xycb,z(ibuf1),npoint,iz(icore+1),ifile(i),lcore)
   IF(npoint == 0) CYCLE
   nz=lcore-5*npoint
   icrq = -nz
   IF(nz <= 0) GO TO 980
   
!         OPEN INPUT FILE
   
   FILE=ifile(i)
   CALL OPEN(*1000,FILE,z(ibuf1),0)
   ip=icore+1
   ndone=0
   oldld=0
   ics=icore+5*npoint
   llist=5*npoint
   
!         COMPUTE NUMBER OF POINTS TO DO AT SAME TIME
   
   13 ndo = MIN0(npoint-ndone,nz/mincr)
   icrq = MAX0(npoint-ndone,mincr)
   IF(ndo == 0) GO TO 980
   llists = llist
   icdone=0
   ipsave=ip
   nload =0
!         GET READY TO OBTAIN FIRST VALUE
   
   15 CALL rand2 (ifile(i),iz(ip),load,IF,LEN,llist)
   IF(load == 0) GO TO 159
   
!         CHECK FOR NEW LOAD
   
   IF(load == oldld) GO TO 50
   
!         NEW LOAD -- SEE IF WANTED
   
   DO   kk=1,nunq
     l=iload+kk
     IF(load == iz(l)) GO TO 20
   END DO
   
!         REJECT LOAD -- NOT NEEDED
   
   GO TO 15
   
!         GOOD LOAD -- SAVE DATA
   
   20 oldld=load
   
!         BRING DATA INTO KK-TH H SAVE AREA
   
   kk = ics +(kk-1)*nfreq*2
   50 IF(LEN > 100) GO TO 970
   DO   j=1,nfreq
     
!     ACCESS DATA FROM FILE  INTO DATA  ARRAY
     
     CALL rand2a( DATA(1))
     ip1=ip
     ii=icdone
     
!         COMPUTE REAL/IMAG OF CURRENT COMPONENT
     
     52 IF( (LEN-2)/2 >= iz(ip1+2)) GO TO 53
     
!     REQUEST OUT OF RANGE
     
     CALL mesage(52,iz(ip1),iz(ip1+1))
     iz(ip1+2) = (LEN-2)/2
     53 jj = iz(ip1+2) +2
     k=jj+LEN/2-1
     IF ( IF <= 0) GO TO 51
     x=DATA(jj)*COS(degrad*DATA(k))
     DATA(k)=DATA(jj)*SIN(degrad*DATA(k))
     DATA(jj)=x
     51 l=kk+j*2-1+ii*mincr
     z(l)=DATA(jj)
     z(l+1)=DATA(k)
     
!         TEST FOR CORE OVERFLOW
     
     IF(ii == ndo-1) CYCLE
     
!         IS NEXT REQUEST FROM SAME POINT
     
     IF(iz(ip1) /= iz(ip1+5) .OR. iz(ip1+1) /= iz(ip1+6)) CYCLE
     ii=ii+1
     ip1=ip1+5
     GO TO 52
   END DO
   icdone=ii+1
   ip=ip1+5
   llist=llist-5*icdone
   
!         HAVE I DONE ALL REQUESTS (IN CURRENT CORE)
   
   IF(icdone /= ndo) GO TO 15
   
!         HAVE I ADDED IN ALL LOADS
   nload = nload +1
   ip=ipsave
   IF(nload == nunq ) GO TO 100
   
!         START AGAIN ON NEXT LOAD
   llist = ndo*5
   icdone=0
   GO TO 15
   
!         ALL LOADS FOR CURRENT BUNCH DONE
!              COMPUTE SJ-S
   
!              ZERO ALL SJ-S
   
   100 DO   j=1,ndo
     k=ics+j*mincr-nfreq
     DO   l=1,nfreq
       jj=k+l
       z(jj)=0.0
     END DO
   END DO
   
!         FOR EACH PSDL CARD  1. EVALUATE SAB
!              FOR EACH POINT
!                   IN CORE   2. COMPUTE 2*RE(HI*SIJ*HJBAR)
!                             3. ADD TO SJ AT EACH FREQ.
   
   DO   j=1,npsdl
     
!         EVALUATE SAB
     
     two=2.0
     l=nfreq+(j-1)*5
     IF(iz(l+1) == iz(l+2)) two=1.0
     q(1) = z(l+3)
     r=z(l+4)
     DO   k=1,nfreq
       jj=isaa+k
       
       
!                TAB     X    F(X)
       CALL tab (iz(l+5),z(k),z(jj))
       IF(iz(l+5) == 0) z(jj) =1.0
     END DO
     
!         FOR EACH POINT IN CORE
     
     DO   k=1,ndo
       l2=ics+k*mincr-nfreq
       l1=ics+(k-1)*mincr-1-nfreq*2
       DO   m=1,nfreq
         ih1=iz(l+1)*nfreq*2 +l1  +2*m
         ih2=iz(l+2)*nfreq*2 +l1  +2*m
         jj=isaa+m
         isj=l2+m
         z(isj)=z(isj)+z(jj)*two*((z(ih1)*q(1)-z(ih1+1)*r)*z(ih2)  &
             +(z(ih1+1)*q(1)+z(ih1)*r)*z(ih2+1))
       END DO
     END DO
   END DO
   
!         OUTPUT STUFF IN CORE
   
   jj=ip
   j=ndo*5+jj-1
   l=ics-nfreq
   DO   k=jj,j,5
     l=l+mincr
     
!         CONVERT SJ TO ABSOLUTE VALUE
     
     DO   ll=1,nfreq
       kk=l+ll
       z(kk)=ABS(z(kk))
     END DO
     
!         COMPUTE MEAN RESPONSE
     
     CALL rand3 (z(1),z(l+1),q,nfreq)
     IF(iz(k+3) == 2) GO TO 155
     
!         PSDF REQUESTED -- PUT OUT ID
     
     mcb1(7)=mcb1(7)+1
     CALL rand1(psdf,ipsdf,iz(k),iz(k+1),iz(k+4),q)
     
!         PUT OUT DATA RECORDED
     
     DO   ll=1,nfreq
       kk=l+ll
       CALL WRITE (psdf,z(ll),1,0)
       CALL WRITE (psdf,z(kk),1,0)
     END DO
     CALL WRITE (psdf,0,0,1)
     155 IF(iz(k+3) == 1) CYCLE
     
!         AUTO CORRELATION REQUESTED
     
     IF(ntau == 0) CYCLE
     CALL rand1(auto,iauto,iz(k),iz(k+1),iz(k+4),q)
     mcb2(7)=mcb2(7)+1
     
!         PUT OUT DATA RECORD
     
     DO   ll=1,ntau
       kk=itau+ll
       CALL WRITE (auto,z(kk),1,0)
       
!         COMPUTE AUTO
       
       CALL rand4 (z(1),z(l+1),z(kk),r,nfreq)
       CALL WRITE (auto,r,1,0)
     END DO
     CALL WRITE (auto,0,0,1)
   END DO
   
!         END CORE LOAD
   
   CALL REWIND (ifile(i))
   ndone=ndone+ndo
   IF(ndone /= npoint) GO TO 200
   
!         FINISHED WITH FILE
   
   159 CALL CLOSE(ifile(i),1)
   CYCLE
   
!         SPILL ON POINT LISTS -- GO AGAIN
   
   200 oldld=0
   llist=llists-5*ndo
   ip=ipsave+5*ndo
   GO TO 13
   1000 CONTINUE
 END DO
 
!         ALL STUFF DONE -- GET OUT
 
 CALL CLOSE (psdf,1)
 CALL CLOSE (auto,1)
 CALL wrttrl(mcb1)
 CALL wrttrl(mcb2)
 RETURN
 
!         FILE + MISC ERRORS
 
 901 CALL mesage (ip1,FILE,NAME)
 RETURN
 970 ip1=-7
 GO TO 901
 980 ip1=-8
 FILE = icrq
 GO TO 901
END SUBROUTINE rand8
