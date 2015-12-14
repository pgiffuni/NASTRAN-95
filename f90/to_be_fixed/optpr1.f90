SUBROUTINE optpr1
     
!     THIS ROUTINE IS THE DRIVER FOR PROPERTY OPTIMIZATION, PHASE 1.
 
 
!     OPTPR1  MPT,EPT,ECT,DIT,EST/OPTP1/V,N,PRINT/V,N,TSTART/
!                                       V,N,COUNT $
 
!     WHERE PRINT  = OUTPUT, INTEGER = 1
!           TSTART = OUTPUT, INTEGER = TIME AT EXIT OF OPTPR1.
!           COUNT  = OUTPUT, INTEGER =-1 NOT PROPERTY OPTIMIZATION.
!                                    = 1 IS  PROPERTY OPTIMIZATION.
!     CRITERIA FOR OPTIMIZATION
 
!        1. OUTPUT FILE NOT PURGED.
!        2. BULK DATA CARD -POPT IS PRESENT.
!           AFTER THESE TESTS ALL ERRORS ARE FATAL.
 
 
!      SUBROUTINES USED
 
!      OPTP1A - READS ELEMENT DATA INTO CORE (NWDSE PER ELEMENT).
!      OPTP1B - READS PROPERTY IDS INTO CORE AND SETS ELEMENT DATA
!               POINTER (V1) TO ITS LOCATION. (NWDSP PER PROPERTY).
!      OPTP1C - READS DESIGN PROPERTIES INTO CORE.
!      OPTP1D - READS PLIMIT DATA INTO CORE AND SETS PROPERTY DATA
!               POINTER (PLIM) TO ITS LOCATION. (NWDSK PER LIMIT)
 
 
!     LOGICAL         DEBUG
 INTEGER :: dattyp(21),datdty(90),dtyp(90),sysbuf,b2,b1p1,  &
     NAME(2),crew,FILE,ycor,pcor1,ecor1,prcor1,fnam(2),  &
     PRINT,count,poph(2),hpop(2),plmh(2),NONE(2),  &
     ept,ect,dit,est,optp1,outtap,y(1),scrth1,zcor, pcor2,tstart
 REAL :: x(7)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /BLANK / PRINT,tstart,count,skp(2),ycor,b1p1,npow,  &
     nelw,nwdse,nprw,nwdsp,nklw,mpt,ept,ect,dit,est, optp1,scrth1,neltyp,itype(21)
 COMMON /optpw1/ zcor,z(100)
 COMMON /zzzzzz/ core(1)
 COMMON /names / nrd,nrrew,nwrt,nwrew,crew
 COMMON /system/ sysbuf,outtap
 COMMON /gpta1 / ntypes,last,incr,NE(1)
 EQUIVALENCE     (x(1),core(1)), (x(7),y(1))
!     DATA    DEBUG / .FALSE. /
 DATA    poph  , plmh / 404,4, 304,3 /,  NAME / 4H opt,3HPR1 /,  &
     hpop  / 4H   p,4HOPT  /      ,  NONE / 4H (no,4HNE) /,  &
     ltype / 90 /  ,numtyp / 20  /
 
!     NELTYP      = NO. ELEMENT TYPES THAT MAY BE OPTIMIZED
!     LTYPE       = DIMENSION OF DATDTY AND DTYP
!     DATTYP/DTYP = ARRAY TO GIVE RELATIVE LOCATIONS OF ELEMENTS IN
!                   /GPTA1/
 
 DATA    dattyp/ 34, 81, 80, 16, 62, 63, 15, 19, 18, 1,  4,  7,  6,  17,
!             BR  EB  IS  QM  M1  M2  QP  Q1  Q2  RD  SH  TB  T1  T2  &
 73,  9,  8,  3, 64, 83,  0 /
!             T6  TM  TP  TU  Q4  T3
 
!     SETUP DATDYP/DTYP IN ALPHABETICAL ORDER AND IN /GPTA1/ SEQUENCE
 
 DATA    datdty  / 10, 0, 18, 11,  0, 13, 12, 17, 16,  0
!             ELEMENT   RD  2  TU  SH   5  T1  TB  TP  TM  10  &
 ,                 4*0  ,  7,  4, 14,  9,  8,  0
!             ELEMENT   11-14  QP  QM  T2  Q2  Q1  20  &
 ,                 10*0
!             ELEMENT   21-30  &
 ,                 3*0  ,  1,   6*0
!             ELEMENT   31-33  BR   35-40  &
 ,                 10*0
!             ELEMENT   41-50  &
 ,                 10*0
!             ELEMENT   51-60  &
 ,                  0, 5,  6,  19,  6*0
!             ELEMENT   61 M1  M2   Q4  65-70  &
 ,                 2*0,   15,  6*0,   3
!             ELEMENT   71-72  T6  74-79  D8  &
 ,                  2, 0, 20,  7*0 /
!             ELEMENT   EB 82  T3  84-90
 
!     SET UP ELEMENT TYPES
 
 neltyp = numtyp
 DO  i = 1,21
   IF (ntypes > ltype) GO TO 140
   1 itype(i) = dattyp(i)
 END DO
 DO  i = 1,ntypes
   dtyp(i) = datdty(i)
 END DO
 
 
 zcor  = 100
 mpt   = 101
 ept   = 102
 ect   = 103
 dit   = 104
 est   = 105
 optp1 = 201
 scrth1= 301
 
!     STEP 1.  INITIALIZE AND CHECK FOR OUTPUT FILE
 
 count = 0
 PRINT = 1
 CALL fname (optp1,fnam)
 IF (fnam(1) == NONE(1) .AND. fnam(2) == NONE(2)) GO TO 120
 
 b1p1  = korsz(core(1)) - sysbuf
 b2    = b1p1 - sysbuf
 ycor  = b2 - 7
 pcor1 =-1
 ecor1 =-1
 prcor1=-1
 kcor1 =-1
 nwdse = 5
 nwdsp = 6
 npow  = neltyp
 CALL delset
 
!     STEP 2.  FIND POPT CARD
 
 CALL preloc (*120,x(b1p1),mpt)
 CALL locate (*110,x(b1p1),poph,i)
 CALL READ (*10,*30,mpt,x,7,1,nwds)
 
!     ILLEGAL NUMBER OF WORDS
 
 10 CALL page2 (-2)
 WRITE  (outtap,20) sfm,NAME,nwds,hpop
 20 FORMAT (a25,' 2288, ',2A4,'READ INCORRECT NUMBER WORDS (',i2,2A4, 2H).)
 GO TO 80
 
 30 IF (nwds /= 6) GO TO 10
 
!     STEP 2A.  PROCESS PLIMIT CARDS ON SCRATCH FILE
 
 IF (ycor <= 11) GO TO 60
 nklw = 0
 CALL locate (*40,x(b1p1),plmh,i)
 CALL gopen (scrth1,x(b2),nwrew)
 CALL optpx (dtyp)
 CALL CLOSE (scrth1,crew)
 40 CALL CLOSE (mpt,crew)
 IF (nklw    < 0) GO TO 60
 IF (count+1 == 0) GO TO 80
 
!     STEP 3.  LOAD MATERIAL DATA
 
 CALL premat (y(1),y(1),x(b1p1),ycor,mcor,mpt,dit)
 pcor1 = mcor  + 1
 pcor2 = pcor1 + ntypes
 ecor1 = pcor2 + 2*(npow+1)
 ycor  = ycor  - ecor1
 IF (ycor < (nwdse+nwdsp)) GO TO 60
 
!     STEP 4.  READ ELEMENTS INTO CORE
 
 CALL gopen (est,x(b2),0)
 CALL optp1a (y(pcor1),y(pcor2),y(ecor1),dtyp)
 CALL CLOSE (est,crew)
 IF (count+1 == 0) GO TO 80
 IF (nelw    <= 0) GO TO 60
 
!     STEP 5.  READ IN PROPERTIES IDS, SET V1.  SECOND BUFFER NOT NEEDED
 
 prcor1 = ecor1 + nelw
 ycor   = ycor  - nelw + sysbuf
 IF (ycor < nwdsp) GO TO 60
 FILE = ect
 CALL preloc (*90,x(b1p1),ect)
 CALL optp1b (y(pcor1),y(pcor2),y(ecor1),y(prcor1))
 CALL CLOSE (ect,crew)
 IF (count+1 == 0) GO TO 60
 IF (nprw    <= 0) GO TO 80
 
!     STEP 6.  READ PROPERTY DATA INTO CORE
 
 kcor1 = prcor1 + nprw
 ycor  = ycor   - nprw
 
 FILE = ept
 CALL preloc (*90,x(b1p1),ept)
 CALL optp1c (y(pcor1),y(pcor2),y(prcor1))
 CALL CLOSE (ept,crew)
 IF (count+1 == 0) GO TO 80
 
!     STEP 7.  PROCESS PLIMIT CARDS
 
 IF (nklw <= 0) GO TO 50
 IF (ycor < 4) GO TO 60
 CALL gopen (scrth1,x(b1p1),nrrew)
 CALL optp1d (y(pcor2),y(prcor1),y(kcor1))
 CALL CLOSE (scrth1,crew)
 IF (nklw    < 0) GO TO 60
 IF (count+1 == 0) GO TO 80
 
!     STEP 7.  COUNT=0, OUTPUT FILE OPTPR1
 
 50 FILE = optp1
 CALL OPEN  (*90,optp1,x(b1p1),nwrew)
 CALL WRITE (optp1,fnam,2,0)
 CALL WRITE (optp1,x(1),6,1)
 
 CALL WRITE (optp1,y(pcor1),ntypes,0)
 CALL WRITE (optp1,npow,1,0)
 CALL WRITE (optp1,y(pcor2),2*(npow+1),1)
 CALL WRITE (optp1,y(ecor1),nelw,1)
 CALL WRITE (optp1,y(prcor1),nprw,1)
 CALL WRITE (optp1,y(kcor1),nklw,1)
 CALL eof   (optp1)
 j      = 0
 y(j+1) = optp1
 y(j+2) = 0
 y(j+3) = nelw
 y(j+4) = nprw
 y(j+5) = nklw
 y(j+6) = 0
 y(j+7) = ntypes
 CALL wrttrl (y(1))
 CALL CLOSE (optp1,crew)
 GO TO 130
 
!     ERROR MESSAGES - FILE NOT CREATED
 
!     INSUFFICIENT CORE
 
 60 CALL page2 (-3)
 WRITE  (outtap,70) ufm,NAME,b1p1,pcor1,ecor1,prcor1,kcor1
 70 FORMAT (a23,' 2289, ',2A4,'INSUFFICIENT CORE (',i10,2H ), /9X,i9,  &
     ' = MATERIAL',i9,' = POINTERS',i9,' = ELEMENTS',i9, ' = PROPERTIES')
 80 CALL mesage(-61,ept,NAME)
 
!    INPUT FILE PURGED - ILLEGALLY
 
 90 CALL mesage (-1,FILE,NAME)
 
!    OPTPR1 NOT CREATED
 
 110 CALL CLOSE (mpt,crew)
 120 count = -1
 
!     OPTPR1 CREATED
 
 130 CONTINUE
 CALL klock (tstart)
 RETURN
 
!     ERROR MESSAGE
 
 140 WRITE  (outtap,150) sfm
 150 FORMAT (a25,', DATDTY AND DTYP ARRAYS TOO SMALL')
 CALL mesage (-37,0,NAME)
END SUBROUTINE optpr1
