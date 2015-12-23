SUBROUTINE ifp1c (i81,nz)
 
    INTEGER, INTENT(IN OUT)                  :: i81
    INTEGER, INTENT(OUT)                     :: nz
    LOGICAL :: bit64
    INTEGER :: core(1),corey(401),scr1,thru,otpe,exce,BLANK, nifp1c(2)
    COMMON /system/ sysbuf,otpe,nogo,intp,mpcn,spcn,method,loadnn,  &
        nlpp,stftem,ipage,line,tline,maxlin,date(3),tim, iecho,splots,skip(65),intra
    COMMON /ifp1a / scr1,casecc,is,nwpc,ncpw4,nmodes,icc,nset,nsym,  &
        zzzzbb,istr,isub,lencc,iben,equal,IEOR
    COMMON /xifp1 / BLANK,bit64
    COMMON /zzzzzz/ corex(1)
    EQUIVALENCE     (corex(1),corey(1)), (core(1),corey(401))
    DATA   thru   / 4HTHRU/,exce / 4HEXCE/
    DATA   nifp1c / 4H ifp,4H1C  /
 
    i81o = i81
    core(i81+2) = isub
    IF (core(i81+3) /= -1) GO TO 260
    core(i81) = core(i81+4)
    ilset = i81 + 1
    core(ilset) = 0
 
    !     FIND BEGINNING OF SET LIST
 
    i81 = i81 + 5
    IF (core(i81) == IEOR) GO TO 270
    ireal = 0
    IF (core(i81) >    1) GO TO 200
    i81 = i81 + 3
    IF (core(i81) == IEOR) GO TO 270
    iput  = ilset + 2
20  ithru = 0
    iexcpt= 0
30  ASSIGN 20 TO iret
    IF (core(i81) < 0.0) THEN
        GO TO    40
    ELSE IF (core(i81) == 0.0) THEN
        GO TO    60
    ELSE
        GO TO    80
    END IF
40  ithru  = 0
    iexcpt = 0
50  IF (IABS(core(i81)) /= 1) ireal = 1
    core(iput) = core(i81+1)
    ibk1 = IABS(core(i81+1))
    i81  = i81  + 2
    iput = iput + 1
    core(ilset) = core(ilset) + 1
    GO TO 30
 
    !     CONTINUATION CARD
 
    ! ... ALLOW ON-LINE READ IF INTRA IS .GT. ZERO, SET BY ONLINS
 
60  IF (intra <= 0) GO TO 65
    CALL xread (*240,core(1))
    icc = icc + 1
    GO TO 67
65  CALL READ (*240,*240,scr1,core(1),nwpc,0,flag)
    WRITE (otpe,250) icc,(core(i),i=1,nwpc)
    icc  = icc  + 1
    line = line + 1
    IF (line >= nlpp) CALL page
67  i81 = iput
    nz  = nz - core(ilset)
    CALL xrcard (core(i81),nz,core(1))
    GO TO iret, (20,120)
 
    !     END OF RECORD
 
70  i81 = iput
    IF (core(ilset)-1 < 0.0) THEN
        GO TO   200
    ELSE IF (core(ilset)-1 == 0.0) THEN
        GO TO   230
    END IF
71 CONTINUE
   IF (ireal == 1) GO TO 230
 
   !     SORT LIST
 
   iset = core(ilset)
   CALL ifp1s (core(ilset+2),core(i81),core(ilset))
 
   !     CORRECT FOR DELETIONS
 
   i81 = i81 + core(ilset) - iset
   GO TO 230
 
   !     THRU AND EXCEPT
 
80 IF (core(i81) == IEOR) GO TO 70
   IF (ireal == 1) CALL ifp1d(-622)
   IF (bit64) CALL mvbits (BLANK,0,32,core(i81+1),0)
   IF (core(i81+1)  /= thru) GO TO 90
   IF (core(ilset)  == 0) GO TO 200
   IF (core(iput-1) < 0) GO TO 280
   i81 = i81 +3
   IF (core(i81) == IEOR) GO TO 270
   ibk  = ibk1
   ifwd = core(i81+1)
   ifwd1= ifwd
   IF (ibk >= ifwd) GO TO 200
   ithru = 1
   !     TEST FOR DEGENERATE THRU INTERVAL
   IF (ifwd-ibk == 1) GO TO 50
   core(i81+1) = -core(i81+1)
   GO TO 50
 
   !     EXCEPT
 
90 IF (core(i81+1) /= exce) GO TO 200
   IF (ithru == 1) GO TO 110
 
   !     EXCEPT WITHOUT THRU
 
   CALL ifp1d (-613)
   GO TO 220
 
   !     PROCESS EXCEPT CANDIDATES
 
110 i81 = i81 + 3
   IF (core(i81) == IEOR) GO TO 270
   IF (iexcpt == 1) GO TO 280
   iexcpt = 1
   jexcpt = 0
120 ASSIGN 120 TO iret
   IF (core(i81) < 0.0) THEN
       GO TO   130
   ELSE IF (core(i81) == 0.0) THEN
       GO TO    60
   ELSE
       GO TO    80
   END IF
130 IF (core(i81+1) > ifwd1) GO TO 20
   IF (core(i81+1) <   ibk) GO TO 200
   IF (core(i81+1) <= core(i81-1) .AND. jexcpt == 1 .AND.  &
       (core(i81+2) <= 0 .OR. core(i81+2) == IEOR)) GO TO 160
   jexcpt = 1
   IF (core(i81+1) ==   ibk) GO TO 290
   IF (core(i81+1) ==  ifwd) GO TO 300
   IF (core(i81+1) == ifwd1) GO TO 310
   IF (core(i81+1)-1 == ibk) GO TO 140
   IF (core(i81+1)+1 == ifwd) GO TO 180
   !     EXCEPT IN MIDDLE OF INTERVAL
   core(iput-1) = -core(i81+1) + 1
   iacip = IABS(core(iput-1))
   IF (iacip-ibk == 1) core(iput-1) = iacip
   core(iput  ) = core(i81+1)+1
   core(iput+1) = -ifwd
   IF (ifwd-core(iput) == 1) core(iput+1) = ifwd
   ibk = core(iput)
   i81 = i81  + 2
   iput= iput + 2
   core(ilset) = core(ilset) + 2
   GO TO 120
   !     EXCEPT ADJACENT TO BOTTOM OF INTERVAL
140 il1 = core(iput-1)
   ibk = ibk + 2
   core(iput-1) = ibk
   ial1 = IABS(il1)
   IF (ial1-ibk == 1) il1 = ial1
   core(iput) = il1
   IF (ibk /= ial1) GO TO 150
   ibk  = 0
   ifwd = 0
   i81  = i81 + 2
   GO TO 120
150 iput = iput + 1
   i81  = i81  + 2
   core(ilset) = core(ilset) + 1
   GO TO 120
160 CALL ifp1d (-626)
   i81 = i81 + 2
   GO TO 120
   !     EXCEPT ADJACENT TO TOP OF INTERVAL
180 core(iput) = IABS(core(iput-1))
   ifwd = ifwd - 2
   core(iput-1) = -ifwd
   IF (ifwd-ibk == 1) core(iput-1) = ifwd
   GO TO 150
 
   !     FOULED UP SET
 
200 CALL ifp1d (-614)
220 i81  = i81o
   nset = nset - 1
230 RETURN
240 CALL mesage (-1,scr1,nifp1c)
   GO TO 240
250 FORMAT (11X,i8,6X,20A4)
 
   !     NO NAME FOR SET
 
260 CALL ifp1d (-615)
   GO TO 220
 
   !     UNEXPECTED END OF RECORD
 
270 CALL ifp1d (-623)
   GO TO 220
 
   !     EXCEPT FOLLOWED BY THRU
 
280 CALL ifp1d (-616)
   GO TO 220
 
   !     EXCEPTING BEGINNING OF INTERVAL
 
290 ibk = ibk + 1
   core(iput-2) = ibk
   i81 = i81 + 2
   IF (ifwd-ibk == 1) core(iput-1) = ifwd
   IF (ibk /= ifwd) GO TO 120
   iput = iput - 1
   core(ilset) = core(ilset) - 1
   ibk  = 0
   ifwd = 0
   GO TO 120
 
   !     EXCEPT END OF INTERVAL
 
300 ifwd = ifwd - 1
   core(iput-1) = -ifwd
   i81 = i81 + 2
   IF (ifwd-ibk == 1) core(iput-1) = ifwd
   IF (ibk /= ifwd) GO TO 20
   iput = iput - 1
   core(ilset) = core(ilset) - 1
   GO TO 20
 
   !     EXCEPT PAST OLD END OF INTERVAL
 
310 i81  = i81 + 2
   iput = iput- 1
   core(ilset) = core(ilset) - 1
   GO TO 20

   END SUBROUTINE ifp1c
