SUBROUTINE rcovls (lastss)
     
!     THIS ROUTINE CREATES THE SOLN ITEM FOR A LOWER LEVEL SUBSTRUCTURE,
!     LASTSS, BY EDITING THAT OF THE SOLUTION SUBSTRUCTURE FSS.
 
 
 INTEGER, INTENT(IN OUT)                  :: lastss(2)
 INTEGER :: dry        ,step       ,fss        ,rfno       ,  &
     uinms      ,ua         , eqss       ,  &
     rc         ,soln       ,srd        ,swrt       ,  &
     iz(3)      ,eoi        ,eog        ,NAME(2)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /BLANK / dry        ,loop       ,step       ,fss(2)     ,  &
     rfno       ,neigv      ,lui        ,uinms(2,5) ,  &
     nosort     ,uthres     ,pthres     ,qthres
 COMMON /rcovcr/ icore      ,lcore      ,buf1       ,buf2       ,  &
     buf3       ,buf4       ,sof1       ,sof2       , sof3
 COMMON /rcovcm/ mrecvr     ,ua         ,pa         ,qa         ,  &
     iopt       ,rss(2)     ,energy     ,uimpro     ,  &
     range(2)   ,ireq       ,lreq       ,lbasic
 COMMON /system/ sysbuf     ,nout
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE     (z(1),iz(1))
 DATA    NAME  / 4HRCOV, 4HLS  /
 DATA    eqss  , soln  / 4HEQSS,4HSOLN /
 DATA    srd   , swrt  / 1,2   /
 DATA    eog   , eoi   / 2,3   /
 
!     CREATE SOLN ITEM FOR THE RECOVERED SUBSTRUCTURE
 
 IF (rfno == 3) GO TO 490
 
!     OBTAIN LIST OF CONTRIBUTING BASIC SUBSTRUCTURES FROM EQSS.
!     STORE IN OPEN CORE AT ICORE.
 
 CALL sfetch (lastss,eqss,srd,rc)
 CALL suread (z(icore),4,nwds,rc)
 nss = iz(icore+2)
 IF (lcore <= icore+2*nss-1) GO TO 9008
 CALL suread (z(icore),2*nss,nwds,rc)
 
!     CONSTRUCT SOLN GROUP 0 IN OPEN CORE AT IG0.  TWO SLOTS FOR THE
!     NUMBER OF LOADS ON EACH SUBSTRUCTURE.  FIRST IS FOR OLD FSS
!     SOLN, SECOND FOR NEW ONE.
 
 ig0 = icore + 2*nss
 CALL sfetch (fss,soln,srd,rc)
 IF (rc == 1) GO TO 462
 CALL smsg (rc-2,soln,fss)
 GO TO 498
 462 CALL suread (z(ig0),5,nwds,rc)
 isol = iz(ig0+2)
 IF (isol /= rfno) GO TO 6369
 ns = iz(ig0+3)
 nc = iz(ig0+4)
 IF (ig0+4+4*ns > lcore) GO TO 9008
 loop465:  DO  i = 1,ns
   CALL suread (z(ig0+1+4*i),3,nwds,rc)
   iz(ig0+4*i+4) = -65535
   DO  j = 1,nss
     IF (iz(ig0+4*i+1) /= iz(icore+2*j-2)) CYCLE
     IF (iz(ig0+4*i+2) /= iz(icore+2*j-1)) CYCLE
     iz(ig0+4*i+4) = iz(ig0+4*i+3)
     CYCLE loop465
   END DO
 END DO loop465
 IF (rfno == 8 .OR. rfno == 9) GO TO 600
 i = 1
 CALL sjump (i)
 
!     STATICS SOLUTION ITEM
 
!     READ ALL GROUPS OF THE OLD FSS SOLN INTO OPEN CORE AT IGS.
!     AS EACH ONE IS READ, ELIMINATE LOAD VECTORS WHICH DO NOT
!     APPLY TO THE NEW SOLN BY SETTING THEIR LOAD VECTOR
!     NUMBERS TO -65535.
!     UPDATE THE NUMBER OF LOAD VECTORS WHICH DO APPLY.
 
 igs = ig0 + 4*ns + 5
 jgs = igs
 DO  i = 1,nc
   CALL suread (z(jgs),1,nwds,rc)
   n = IABS(iz(jgs))
   IF (jgs+n*2 > lcore) GO TO 9008
   CALL suread (z(jgs+1),-1,nwds,rc)
   nl = 0
   IF (n == 0) GO TO 477
   DO  j = 1,n
     lvn = iz(jgs+2*j-1)
     
!     FIND SUBSTRUCTURE WHERE LVN IS APPLIED FOR FSS SOLN ITEM.
     
     l1 = 0
     l2 = 0
     DO  k = 1,ns
       IF (lvn > l2+iz(ig0+4*k+3)) GO TO 468
       IF (iz(ig0+4*k+4) < 0) EXIT
       lvn = lvn - l1
       nl  = nl + 1
       GO TO 472
       468 IF (iz(ig0+4*k+4) < 0) l1 = l1 + iz(ig0+4*k+3)
       l2 = l2 + iz(ig0+4*k+3)
     END DO
     471 lvn = -65535
     472 iz(jgs+2*j-1) = lvn
   END DO
   IF (iz(jgs) < 0) nl = -nl
   iz(jgs) = nl
   477 jgs = jgs + 2*n + 1
 END DO
 
!     WRITE THE NEW SOLN FOR THE RECOVERED SUBSTRUCTURE ON THE SOF.
!     IN CASE USER FORGOT TO EDIT OUT THIS SOLN FROM A PREVIOUS
!     RUN, DELETE IT TO AVOID LOSING OR SCREWING UP THE RECOVERED
!     DISPLACEMENTS.
 
 CALL DELETE (lastss,soln,rc)
 iz(ig0+3) = nss
 rc = 3
 CALL sfetch (lastss,soln,swrt,rc)
 CALL suwrt (z(ig0),5,1)
 DO  i = 1,ns
   IF (iz(ig0+4*i+4) < 0) CYCLE
   CALL suwrt (z(ig0+4*i+1),2,1)
   CALL suwrt (z(ig0+4*i+4),1,1)
 END DO
 CALL suwrt (0,0,eog)
 jgs = igs
 DO  i = 1,nc
   k  = 0
   nl = iz(jgs)
   jgs= jgs + 1
   CALL suwrt (nl,1,1)
   IF (nl == 0) GO TO 485
   nl = IABS(nl)
   482 IF (iz(jgs) == -65535) GO TO 484
   CALL suwrt (iz(jgs),2,1)
   k = k + 1
   484 jgs = jgs + 2
   IF (k < nl) GO TO 482
   485 IF (iz(jgs) /= -65535) GO TO 486
   jgs = jgs + 2
   GO TO 485
   486 CALL suwrt (0,0,eog)
 END DO
 CALL suwrt (0,0,eoi)
 GO TO 498
 
!     MODAL SOLUTION ITEM
 
!     FOR MODAL COPY THE SOLN UNCHANGED.  IN CASE THE USER FORGOT
!     TO EDIT OUT THIS SOLN FROM A PREVIOUS RUN, DELETE IT TO AVOID
!     LOSING OR SCREWING UP THE RECOVERED DISPLACEMENTS.
 
 490 CALL DELETE (lastss,soln,rc)
 CALL sfetch (fss,soln,srd,rc)
 CALL suread (iz(icore),-1,nwds,rc)
 isol = iz(icore+2)
 IF (isol /= rfno) GO TO 6369
 IF (iz(icore+3) > 0) GO TO 492
 rc = 3
 CALL sfetch (lastss,soln,swrt,rc)
 CALL suwrt (z(icore),4,eog)
 CALL suwrt (0,0,eoi)
 GO TO 498
 492 IF (lcore < icore+7*iz(icore+3)+3) GO TO 9008
 CALL suread (iz(icore+4),-1,nwds,rc)
 rc = 3
 CALL sfetch (lastss,soln,swrt,rc)
 CALL suwrt (z(icore),4,eog)
 CALL suwrt (z(icore+4),7*iz(icore+3),eog)
 CALL suwrt (0,0,eoi)
 GO TO 498
 
!     DYNAMIC SOLUTION ITEM
 
!     READ IN STATIC LOAD SETS
 
 600 incr = 1
 IF (rfno == 8) incr = 2
 igs = ig0 + 4*ns + 5
 CALL suread (z(igs),1,nwds,rc)
 nsl = iz(igs)
 lsl = nsl*incr
 nsll= 0
 IF (nsl == 0) GO TO 660
 IF (igs+nsl > lcore) GO TO 9008
 CALL suread (z(igs+1),nsl,nwds,rc)
 
!     FLAG THOSE STATIC LOAD IDS THAT ARE NOT IN THE LOWER LEVEL
!     SUBSTRUCTURE AND RENUMBER THOSE THAT ARE LEFT
 
 DO  j = 1,nsl
   lvn = iz(igs+j)
   l1  = 0
   l2  = 0
   DO  k = 1,ns
     IF (lvn > l2+iz(ig0+4*k+3)) GO TO 610
     IF (iz(ig0+4*k+4) < 0) EXIT
     lvn  = lvn - l1
     nsll = nsll + 1
     GO TO 640
     610 IF (iz(ig0+4*k+4) < 0) l1 = l1 + iz(ig0+4*k+3)
     l2 = l2 + iz(ig0+4*k+3)
   END DO
   630 lvn = -65535
   640 iz(igs+j) = lvn
 END DO
 
!     COPY THE FREQUENCY OR TIME STEP RECORD INTO CORE
 
 660 i = 1
 CALL sjump (i)
 istep = igs + nsl + 1
 IF (istep+nc > lcore) GO TO 9008
 CALL suread (iz(istep),-1,nwds,rc)
 
!     COPY IN ALL LOAD FACTOR DATA
 
 IF (nsll == 0) GO TO 675
 jgs = istep + nc
 DO  i = 1,nc
   IF (jgs+lsl > lcore) GO TO 9008
   CALL suread (z(jgs),-1,nwds,rc)
   jgs = jgs + lsl
 END DO
 
!     WRITE THE NEW SOLN ITEM FOR THE RECOVERED SUBSTRUCTURE.  IN CASE
!     THE USER FORGOT TO EDIT OUT THIS SOLN FROM A PREVIOUS RUN,
!     DELETE IT TO AVOID LOSING OR SCREWING UP THE RECOVERED
!     DISPLACEMENTS
 
 675 CALL DELETE (lastss,soln,rc)
 rc = 3
 CALL sfetch (lastss,soln,swrt,rc)
 iz(ig0+3) = nss
 CALL suwrt (z(ig0),5,1)
 DO  i = 1,ns
   IF (iz(ig0+4*i+4) < 0) CYCLE
   CALL suwrt (z(ig0+4*i+1),2,1)
   CALL suwrt (z(ig0+4*i+4),1,1)
 END DO
 CALL suwrt (nsll,1,1)
 IF (nsll == 0) GO TO 700
 DO  i = 1,nsl
   IF (z(igs+i) < 0) CYCLE
   CALL suwrt (z(igs+i),1,1)
 END DO
 700 CALL suwrt (0,0,eog)
 
!     COPY THE TIME OR FREQUENCY STEP INFO TO SOF.
 
 CALL suwrt (z(istep),nc,eog)
 
!     COPY LOAD FACTORS FOR EACH STEP TO SOF EDITING OUT THOSE
!     THAT NO LONGER PARTICIAPTE
 
 IF (nsll == 0) GO TO 730
 kgs = istep + nc
 DO  i = 1,nc
   k = 1
   DO  j = 1,nsl
     IF (z(igs+j) < 0) CYCLE
     CALL suwrt (z(kgs+k-1),incr,1)
     k = k + incr
   END DO
   CALL suwrt (0,0,eog)
   kgs = kgs + lsl
 END DO
 
 730 CALL suwrt (0,0,eoi)
 GO TO 498
 
!     NORMAL RETURN
 
 498 RETURN
 
!     ERROR PROCESSING
 
 6369 WRITE (nout,63690) ufm,isol,rfno
 GO TO 9100
 9008 n = 8
 CALL mesage (n,0,NAME)
 9100 iopt = -1
 CALL sofcls
 RETURN
 
 63690 FORMAT (a23,' 6369.  SOLN ITEM HAS INCORRECT RIGID FORMAT NUMBER',  &
     /31X,'SOLUTION RIGID FORMAT WAS',i5,  &
     ' AND CURRENT NASTRAN EXECUTION RIGID FORMAT IS',i5)
END SUBROUTINE rcovls
