SUBROUTINE rcovb
     
!     RCOVB PERFORMS THE BACK-SUBSTITUTIONS TO OBTAIN THE G-SET
!     DISPLACEMENTS OF A SUBSTRUCTURE WHOSE LEVEL IS LOWER THAN OR
!     EQUAL TO THAT OF THE FINAL SOLUTION STRUCTURE (FSS).
!     FOR EACH SUBSTRUCTURE WHOSE DISPLACEMENTS ARE RECOVERED,
!     AN SOLN ITEM IS CREATED BY EDITING THE SOLN ITEM OF THE FSS.
 
 EXTERNAL         andf
 LOGICAL :: modal
 INTEGER :: mcbtrl(7)  ,dry        ,step       ,fss       ,  &
     rfno       ,uinms      ,schk       ,ua        ,  &
     ssnm1      ,sysbuf     ,rsp        ,rdp       ,  &
     rect       ,upper      ,lower      ,sym       ,  &
     hmcb       ,ubmcb      ,uaomcb     ,uamcb     ,  &
     tflag      ,signab     ,signc      ,scrm      ,  &
     ugv        ,ui(5)      ,scr2       ,scr3      ,  &
     scr5       ,NAME(2)    ,BLANK      ,uvec      ,  &
     pove       ,horg       ,scr1       ,gmask     ,  &
     pao        ,ub         ,iz(1)      ,sofsiz    ,  &
     sof1       ,sof2       ,sof3       ,buf1      ,  &
     buf2       ,rc         ,ssnm(2)    ,eqss      ,  &
     buf(1)     ,andf       ,uimpro     ,energy    ,  &
     rmask      ,FILE       ,rd         ,rdrew     ,  &
     wrt        ,wrtrew     ,rew        ,eofnrw    , buf3       ,buf4
!     INTEGER          SCR6       ,SCR7       ,SRD        ,SWRT
 DOUBLE PRECISION :: dz(1)
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm        ,uwm        ,uim        ,sfm       , swm
 COMMON /BLANK /  dry        ,loop       ,step       ,fss(2)    ,  &
     rfno       ,neigv      ,lui        ,uinms(2,5),  &
     nosort     ,uthres     ,pthres     ,qthres
 COMMON /rcovcr/  icore      ,lcore      ,buf1       ,buf2      ,  &
     buf3       ,buf4       ,sof1       ,sof2      , sof3
 COMMON /rcovcm/  mrecvr     ,ua         ,pa         ,qa        ,  &
     iopt       ,ssnm1(2)   ,energy     ,uimpro    ,  &
     range(2)   ,ireq       ,lreq       ,lbasic
 COMMON /system/  sysbuf     ,nout
 COMMON /names /  rd         ,rdrew      ,wrt        ,wrtrew    ,  &
     rew        ,norew      ,eofnrw     ,rsp       ,  &
     rdp        ,csp        ,cdp        ,square    ,  &
     rect       ,diag       ,upper      ,lower     , sym
 COMMON /mpyadx/  hmcb(7)    ,ubmcb(7)   ,uaomcb(7)  ,uamcb(7)  ,  &
     mpyz       ,tflag      ,signab     ,signc     , mprec      ,scrm
 COMMON /zzzzzz/  z(1)
 EQUIVALENCE      (buf(1)    ,z(1))
 EQUIVALENCE      (z(1)      ,iz(1)      ,dz(1))
 DATA    NAME  /  4HRCOV,4HB          /
 DATA    ugv   ,  scr1,scr2,scr3,scr5 / 106   ,  301, 302, 303, 305  /
 DATA    ui    /  204, 205, 206, 207, 208 /
 DATA    uvec  ,  pove,horg,eqss / 4HUVEC,4HPOVE,4HHORG,4HEQSS /
 DATA    ib    ,  schk / 1,  3   /
 DATA    scr6  ,  scr7,srd,swrt  / 306,307, 1,2 /
 DATA    rmask /  469762048  /
 DATA    gmask /  268435456  /
 DATA    mmask /  134217728  /
 DATA    BLANK /  4H         /
 
!     INITIALIZE
 
 lcorez= korsz(z) - lreq
 sof1  = lcorez - sysbuf + 1
 sof2  = sof1 - sysbuf - 1
 sof3  = sof2 - sysbuf
 buf1  = sof3 - sysbuf
 buf2  = buf1 - sysbuf
 buf3  = buf2 - sysbuf
 buf4  = buf3 - sysbuf
 lcore = buf4 - 1
 IF (lcore <= 0) GO TO 9008
 CALL sofopn (z(sof1),z(sof2),z(sof3))
 ua    = 0
 pao   = 0
 tflag = 0
 signab= 1
 signc = 1
 mprec = 0
 scrm  = scr5
 
!     FIND OUT HOW MANY UI FILES THERE ARE AND WHICH ONES
 
 DO  i = 1,5
   iz(1) = ui(i)
   CALL rdtrl (iz)
   IF (iz(1) < 0) uinms(1,i) = 0
 END DO
 
!     IF UINMS(1,I) = 0         THEN  FILE UI(I) IS PURGED
!     IF UINMS(1,I) = BLANK     THEN  FILE UI(I) IS AVAILABLE AND NOT
!                                     IN USE
!     IF UINMS(1,I) = OTHER     THEN  FILE UI(I) CONTAINS UGV FOR
!                                     SUBSTRUCTURE -OTHER-
 
 ssnm(1) = ssnm1(1)
 ssnm(2) = ssnm1(2)
 
!     IF SSNM IS THE FINAL SOLUTION STRUCTURE (FSS), NO RECOVERY IS
!     NECESSARY.
 
 IF (ssnm(1) /= fss(1) .OR. ssnm(2) /= fss(2)) GO TO 190
 ua = ugv
 GO TO 508
 
!     SEARCH THE SOF FOR A DISPLACEMENT MATRIX OF SSNM OR A HIGHER
!     LEVEL SUBSTRUCTURE FROM WHICH THE REQUESTED DISPLACEMENTS CAN BE
!     RECOVERED
 
 190 jlvl = 1
 200 CALL softrl (ssnm,uvec,mcbtrl)
 rc  = mcbtrl(1)
 IF (rc == 1) GO TO 270
 IF (rc == 2 .AND. dry < 0) GO TO 270
 IF (rc == 3) GO TO 210
 IF (rc == 5) CALL smsg (3,uvec,ssnm)
 IF (rc == 4) GO TO 500
 WRITE (nout,63070) uwm,ssnm1,ssnm
 GO TO 9200
 
!     NO UVEC AT THIS LEVEL.  SAVE SSNM IN A STACK AT TOP OF OPEN CORE
!     AND SEARCH FOR UVEC OF THE NEXT HIGHER LEVEL
 
 210 lastss = 2*jlvl - 1
 iz(lastss  ) = ssnm(1)
 iz(lastss+1) = ssnm(2)
 jlvl = jlvl + 1
 CALL fndnxl (z(lastss),ssnm)
 IF (ssnm(1) /= BLANK) GO TO 230
 WRITE (nout,63060) uwm,iz(lastss),iz(lastss+1)
 GO TO 9200
 230 IF (ssnm(1) /= iz(lastss) .OR. ssnm(2) /= iz(lastss+1)) GO TO 240
 WRITE (nout,63080) uwm,ssnm1,ssnm
 GO TO 9200
 
!     IF SSNM IS NOT THE FSS, LOOK FOR UVEC ON THE SOF.  IF DRY RUN,
!     EXIT.  IF IT IS THE FSS, SET UA=UGV.  IF UGV IS NOT PURGED GO TO
!     BEGIN BACK-SUBSTITUTION.  OTHERWISE, GIVE IT THE SAME TREATMENT
!     AS IF IT WERE NOT THE FSS.
 
 240 IF (ssnm(1) /= fss(1) .OR. ssnm(2) /= fss(2)) GO TO 200
 IF (dry < 0) GO TO 500
 ua = ugv
 mcbtrl(1) = ua
 CALL rdtrl (mcbtrl)
 IF (mcbtrl(1) > 0) GO TO 340
 GO TO 200
 
!     FOUND A UVEC ON SOF FOR THIS LEVEL.  SEE IF IT HAS ALREADY BEEN
!     PUT ON A UI FILE.  (IF DRY RUN, EXIT)
 
 270 IF (dry < 0 .OR. jlvl == 1) GO TO 500
 DO  i = 1,5
   ua = ui(i)
   IF (ssnm(1) == uinms(1,i) .AND. ssnm(2) == uinms(2,i)) GO TO 340
 END DO
 
!     DATA BLANK /4H    /
 
!     IT DOES NOT RESIDE ON ANY UI FILE.  FIND A UI FILE TO USE.
 
 j = 0
 DO  i = 1,5
   IF (uinms(1,i) == 0) CYCLE
   j = j + 1
   IF (uinms(1,i) == BLANK) GO TO 310
 END DO
 GO TO 297
 
!     ALL UI FILES SEEM TO BE IN USE.  DO ANY REALLY EXIST
 
 297 IF (j == 0) GO TO 320
 
!     AT LEAST ONE EXISTS.  RE-USE THE ONE WITH OLDEST DATA
 
 i = lui + 1
 IF (i > 5) i = 1
 j = i
 300 IF (uinms(1,i) /= 0) GO TO 310
 
!     NO FILE THERE.  TRY NEXT ONE.
 
 i = i + 1
 IF (i > 5) i = 1
 IF (i == j) GO TO 320
 GO TO 300
 
!     FOUND A UI FILE TO USE
 
 310 lui = i
 ua  = ui(i)
 uinms(1,i) = ssnm(1)
 uinms(2,i) = ssnm(2)
 GO TO 330
 
!     ALL UI FILES ARE PURGED.  USE SCR1 INSTEAD
 
 320 ua = scr1
 
!     COPY UVEC FROM SOF TO UA
 
 330 CALL mtrxi (ua,ssnm,uvec,0,rc)
 
!     TOP OF BACK-SUBSTITUTION LOOP
 
 340 ub = ua
 uaomcb(1) = 0
 icore = lastss  + 2
 idpcor= icore/2 + 1
 
!     CHECK IF THE EQSS ITEM IS THERE FOR THIS SUBSTRUCTURE
 
 CALL sfetch (z(lastss),eqss,schk,rc)
 IF (rc /= 1) GO TO 6317
 
!     COMPUTE TIME TO RECOVER THIS LEVEL AND CHECK TIME-TO-GO
 
!     (A DETAILED TIME CHECK SHOULD BE CODED LATER.  FOR THE PRESENT,
!     JUST CHECK TO SEE IF TIME HAS RUN OUT NOW.)
 
 CALL tmtogo (i)
 IF (i <= 0) GO TO 6309
 
!     CHECK REMAINING SPACE ON SOF.  FIRST CALCULATE HOW MUCH SPACE
!     THE RECOVERED DISPLACEMENT MATRIX WILL TAKE (ASSUMING IT IS FULL).
 
 mcbtrl(1) = ub
 CALL rdtrl (mcbtrl)
 i = mcbtrl(2)
 
!     NO. OF COLUMNS IN DISPLACEMENT MATRIX IN I
 
 CALL softrl (z(lastss),horg,mcbtrl)
 rc = mcbtrl(1)
 item = horg
 IF (rc > 1) GO TO 6317
 nrow = mcbtrl(3)
 j = i*nrow
 
!     NOW CHECK SPACE
 
 IF (sofsiz(i) < j) GO TO 6310
 
!     CREATE THE SOLUTION ITEM FOR THE RECOVERED SUBSTRUCTURE.
 
 CALL rcovls (z(lastss))
 IF (iopt < 0) GO TO 9000
 
!     FIND A UI FILE FOR DISPLACEMENTS
 
 j = 0
 DO  i = 1,5
   IF (uinms(1,i) == 0) CYCLE
   j = j + 1
   IF (uinms(1,i) == BLANK) GO TO 440
 END DO
 
!     NO UNUSED UI FILES ARE AVAILABLE.  IF TWO OR MORE UI FILES ARE
!     NOT PURGED, USE THE ONE WITH OLDEST DATA.  OTHERWISE, USE SCR2.
!     MAKE SURE WE DON T ASSIGN THE SAME FILE AS THE HIGHER
!     LEVEL DISPLACEMENTS ARE ON (UB)
 
 IF (j < 2) GO TO 450
 i = lui + 1
 IF (i > 5) i = 1
 j = i
 430 IF (uinms(1,i) /= 0 .AND. ui(i) /= ub) GO TO 440
 i = i + 1
 IF (i > 5) i = 1
 IF (i == j) GO TO 450
 GO TO 430
 
!     FOUND A UI FILE
 
 440 lui = i
 ua  = ui(i)
 uinms(1,i) = iz(lastss  )
 uinms(2,i) = iz(lastss+1)
 GO TO 455
 450 ua = scr2
 
!     IF THE RECOVERED SUBSTRUCTURE WAS NOT REDUCED GENERATE THE
!     DISPLACEMENTS DIRECTLY.
!     IF THE SUBSTRUCTURE WAS REDUCED AND THE UIMPROVED FLAG IS SET
!     AND THIS IS A NON-STATICS RUN GENERATE THE IMPROVED DISPLACEMENTS.
!     IF THE SUBSTRUCTURE WAS IN A GUYAN REDUCTION AND THIS IS A
!     STATICS RUN GENERATE THE LOADS ON THE OMMITED POINTS.
 
!     INCLUDE THE CHECK ON THE POVE ITEM ALSO TO BE COMPATABLE WITH
!     PREVIOUS SOFS WITH NO TYPE BITS
 
 455 CALL softrl (z(lastss),pove,mcbtrl)
 ipove = mcbtrl(1)
 CALL fdsub (ssnm,idit)
 rc = 4
 IF (idit < 0) GO TO 6317
 CALL fmdi (idit,imdi)
 modal = .false.
 IF (andf(buf(imdi+ib),mmask) /= 0) modal = .true.
 IF (andf(buf(imdi+ib),rmask) /= 0 .AND. uimpro /= 0 .AND.  &
     rfno > 2) GO TO 470
 IF (andf(buf(imdi+ib),gmask) /= 0 .AND. rfno <= 2) GO TO 480
 IF (andf(buf(imdi+ib),rmask) == 0 .AND. ipove == 1 .AND.  &
     rfno <= 2) GO TO 480
 GO TO 490
 
!     IF THE USER REQUESTED AN IMPROVED VECTOR AND THIS IS A NONSTATICS
!     RUN THEN GENERATE IT.
 
 470 CALL rcovui (ub,z(lastss),modal)
 IF (iopt < 0) GO TO 9000
 GO TO 495
 
!     GENERATE THE LOADS ON THE OMITED POINTS FOR REDUCED SUBSTRUCTURES
!     IF THIS IS A STATICS RUN
 
 480 CALL rcovuo (0,uaomcb(1),z(lastss))
 IF (iopt < 0) GO TO 9000
 
!     MULIPLY AND ADD TO GET DISPLACEMENTS OF LOWER-LIVEL SUBSTRUCTURE.
 
!     COPY H OR G TRANSFORMATION MATRIX TO SCR3
 
 CALL sofopn (z(sof1),z(sof2),z(sof3))
 490 item = horg
 CALL mtrxi (scr3,z(lastss),horg,0,rc)
 IF (rc /= 1) GO TO 6317
 
!     SETUP FOR MPYAD
 
 CALL sofcls
 hmcb(1) = scr3
 ubmcb(1)= ub
 CALL rdtrl (hmcb)
 CALL rdtrl (ubmcb)
 IF (uaomcb(1) /= 0) CALL rdtrl (uaomcb)
 CALL makmcb (uamcb,ua,hmcb(3),rect,ubmcb(5))
 mpyz = lcorez - icore - 7
 CALL mpyad (dz(idpcor),dz(idpcor),dz(idpcor))
 CALL wrttrl (uamcb)
 
!     COPY RECOVERED DISPLACEMENTS TO SOF
 
 495 CALL sofopn (z(sof1),z(sof2),z(sof3))
 CALL mtrxo (ua,z(lastss),uvec,0,rc)
 
!     END OF BACK-SUBSTITUTION LOOP
!     CLOSE AND REOPEN THE SOF TO GET ANY CONTROL BLOCKS WRITTEN TO
!     FILE
 
 CALL sofcls
 CALL sofopn (z(sof1),z(sof2),z(sof3))
 ssnm(1) = iz(lastss)
 ssnm(2) = iz(lastss+1)
 lastss  = lastss - 2
 jlvl    = jlvl - 1
 WRITE (nout,63120) uim,jlvl,ssnm
 IF (jlvl > 1) GO TO 340
 
!     NORMAL COMPLETION OF MODULE EXECUTION
 
 508 CONTINUE
 500 DO  i = 1,5
   IF (uinms(1,i) == 0) uinms(1,i) = BLANK
 END DO
 CALL sofcls
 RETURN
 
!     ERROR PROCESSING
 
 6309 WRITE (nout,63090) sfm,iz(lastss),iz(lastss+1),ssnm,ssnm1
 n = -37
 GO TO 9100
 6310 WRITE (nout,63100) swm,iz(lastss),iz(lastss+1),ssnm,ssnm1
 GO TO 9200
 6317 IF (rc == 2) rc = 3
 CALL smsg (rc-2,item,z(lastss))
 9000 WRITE (nout,63170) swm,ssnm1
 GO TO 9200
 9008 n = 8
 9100 CALL sofcls
 CALL mesage (n,FILE,NAME)
 9200 iopt = -1
 GO TO 500
 
!     FORMAT STATEMENTS
 
 63060 FORMAT (a25,' 6306, ATTEMPT TO RECOVER DISPLACEMENTS FOR NON-',  &
     'EXISTANT SUBSTRUCTURE ',2A4)
 63070 FORMAT (a25,' 6307, WHILE ATTEMPTING TO RECOVER DISPLACEMENTS ',  &
     'FOR SUBSTRUCTURE ',2A4,1H,, /32X,'THE DISPLACEMENTS FOR ',  &
     'SUBSTRUCTURE ',2A4,' WERE FOUND TO EXIST IN DRY RUN ', 'FORM ONLY.')
 63080 FORMAT (a25,' 6308, NO SOLUTION AVAILABLE FROM WHICH DISPLACE',  &
     'MENTS FOR SUBSTRUCTURE ',2A4, /32X,'CAN BE RECOVERED.  ',  &
     'HIGHEST LEVEL SUBSTRUCTURE FOUND WAS ',2A4)
 63090 FORMAT (a25,' 6309, INSUFFICIENT TIME REMAINING TO RECOVER DIS',  &
     'PLACEMENTS OF SUBSTRUCTURE ',2A4, /32X,'FROM THOSE OF ',  &
     'SUBSTRUCTURE ',2A4,'.  (PROCESSING USER RECOVER REQUEST',  &
     /32X,'FOR SUBSTRUCTURE ',2A4,1H))
 63100 FORMAT (a27,' 6310, INSUFFICIENT SPACE ON SOF TO RECOVER DIS',  &
     'PLACEMENTS OF SUBSTRUCTURE ',2A4, /32X,' FROM THOSE OF ',  &
     'SUBSTRUCTURE ',2A4,' WHILE PROCESSING USER RECOVER ',  &
     'REQUEST', /32X,'FOR SUBSTRUCTURE ',2A4)
 63120 FORMAT (a29,' 6312, LEVEL',i4,' DISPLACEMENTS FOR SUBSTRUCTURE ',  &
     2A4, /36X,'HAVE BEEN RECOVERED AND SAVED ON THE SOF.')
 63170 FORMAT (a25,' 6317, RECOVER OF DISPLACEMENTS FOR SUBSTRUCTURE ',  &
     2A4,' ABORTED.')
END SUBROUTINE rcovb
