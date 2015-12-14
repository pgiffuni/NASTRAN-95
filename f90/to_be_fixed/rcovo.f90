SUBROUTINE rcovo
     
!     RCOVO READS THE CASESS RECOVER RECORD AND PROCESSES ANY
!     OUTPUT REQUESTS FOR THE CURRENT SAVE OR PRINT REQUEST
 
!     THE OUTPUT REQUESTS ARE STORED AT THE BOTTOM OF OPEN CORE IN
!     A TABLE WITH THE FOLLOWING FORM
 
!     BUF(IREQ) -  UVEC  -I
!                  PVEC   I- NONZERO IF ANY REQUEST PRESENT
!                  QVEC  -I
!                  NO. OF POINTS
!                  NO. OF BASICS
!                  BASIC NAME(2) -I
!                  DISP SET       I
!                  OLOAD SET      I
!                  SPCF SET       I - REPEATED FOR EACH BASIC
!                  SUBCASES SET   I   SUBSTRUCTURE
!                  MODES SET      I
!                  RANGE(2)       I
!                  VELO SET       I
!                  ACCE SET       I
!                  STEPS SET      I
!                  GRID OR MODAL -I
 
 EXTERNAL        lshift     ,andf
 LOGICAL :: basic      ,mrecvr
 INTEGER :: rss        ,buf(1)     ,step       ,casess     ,  &
     eqss       ,z          ,buf1       ,sysbuf     ,  &
     sof1       ,sof2       ,sof3       ,recovr     ,  &
     SAVE       ,PRINT      ,srd        ,fss(2)     ,  &
     rd         ,subnam(2)  ,rdrew      ,wrt        ,  &
     wrtrew     ,rew        ,norew      ,rc         ,  &
     REC(3)     ,subc       ,subs       ,mode       ,  &
     all        ,NONE       ,comds(13)  ,energy     ,  &
     uimpro     ,rfno       ,time       ,freq       , andf       ,mrecov
 REAL :: rbuf(1)    ,rrec(3)
 CHARACTER (LEN=1) ::
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm swm*27
 COMMON /xmssg / ufm        ,uwm        ,uim        ,sfm        , swm
 COMMON /BLANK / dry        ,loop       ,step       ,fss        ,  &
     rfno       ,neigv      ,lui        ,uinms(2,5) ,  &
     nosort     ,uthres     ,pthres     ,qthres
 COMMON /rcovcr/ icore      ,lcore      ,buf1       ,buf2       ,  &
     buf3       ,buf4       ,sof1       ,sof2       , sof3
 COMMON /rcovcm/ mrecvr     ,ua         ,pa         ,qa         ,  &
     iopt       ,rss(2)     ,energy     ,uimpro     ,  &
     rang(2)    ,ireq       ,lreq       ,lbasic
 COMMON /zzzzzz/ z(1)
 COMMON /system/ sysbuf     ,nout
 COMMON /names / rd         ,rdrew      ,wrt        ,wrtrew     ,  &
     rew        ,norew      ,eofnrw
 EQUIVALENCE     (buf(1)    ,z(1))
 EQUIVALENCE     (buf(1)    ,rbuf(1))   ,(REC(1)    ,rrec(1))   ,  &
     (iset      ,rset)
 DATA   casess / 101   /, subnam / 4HRCOV,4HO   /
 DATA   eqss   / 4HEQSS/
 DATA   recovr / 4HRECO/, mrecov / 4HMREC       /
 DATA   PRINT  / 4HPRIN/, SAVE   / 4HSAVE       /
 DATA   srd    / 1     /
 DATA   ll     / 2     /
 DATA   ncomds / 13    /
 DATA   comds  / 4HDISP,4HOLOA,4HSPCF,4HMODE,4HRANG,4HSUBC,4HSORT,  &
     4HBASI,4HVELO,4HACCE,4HENER,4HUIMP,4HSTEP       /
 DATA   subc   , subs,   mode,   all,    NONE,   time,   freq    /  &
     4HSUBC , 4HSUBS, 4HMODE, 4HALL , 4HNONE, 4HTIME, 4HFREQ  /
 
!     SET UP BUFFERS
 
 sof1  = 1
 sof2  = sof1 + sysbuf
 sof3  = sof2 + sysbuf + 1
 buf1  = sof3 + sysbuf
 icore = buf1 + sysbuf
 lcore = korsz(z(1)) - icore + 1
 IF (lcore <= 0) GO TO 9008
 
!     FIND RECOVER RECORD IN CASESS
 
 CALL gopen (casess,z(buf1),rdrew)
 IF (step == 1) GO TO 20
 DO  i = 2,step
   CALL fwdrec (*9002,casess)
 END DO
 20 CALL fread (casess,REC,2,0)
 IF (REC(1) /= recovr .AND. REC(1) /= mrecov) GO TO 6305
 mrecvr = .false.
 IF (REC(1) == mrecov) mrecvr = .true.
 
!     GET PRINT OR SAVE OPTION FOR THIS PASS
 
 i = 0
 30 CALL READ (*9002,*600,casess,REC,3,0,nwds)
 IF (REC(1) /= PRINT .AND. REC(1) /= SAVE) GO TO 30
 IF (loop == i) GO TO 40
 i = i + 1
 GO TO 30
 
!     GET NAME OF SUBSTRUCTURE TO BE OPERATED ON
 
 40 rss(1) = REC(2)
 rss(2) = REC(3)
 loop   = loop + 1
 IF (REC(1) == SAVE) GO TO 700
 iopt   = 1
 
!     OPEN SOF AND FETCH EQSS FOR SUBSTRUCTURE TO BE PRINTED
 
 
 CALL sofopn (z(sof1),z(sof2),z(sof3))
 CALL sfetch (rss,eqss,srd,rc)
 SELECT CASE ( rc )
   CASE (    1)
     GO TO 60
   CASE (    2)
     GO TO 50
   CASE (    3)
     GO TO 50
   CASE (    4)
     GO TO 6306
   CASE (    5)
     GO TO 50
 END SELECT
 
!     FETCH ON EQSS WAS UNSUCCESSFUL
 
 50 IF (rc == 2) rc = 3
 CALL smsg (rc-2,eqss,rss)
 GO TO 800
 
!     READ GROUP 0 OF EQSS INTO CORE
 
 60 CALL suread (z(icore),lcore,nwds,rc)
 SELECT CASE ( rc )
   CASE (    1)
     GO TO 9008
   CASE (    2)
     GO TO 65
   CASE (    3)
     GO TO 62
   CASE (    4)
     GO TO 800
 END SELECT
 62 CALL smsg (7,eqss,rss)
 GO TO 800
 
!     DETERMINE SIZE OF OUTPUT REQUEST BLOCK AND ALLOCATE SPACE
!     AT BOTTOM OF OPEN CORE
 
 65 nbs  = z(icore+2)
 np   = z(icore+3)
 lbasic = 13
 lreq = 5 + lbasic*nbs + 2
 IF (lreq > lcore-nwds) GO TO 9008
 nreq = korsz (buf(1))
 ireq = nreq - lreq + 1
 DO  i = ireq,nreq
   buf(i) = 0
 END DO
 
!     MOVE NAMES OF BASICS INTO OUTPUT AREA
 
 buf(ireq+3) = np
 buf(ireq+4) = nbs
 DO  i = 1,nbs
   i1 = ireq  + (i-1)*lbasic + 5
   i2 = icore + (i-1)*2+4
   buf(i1  ) = z(i2  )
   buf(i1+1) = z(i2+1)
 END DO
 
!     INSERT DEFAULTS INTO OUTPUT BLOCK
 
!     MODES    = ALL
!     SUBCASES = ALL
!     RANGE    = -1.0E+35,1.0E+35
!     STEPS    = ALL
 
 energy  = 0
 uimpro  = 0
 rang(1) = -1.0E+35
 rang(2) =  1.0E+35
 buf(ireq  ) = -2
 buf(ireq+1) = -2
 buf(ireq+2) = -2
 DO  i = 1,nbs
   i1 = ireq + (i-1)*lbasic + 5
   buf(i1+ 2) = -2
   buf(i1+ 3) = -2
   buf(i1+ 4) = -2
   buf(i1+ 5) = -1
   buf(i1+ 6) = -1
   rbuf(i1+7) = -1.0E+35
   rbuf(i1+8) =  1.0E+35
   buf(i1+ 9) = -2
   buf(i1+10) = -2
   buf(i1+11) = -1
 END DO
 
!     READ NEXT COMMAND AND PROCESS OUTPUT REQUEST
 
 nss1   = 1
 nss2   = nbs
 irange = 0
 basic  = .false.
 90 CALL READ (*9002,*500,casess,REC,3,0,nwds)
 IF (REC(1) == PRINT .OR. REC(1) == SAVE) GO TO 510
 DO  i = 1,ncomds
   IF (REC(1) == comds(i)) GO TO 110
 END DO
 GO TO 90
 110 CONTINUE
 SELECT CASE ( i )
   CASE (    1)
     GO TO 120
   CASE (    2)
     GO TO 130
   CASE (    3)
     GO TO 140
   CASE (    4)
     GO TO 150
   CASE (    5)
     GO TO 160
   CASE (    6)
     GO TO 165
   CASE (    7)
     GO TO 170
   CASE (    8)
     GO TO 190
   CASE (    9)
     GO TO 230
   CASE (   10)
     GO TO 240
   CASE (   11)
     GO TO 250
   CASE (   12)
     GO TO 260
   CASE (   13)
     GO TO 270
 END SELECT
 
!     DISP REQUEST
 
 120 IF (REC(2) /= NONE) buf(ireq) = 1
 iloc = 2
 GO TO 400
 
!     OLOAD REQUEST
 
 130 IF (REC(2) /= NONE) buf(ireq+1) = 1
 IF (REC(2) == NONE .AND. .NOT.basic) buf(ireq+1) = 0
 iloc = 3
 GO TO 400
 
!     SPCF REQUEST
 
 140 IF (REC(2) /= NONE) buf(ireq+2) = 1
 IF (REC(2) == NONE .AND. .NOT.basic) buf(ireq+2) = 0
 iloc = 4
 GO TO 400
 
!     MODES REQUEST
 
 150 iloc = 6
 GO TO 400
 
!     RANGE REQUEST (IF BEFORE A BASIC COMMAND SAVE IT FOR ENERGY
!                    PROCESSING ALSO)
 
 160 iloc = 7
 IF (MOD(irange,2) == 1) iloc = 8
 irange = irange + 1
 IF (REC(2) /= -2 .AND. REC(3) /= 0) GO TO 450
 IF (basic) GO TO 410
 rang(iloc-6) = rrec(3)
 GO TO 410
 
!     SUBCASES REQUEST
 
 165 iloc = 5
 GO TO 400
 
!     SORT COMMAND - IGNORE COMMAND IF AFTER A BASIC DESIGNATOR
 
 170 IF (basic) GO TO 180
 i = 0
 IF (REC(2) == subc) i = 1
 IF (REC(2) == subs) i = 2
 IF (REC(2) == mode) i = 1
 IF (REC(2) == time) i = 1
 IF (REC(2) == freq) i = 1
 IF (i == 0) GO TO 450
 iopt = i
 GO TO 90
 180 WRITE (nout,63660) uwm
 GO TO 90
 
!     BASIC COMMAND - VERIFY SUBSTRUCTURE NAME
 
 190 DO  i = 1,nbs
   i1 = ireq + (i-1)*lbasic + 5
   IF (buf(i1) == REC(2) .AND. buf(i1+1) == REC(3)) GO TO 210
 END DO
 GO TO 220
 210 nss1  = i
 nss2  = i
 basic = .true.
 GO TO 90
 
!     NAME NOT A BASIC - SKIP TO NEXT BASIC, PRINT OR SAVE COMMAND
 
 220 WRITE (nout,63680) uwm,REC(2),REC(3),rss
 225 CALL READ (*9002,*500,casess,REC,3,0,nwds)
 IF (REC(1) == PRINT .OR. REC(1) == SAVE) GO TO 510
 IF (REC(1) == comds(8)) GO TO 190
 GO TO 225
 
!     VELOCITY REQUEST
 
 230 IF (rfno /= 8 .AND. rfno /= 9) GO TO 90
 IF (REC(2) /= NONE) buf(ireq) = 1
 iloc = 9
 GO TO 400
 
!     ACCELERATION REQUEST
 
 240 IF (rfno /= 8 .AND. rfno /= 9) GO TO 90
 IF (REC(2) /= NONE) buf(ireq) = 1
 iloc = 10
 GO TO 400
 
!     ENERGY REQUEST
 
 250 iloc = -1
 GO TO 400
 
!     UIMPROVED REQUEST
 
 260 uimpro = 1
 GO TO 90
 
!     STEPS REQUEST
 
 270 iloc = 11
 GO TO 400
 
!     CHECK VALIDITY OF SET REQUEST
 
 400 IF (REC(2) == -2) GO TO 450
 410 iset = 1
 IF (REC(2) ==  all) iset = -1
 IF (REC(2) == NONE) iset = 0
 IF (iset   <=  0) GO TO 430
 IF (REC(2) == -2) GO TO 420
 IF (REC(2) /= -1) GO TO 450
 
!     INTEGER VALUE
 
 iset = REC(3)
 GO TO 430
 
!     REAL VALUE
 
 420 rset = rrec(3)
 
!     LOOP OVER APPROPRIATE BASIC AREA AND INSERT REQUEST
 
 430 IF (iloc < 0) GO TO 445
 DO  i = nss1,nss2
   i1 = ireq + (i-1)*lbasic + 5 + iloc
   buf(i1) = iset
 END DO
 GO TO 90
 
 445 energy = iset
 GO TO 90
 
!     ILLEGAL COMMAND FORMAT
 
 450 WRITE (nout,63670) uwm,REC(1)
 GO TO 90
 
 
!     END OF RECORD READING CASESS - THIS IS THEREFORE THE LAST
!     SAVE OR PRINT COMMAND
 
 500 loop = -1
 
!     END OF PROCESSING FO THIS PRINT COMMAND
 
 510 CALL CLOSE (casess,rew)
 
!     DETERMINE IF EACH BASIC IS REALLY A BASIC.  IF NOT THEN THESE
!     WILL BE MODAL POINTS
 
!     BASIC   POINT TYPE = 1
!     MODAL   POINT TYPE = 4
 
 maskll = lshift(1023,20)
 DO  i = 1,nbs
   i1 = ireq + (i-1)*lbasic + 5
   buf(i1+12) = 1
   CALL fdsub (buf(i1),idit)
   IF (idit < 0) CYCLE
   CALL fmdi (idit,imdi)
   IF (andf(buf(imdi+ll),maskll) /= 0) buf(i1+12) = 4
 END DO
 CALL sofcls
 RETURN
 
 
!     NO PRINT OR SAVE COMMAND SPECIFIED - GENERATE A SAVE ON
!     THE SOLUTION SUBSTRUCTURE
 
 600 rss(1) = fss(1)
 rss(2) = fss(2)
 loop   = -1
 GO TO 720
 
!     THIS LOOP IS A SAVE COMMAND - SEE IF ANY OTHER COMMANDS FOLLOW
 
 700 CALL READ (*9002,*710,casess,REC,3,0,nwds)
 IF (REC(1) == PRINT .OR. REC(1) == SAVE) GO TO 720
 GO TO 700
 710 loop = -1
 
!     NO OUTPUT BLOCK IS REQUIRED FOR A SAVE COMMAND
 
 720 CALL CLOSE (casess,rew)
 ireq   = 0
 lreq   = 0
 iopt   = 0
 energy = 0
 uimpro = 0
 RETURN
 
!     ERROR RETURNS
 
 800 CALL sofcls
 iopt = -1
 loop = -1
 CALL CLOSE (casess,rew)
 RETURN
 
 6305 WRITE (nout,63050) swm,step,REC(1)
 GO TO 800
 6306 WRITE (nout,63060) uwm,rss
 GO TO 800
 9002 n = -2
 GO TO 9100
 9008 n = -8
 GO TO 9100
 9100 CALL sofcls
 CALL mesage (n,casess,subnam)
 RETURN
 
!     FORMATS
 
 63050 FORMAT (a27,' 6305, RECORD NUMBER',i5,' IS NOT A RECOVER RECORD.',  &
     '  IT IS A ', a4,' RECORD.')
 63060 FORMAT (a25,' 6306, ATTEMPT TO RECOVER DISPLACEMENTS FOR NON-',  &
     'EXISTANT SUBSTRUCTURE ',2A4)
 63660 FORMAT (a25,' 6366, THE RECOVER OUTPUT COMMAND SORT MUST APPEAR ',  &
     'BEFORE THE FIRST BASIC SUBCOMMAND.', /32X,  &
     'ANY OTHER SORT COMMANDS ARE IGNORED.')
 63670 FORMAT (a25,' 6367, ILLEGAL FORMAT ON THE RECOVER OUTPUT COMMAND',  &
     1X,a4,', COMMAND IGNORED.')
 63680 FORMAT (a25,' 6368, THE SUBSTRUCTURE ',2A4,' APPEARING ON A ',  &
     'BASIC COMMAND IS NOT A COMPONENT OF ',2A4, /32X,  &
     'ALL OUTPUT REQUESTS UNTIL THE NEXT BASIC, PRINT OR SAVE ',  &
     'COMMAND ARE IGNORED.')
END SUBROUTINE rcovo
