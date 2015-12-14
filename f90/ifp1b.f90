SUBROUTINE ifp1b
     
!     THIS ROUTINE DETERMINES THE LOOP CONDITIONS AND CASE CONTROL
!     REQUEST CHNGES
!     LOOP$ -- THE CURRENT PROBLEM WILL LOOP
 
!     LOOP1$-- THE OLD PROBLEM WAS A LOOP AND CASE CONTROL IS CHANGED
!     IN LENGTH
 
!     COMMENTS FROM G.C.  10/92
!     IWORD AND IBIT 200 WORDS EACH CORRESPOND TO 200 WORDS IN CASECC
!     ZERO IN IWORD MEANS NO FURTHER CHECKING
!     INTEGER VALUE IN IBIT POINTS TO RESTART BIT POSITION, AND WILL BE
!     SAVED IN BITS(17) AND BITS(18), BITS FOR LCC. (LBD = 16)
 
!     LAST REVISED  7/91, BY G.CHAN/UNISYS, TO ALLOW HUGE THRU-RANGE ON
!     SET IN CASE CONTROL SECTION FOR PRINTOUT OR PLOTTING
 
 EXTERNAL        andf,orf
 LOGICAL :: NEW,debug
 INTEGER :: NAME(2),optp,casecc,bits,two1,core(2),case,cc,ss,  &
     orf,andf,corey(401)
 DIMENSION       icase(200,2),iword(200),ibit(200)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /two   / two1(32)
 COMMON /system/ ibuf,nout
 COMMON /ifpx0 / lbd,lcc,bits(1)
 COMMON /xifp1 / iblank
 COMMON /zzzzzz/ corex(1)
 COMMON /ifp1a / scr1,casecc,is,nwpc,ncpw4,nmodes,icc,nset,nsym,  &
     zzzzbb,istr,isub,lencc,iben,equal,IEOR
 EQUIVALENCE    (corex(1),corey(1),icase(1,1)),(core(1),corey(401))
 DATA    NAME  / 4HIFP1,  4HB          /
 DATA    case  , cc     / 4HCASE,4HCC  /
 DATA    ss    / 4HSS   /
 DATA    optp  / 4HOPTP /
 DATA    iword /  &
     -1,01,01,01,01,01,01,01,01,00,01,01,01,01,01,-1,00,01,01,00,  &
     01,01,00,01,01,00,01,01,00,01,01,00,01,01,00,01,01,01,-1,-1,  &
     -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,  &
     -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,  &
     -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,  &
     -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,  &
     -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,00,-1,-1,01,01,01,  &
     01,01,01,01,00,00,-1,01,01,-1,00,01,01,00,01,01,00,01,01,01,  &
     -1,-1,01,01,01,-1,00,01,01,00,01,01,00,01,01,00,01,01,01,00,  &
     01,01,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1/
 DATA    ibit  /  &
     00,02,03,04,05,06,07,08,09,10,10,10,13,14,15,00,18,18,18,18,  &
     18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,17,00,00,  &
     00,00,00,00,00,00,00,00,00,00,00,00,00,00,00,00,00,00,00,00,  &
     00,00,00,00,00,00,00,00,00,00,00,00,00,00,00,00,00,00,00,00,  &
     00,00,00,00,00,00,00,00,00,00,00,00,00,00,00,00,00,00,00,00,  &
     00,00,00,00,00,00,00,00,00,00,00,00,00,00,00,00,00,00,00,00,  &
     00,00,00,00,00,00,00,00,00,00,00,00,00,00,16,00,00,20,21,21,  &
     22,22,23,23,18,18,00,24,25,00,10,10,10,10,10,10,10,10,10,27,  &
     00,00,30,26,29,00,18,18,18,18,18,18,10,10,10,10,10,10,33,18,  &
     18,18,00,00,00,00,00,00,00,00,00,00,00,00,00,00,00,00,00,00/
 DATA    NEW,debug   / 2*.false. /
 
 k  = lbd + 1
 ifirod = 0
 ioloop = 0
 iloop  = 0
 ifirst = 0
 ieoptp = 0
 
!     ALLOCATE GINO BUFFERS
 
 nz     = korsz(core)
 ibuf1  = nz - ibuf + 1
 ibuf2  = ibuf1 - ibuf
 nz     = nz - 2*ibuf
 icrq   =-nz
 IF (nz <= 0) GO TO 700
 iecase = 0
 
!     TRY TO FIND CASECC ON OPTP - TRY TO ASSUME PROPER POSITION
 
 CALL OPEN (*560,optp,core(ibuf1),2)
 iopn = 0
 
!     FIND CASECC
 
 10 CALL READ (*600,*610,optp,core(1),2,1,iflag)
 IF (core(1) == case .AND. core(2) == cc) GO TO 30
 IF (core(1) == case .AND. core(2) == ss) GO TO 20
 IF (iopn == 0) CALL REWIND (optp)
 iopn = 1
 CALL skpfil (optp,1)
 GO TO 10
 
!     CASESS FOUND ON OPTP - SKIP TO CASECC
 
 20 CALL READ (*600,*610,optp,core(1),2,1,iflag)
 IF (core(1) /= case .OR. core(2) /= cc) GO TO 20
 
!     CASECC FOUND ON OLD PROB TAPE
 
!     OPEN CASECC AND SKIP CASESS IF PRESENT
 
 30 CALL OPEN (*650,casecc,core(ibuf2),0)
 40 CALL READ (*670,*680,casecc,core(1),2,1,iflag)
 IF (core(1) /= case .OR. core(2) /= cc) GO TO 40
 ASSIGN 50 TO ihop
 50 CALL READ (*550,*610,optp,icase(1,2),lencc,0,iflag)
 IF (icase(16,2) == 0) GO TO 60
 CALL fwdrec (*600,optp)
 GO TO 50
 60 IF (icase(lencc,2) == 0) GO TO 70
 lsym = icase(lencc,2)
 CALL READ (*600,*610,optp,core(1),-lsym,0,iflag )
 70 CALL READ (*600,*80 ,optp,core(1),   nz,1,ifoptp)
 icrq = nz
 GO TO 700
 80 CALL READ (*510,*680,casecc,icase(1,1),lencc,0,iflag)
 IF (icase(16,1) == 0) GO TO 90
 CALL fwdrec (*670,casecc)
 GO TO 80
 90 IF (icase(lencc,1) == 0) GO TO 100
 lsym = icase(lencc,1)
 CALL READ (*670,*680,casecc,core(ifoptp+1),-lsym,0,iflag)
 100 CALL READ (*670,*110,casecc,core(ifoptp+1),nz-ifoptp,1,ifcase)
 icrq = nz - ifoptp
 GO TO 700
 
!     CHECK FOR LOOPING PROBLEM
 
 110 IF (ifirst /= 0) GO TO 120
 ifirst= 1
 ispc  = icase(  3,1)
 impc  = icase(  2,1)
 imtd  = icase(  5,1)
 ifreq = icase( 14,1)
 itfl  = icase( 15,1)
 ik1   = icase(139,1)
 ik2   = icase(140,1)
 im1   = icase(141,1)
 im2   = icase(142,1)
 ib1   = icase(143,1)
 ib2   = icase(144,1)
 IF (icase(165,1) > 0 .OR. icase(164,1) > 0) GO TO 130
 GO TO 140
 120 IF (icase(  3,1) /= ispc) GO TO 130
 IF (icase(  2,1) /= impc) GO TO 130
 IF (icase(  5,1) /= imtd) GO TO 130
 IF (icase(139,1) /= ik1 .OR. icase(140,1) /= ik2) GO TO 130
 IF (icase(141,1) /= im1 .OR. icase(142,1) /= im2) GO TO 130
 IF (icase(143,1) /= ib1 .OR. icase(144,1) /= ib2) GO TO 130
 IF (icase( 15,1) /= itfl ) GO TO 130
 IF (icase( 14,1) /= ifreq) GO TO 130
 IF (icase(138,1) > 0) GO TO 130
 IF (icase( 38,1) /= 0) GO TO 130
 GO TO 140
 
!     SET LOOP$
 
 130 bits(k) = orf(bits(k),two1(11))
 iloop = 1
 140 CONTINUE
 
!     DETERMINE IF OLD PROBLEM WOULD HAVE LOOPED
 
 IF (ifirod /= 0) GO TO 150
 ifirod= 1
 ispc1 = icase(  3,2)
 impc1 = icase(  2,2)
 imtd1 = icase(  5,2)
 ik11  = icase(139,2)
 ik21  = icase(140,2)
 im11  = icase(141,2)
 im21  = icase(142,2)
 ib11  = icase(143,2)
 ib21  = icase(144,2)
 itfl1 = icase( 15,2)
 ifreq1= icase( 14,2)
 IF (icase(164,2) > 0 .OR. icase(165,2) > 0) GO TO 160
 GO TO 170
 
!     SECOND RECORD APPLY LOOP RULES
 
 150 IF (icase(  3,2) /= ispc1) GO TO 160
 IF (icase(  2,2) /= impc1) GO TO 160
 IF (icase(  5,2) /= imtd1) GO TO 160
 IF (icase(139,2) /= ik11 .OR. icase(140,2) /= ik21) GO TO 160
 IF (icase(141,2) /= im11 .OR. icase(142,2) /= im21) GO TO 160
 IF (icase(143,2) /= ib11 .OR. icase(144,2) /= ib21) GO TO 160
 IF (icase(138,2) > 0) GO TO 160
 IF (icase( 38,2) /= 0) GO TO 160
 IF (icase( 15,2) /= itfl1 ) GO TO 160
 IF (icase( 14,2) /= ifreq1) GO TO 160
 GO TO 170
 160 ioloop = 1
 170 CONTINUE
 IF (iecase /= 1) GO TO 180
 IF (ioloop == 1) GO TO 530
 GO TO 520
 
!     CHECK FOR CHANGES -
 
 180 IF (ieoptp == 1) ieoptp = 2
 DO  i = 1,lencc
   IF (ibit(i) == 0) CYCLE
   l = ibit(i)
   IF (l <= 32 .AND. andf(bits(k),two1(l)) /= 0) CYCLE
   IF (iword(i) == 0) GO TO 210
   IF (icase(i,1) == icase(i,2)) CYCLE
   190 IF (l > 32) GO TO 200
   bits(k) = orf(bits(k),two1(l))
   CYCLE
   
!     SECOND CASECC WORD
   
   200 l = l - 31
   bits(k+1) = orf(bits(k+1),two1(l))
   CYCLE
   
!     CHECK FOR PRESENCE OF PRINT AND PLOT REQUESTS
   
   210 IF (i /= 135 .AND. icase(i,1) == 0) CYCLE
   IF (i /= 135) GO TO 220
   IF (icase(i,1) /= 0 .OR. icase(i,2) /= 0) GO TO 190
   CYCLE
   220 IF (ibit(i) == 18) bits(k+1) = orf(bits(k+1),two1(3))
   IF (ibit(i) == 10) bits(k+1) = orf(bits(k+1),two1(4))
   IF (ieoptp  ==  2) GO TO 190
   IF (icase(i,1) < 0 .AND. icase(i,2) < 0) CYCLE
   IF (icase(i,1) < 0 .AND. icase(i,2) >= 0) GO TO 190
   IF (icase(i,1) > 0 .AND. icase(i,2) <= 0) GO TO 190
   ipcase = ifoptp + 1
   230 IF (ipcase > ifoptp+ifcase) GO TO 690
   IF (core(ipcase) == icase(i,1)) GO TO 240
   ipcase = ipcase + core(ipcase+1) + 2
   GO TO 230
   240 ipoptp = 1
   250 IF (ipoptp > ifoptp) GO TO 620
   IF (core(ipoptp) == icase(i,2)) GO TO 260
   ipoptp = ipoptp + core(ipoptp+1) + 2
   GO TO 250
   260 iqcase = ifoptp + ifcase + 1
   ix = ipcase
   iy = iqcase
   ASSIGN 280 TO jump
   IF (debug) WRITE (nout,270)
   270 FORMAT (/,' ------ NPTP PASS ------')
   GO TO 360
   280 iqoptp = iy
   ix = ipoptp
   ASSIGN 300 TO jump
   IF (debug) WRITE (nout,290)
   290 FORMAT (/,' ------ OPTP PASS ------')
   GO TO 360
   300 leng1 = iqoptp - iqcase
   leng2 = iy - iqoptp
   IF (debug) WRITE (nout,310) core(ipcase),leng1,leng2,iy,iqoptp, iqcase
   310 FORMAT (//,' IFP1B/@310  CHECKING SETS',i9,' FROM NPTP AND OPTP',  &
       /5X,'LENG1,LENG2, IY,IQOPTP,IQCASE =', 2I5,3I7)
   IF (leng1 /= leng2) GO TO 340
   DO  mm = 1,leng1
     IF (core(iqcase+mm-1) /= core(iqoptp+mm-1)) GO TO 340
   END DO
   IF (debug) WRITE (nout,330) core(ipcase)
   330 FORMAT (' ... NO DIFFERENCES IN SET',i8)
   CYCLE
   340 WRITE  (nout,350) uim,core(ipcase)
   350 FORMAT (a29,', SET',i9,' DEFINITION HAS BEEN CHANGED IN RESTART')
   GO TO 190
   
!     A NEW NON-EXPANDING METHOD IS IMPLEMENTED HERE BY  G.CAHN/UNISYS
!     8/91, IN CASE THE ORIGINAL LOGIC RUNS OUT OF CORE SPACE
   
!     THE NEW METHOD WILL CONCATINATE VARIATIONS OF SET DEFINITION TO
!     THE SIMPLEST FORM.  E.G. THE NEXT 3 LINES SPECIFY THE SAME SET
!     10 THRU 400000  (THIS IS THE SIMPLEST FORM)
!     10, 11, 12 THRU 400008, 400009, 400000
!     10 THRU 20, 21, 22, 23 THRU 200, 201 THRU 500, 501 502 THRU 400000
   
   360 IF (NEW) GO TO 420
   in = core(ix+1)
   ix = ix + 2
   m  = 0
   370 m  = m + 1
   IF (m-in < 0) THEN
     GO TO   380
   ELSE IF (m-in == 0) THEN
     GO TO   400
   ELSE
     GO TO   490
   END IF
   380 IF (core(ix+m) > 0) GO TO 400
   m1 = core(ix+m-1)
   m2 =-core(ix+m  )
   icrq = iy + m2 - m1 - nz
   IF (icrq > 0) GO TO 410
   DO  mm = m1,m2
     core(iy) = mm
     iy = iy + 1
   END DO
   m = m + 1
   GO TO 370
   400 icrq = iy - nz
   IF (iy > nz) GO TO 700
   core(iy) = core(ix+m-1)
   iy = iy + 1
   GO TO 370
   
!     INSUFFICIENT CORE SPACE, SWITCH TO NEW METHOD
   
   410 NEW = .true.
   GO TO 260
   
!     NEW LOGIC WITHOUT THRU RANGE EXPANSION
   
   420 in = core(ix+1)
   ix = ix + 2
   m0 = iy
   core(iy) = core(ix)
   iy = iy + 1
   IF (in == 1) GO TO 490
   core(iy) = core(ix+1)
   IF (core(iy) == core(iy-1)+1) core(iy) = -core(iy)
   iy = iy + 1
   m  = 1
   430 m  = m + 1
   IF (m >= in) GO TO 470
   m1 = core(ix+m)
   m2 = IABS(m1)
   IF (debug) WRITE (nout,440) m,in,ix,iy,m1,core(iy-1)
   440 FORMAT (' @440   M,IN,IX,IY,M1,CORE(IY-1) =',6I8)
   IF (m1 < 0) GO TO 450
   IF (m1 /= 1-core(iy-1)) GO TO 460
   core(iy-1) = -m2
   GO TO 430
   450 IF (core(iy-1) > 0) GO TO 460
   core(iy-1) = -m2
   GO TO 430
   460 core(iy) = m1
   IF (m1 == core(iy-1)+1) core(iy) = -m2
   iy = iy + 1
   GO TO 430
   470 icrq = iy - nz
   IF (iy > nz) GO TO 700
   m1 = iy - 1
   IF (debug) WRITE (nout,480) core(ix-2),(core(j),j=m0,m1)
   480 FORMAT (/,' IFP1B/@480    SET',i8, /,(2X,15I8))
   
   490 GO TO jump, (280,300)
   
 END DO
 GO TO ihop, (50,80)
 
!     EOF ON CASECC
 
 510 CALL CLOSE (casecc,1)
 IF (ieoptp /= 0) GO TO 530
 iecase = 1
 GO TO 150
 520 CALL READ (*530,*610,optp,icase(1,2),lencc,1,iflag)
 IF (icase(16,2) /= 0) GO TO 520
 GO TO 150
 530 CALL CLOSE (optp,2)
 IF (ieoptp == 1 .OR. ioloop == 0) GO TO 540
 
!     SET LOOP1  THIS SHOULD REEXECUTE THE ENTIRE LOOP
 
 bits(k) = orf(bits(k),two1(12))
 
!     CHECK FOR LOOP$ IF NOT ON SET NOLOOP$
 
 540 IF (iloop == 0) bits(k) = orf(bits(k),two1(32))
 RETURN
 
!     EOF ON  OPTP
 
 550 ASSIGN 80 TO ihop
 ieoptp = 1
 GO TO 80
 
!     ERROR MESSAGES
 
 560 ip1 = -1
 570 ip2 = optp
 580 CALL mesage (ip1,ip2,NAME)
 RETURN
 
 600 ip1 = -2
 GO TO 570
 610 ip1 = -3
 GO TO 570
 620 core(1) = optp
 core(2) = iblank
 630 WRITE  (nout,640) sfm,core(1),core(2)
 640 FORMAT (a25,' 651, LOGIC ERROR IN SUBROUTINE IFP1B WHILE ',  &
     'PROCESSING SET DATA ON ',2A4,' FILE.')
 ip1 = -37
 GO TO 580
 650 ip1 = -1
 660 ip2 = casecc
 GO TO 580
 670 ip1 = -2
 GO TO 660
 680 ip1 = -3
 GO TO 660
 690 core(1) = case
 core(2) = cc
 GO TO 630
 700 ip1 = -8
 ip2 = icrq
 GO TO 580
END SUBROUTINE ifp1b
