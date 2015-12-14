SUBROUTINE curv2
!*****
! PASSES NEXT SUBCASE OF ELEMENT STRESS OR STRAIN DATA IN OES1
! AND OUTPUTS OES1M FOR THIS SUBCASE. SETS UP FILES AND TABLES
! FOR -CURV3- IF OES1G IS TO BE FORMED.
 
!     OPEN CORE MAP DURING -CURV2- EXECUTION.
!     =======================================
 
!     FROM-------+------------+
!     CURV1      I  Z(IELTYP) I  MASTER LIST OF ELEMENT TYPES THAT
!     EXECUTION  I    THRU    I  EXIST ON ESTX(SCR1)
!                I  Z(NELTYP) I
!                +------------+
!                I  Z(IMCSID) I  MASTER LIST OF MCSIDS ELEMENTS IN
!                I    THRU    I  PROBLEM REFERENCE, WITH COUNTS OF
!                I  Z(NMCSID) I  OES1M ELEMENTS FOR CURRENT SUBCASE.
!                +------------+
!                I  Z(ICSTM)  I  CSTM FOR EACH -MCSID- IN ABOVE LIST.
!                I    THRU    I  14 WORD ENTRIES. (USER MAY NOT HAVE
!                I  Z(NCSTM)  I  SUPPLIED ALL, BUT MAY BE OK.)
!     FROM-------+------------+
!     AND DURING I  Z(IESTX)  I  SPACE FOR ONE ENTRY OF ESTX(SCR1)
!     CURV2      I    THRU    I  ENTRIES.  (SIZE IS ELEMENT DEPENDENT)
!     EXECUTION  I  Z(NESTX)  I
!                +------------+
!                I  Z(IOES1M) I  TABLE OF INCR-WORD .RIES FOR 1 ELEMENT
!                I    THRU    I  TYPE.  CONTAINS ELEMENT-ID,MCSID,XY-
!                I  Z(NOES1M) I  COMPONENT-CODE, AND (INCR-3) SIGMAS.
!                +------------+  INCR =  9 FOR REAL STRESS
!                I     .      I       = 15 FOR COMPLEX STRESS
!                I     .      I  AVAILABLE CORE.
!                I     .      I
!                I     .      I
!                I     .      I
!                I  Z(JCORE)  I
!                +------------+
!                I  Z(IBUF4)  I  GINO-BUFFER(OES1M)
!                I            I
!                +------------+
!                I  Z(IBUF3)  I  GINO-BUFFER(SCR2 AND SCR3)
!                I            I
!                +------------+
!                I  Z(IBUF2)  I  GINO-BUFFER(SCR1)
!                I            I
!                +------------+
!                I  Z(IBUF1)  I  GINO-BUFFER(OES1)
!                I            I
!                I  Z(LCORE)  I
!                +------------+
 
 
!*****
 REAL :: z(1)     ,rbuf(100),u(9)
 
 INTEGER :: cstms    ,scr1     ,scr2     ,scr3     ,mcb(7)
 INTEGER :: scr4     ,oes1m    ,oes1g    ,oes1     ,scr5
 INTEGER :: cstm     ,est      ,sil      ,gpl
 INTEGER :: eltype   ,subcas   ,FILE     ,estwds
 INTEGER :: ewords   ,owords   ,depts    ,cstype
 INTEGER :: device   ,oldid    ,buf      ,sbuf
 INTEGER :: rd       ,rdrew    ,wrt      ,wrtrew
 INTEGER :: cls      ,clsrew   ,eor      ,sysbuf
 
 LOGICAL :: any      ,eofos1   ,first    ,anyout
 LOGICAL :: foes1g   ,strain   ,any1m    ,any1g
 
 COMMON/BLANK /     ip1      ,ip2      ,icmplx   ,zdum(3)
 
 COMMON/system/     sysbuf   ,ioutpt
 
 COMMON/names /     rd       ,rdrew    ,wrt      ,wrtrew ,clsrew   ,cls
 
 COMMON/condas/     valpi    ,val2pi   ,raddeg   ,degrad ,s4pisq
 
 COMMON/zzzzzz/     iz(1)
 
 COMMON/curvc1/     lsbuf    ,sbuf(10)
 
 COMMON/curvc2/     lbuf     ,buf(100)
 
 COMMON/curvc3/     vec(3)   ,vmax(3)  ,vmin(3)  ,idrec(146)
 
 COMMON/curvtb/     imid     ,nmid     ,lmid     ,nmids  &
     ,ieltyp   ,neltyp   ,jeltyp   ,icstm ,ncstm    ,cstms    ,lcstm    ,iestx  &
     ,nestx    ,imcsid   ,nmcsid   ,lmcsid ,mcsids   ,jmcsid   ,kmcsid   ,isil  &
     ,nsil     ,lsil     ,jsil     ,ioes1m ,noes1m   ,loes1m   ,idep     ,ndep  &
     ,iindep   ,nindep   ,jindep   ,isigma ,nsigma   ,igmat    ,ngmat    ,iext  &
     ,next     ,lext     ,scr1     ,scr2 ,scr3     ,scr4     ,oes1m    ,oes1g  &
     ,oes1     ,mpt      ,cstm     ,est ,sil      ,gpl      ,jcore    ,lcore  &
     ,ibuf1    ,ibuf2    ,ibuf3    ,ibuf4 ,i        ,j        ,k        ,l  &
     ,k1       ,k2       ,ixyz1    ,ixyz2 ,lx1      ,lx2      ,eltype   ,mcsid  &
     ,idscr1   ,idoes1   ,npts     ,npts4 ,iwords   ,nwords   ,subcas   ,kount  &
     ,isig1    ,isig2    ,loc      ,FILE
 COMMON/curvtb/     imsg     ,nelems   ,imatid   ,icomp  &
     ,estwds   ,ewords   ,jp       ,owords  &
     ,matid    ,depts    ,indpts   ,ictype ,ivmat    ,itran    ,cstype   ,ising  &
     ,device   ,oldid    ,any      ,eofos1  &
     ,first    ,anyout   ,foes1g   ,strain ,logerr   ,any1m    ,any1g    ,scr5
 
 EQUIVALENCE        (z(1),iz(1)), (buf(1),rbuf(1))
 EQUIVALENCE        (noeor,rdrew), (eor,cls)
 
 DATA mcb/ 7*1 /
 
!  OPEN OES1M FOR ANY POSSIBLE OUTPUTS DURING THIS SUBCASE PASS.
 
 isig1 = 3
 isig2 = 11
 FILE = oes1m
 loc = 60
 CALL OPEN(*9001,oes1m,iz(ibuf3),wrt)
 
!  OPEN OES1 NOREWIND TO CONTINUE
 
 first = .true.
 any = .false.
 FILE = oes1
 CALL OPEN(*9001,oes1,iz(ibuf1),rd)
 FILE = scr1
 CALL OPEN(*9001,scr1,iz(ibuf2),rdrew)
 
!  ZERO ELEMENT COUNTS FOR EACH -MCSID- THIS SUBCASE MAY REFERENCE.
 
 DO  i = imcsid,nmcsid,2
   iz(i+1) = 0
 END DO
 
!  READ NEXT ID-RECORD
 
 100 FILE = oes1
 loc = 100
 CALL READ(*300,*9003,oes1,idrec(1),146,eor,nwords)
 
!  CHECK IF STILL SAME SUBCASE UNLESS THIS IS THE FIRST ID-RECORD OF A
!  SUBCASE GROUP.
 
 IF( .NOT. first ) GO TO 200
 
!  YES THIS IS FIRST ID-RECORD OF A SUBCASE GROUP.
!  SET SUBCASE IDENTIFIERS.
 
 subcas = idrec(4)
 first = .false.
 GO TO 500
 
!  CHECKING FOR CHANGE IN SUBCASE
 
 200 IF( subcas == idrec(4) ) GO TO 500
 
!  CHANGE IN SUBCASE THUS BACK RECORD OVER THIS ID-RECORD CLOSE
!  OES1, AND WRAP UP OPERATIONS ON OES1M FOR CURRENT SUBCASE.
 
 CALL bckrec( oes1 )
 CALL CLOSE( oes1, cls )
 
!  CLOSE ESTX(SCR1) AND ESTXX(SCR2).
 
 250 CALL CLOSE( scr1, clsrew )
 CALL CLOSE( scr2, clsrew )
 GO TO 5000
 
!  END OF FILE ON OES1. SET EOF FLAG AND WRAP UP CURRENT OPERATIONS
!  ON OES1M.
 
 300 eofos1 = .true.
 CALL CLOSE( oes1, clsrew )
 CALL CLOSE(oes1m,clsrew)
 GO TO 250
 
!  ID RECORD ON OES1 WILL BE FOR SOME KIND OF ELEMENT.
!  CHECK TO SEE IF ITS TYPE IS IN THE LIST OF TYPES NOW ON SCR1
!  WHICH IS THE ABBREVIATED EST. IF NOT THEN SKIP THE DATA RECORD
!  AND GO TO NEXT ID RECORD.
 
 500 eltype = idrec(3)
 iformt = idrec(9)
 owords = idrec(10)
 DO  i = ieltyp,neltyp
   IF( eltype == iz(i) ) GO TO 600
 END DO
 CALL fwdrec(*300,oes1)
 GO TO 100
 
!  POSITION TO SCR1 RECORD FOR THIS ELEMENT TYPE. IF IT CAN NOT BE
!  FOUND BY FORWARD SEARCH THERE IS A LOGIC ERROR, OR OES1 ELEMENT
!  TYPES ARE NOT IN SAME ORDER AS EST ELEMENT TYPES.
 
 600 FILE = scr1
 loc = 600
 CALL REWIND (scr1)
 640 CALL READ(*9002,*9003,scr1,buf(1),3,noeor,nwords)
 IF( buf(1) == eltype ) GO TO 650
 CALL fwdrec(*9002,scr1)
 GO TO 640
 
!  NOW POSITIONED TO READ ELEMENT ENTRIES FROM ESTX(SCR1) WHICH
!  ARE OK FOR INCLUSION IN OES1M AND OES1G PROCESSING.
 
!  ALSO POSITIONED TO READ OUTPUT STRESS/STRAIN ENTRIES FROM OES1.
!  HOWEVER, ONLY THOSE ALSO ON ESTX(SCR1) WILL BE PULLED.
 
 650 anyout = .false.
 ewords = buf(2)
 npts = buf(3)
 npts4 = 4*npts
 iestx = ncstm + 1
 nestx = ncstm + ewords
 ioes1m = nestx + 1
 noes1m = nestx
 loc = 650
 icrq = noes1m - jcore
 IF( noes1m > jcore ) GO TO 9008
 idscr1 = 0
 
!  READ NEXT OES1 ENTRY AND SET IDOES1  (STRIPPING OFF DEVICE CODE)
 
 670 FILE = oes1
 loc = 670
 CALL READ(*9002,*900,oes1,buf(1),owords,noeor,nwords)
 idoes1 = buf(1) / 10
 IF( idoes1  > 0) THEN
   GO TO   700
 ELSE
   GO TO  9000
 END IF
 
!  READ NEXT SCR1 ENTRY AND SET IDSCR1
 
 680 FILE = scr1
 loc = 680
 CALL READ(*9002,*950,scr1,iz(iestx),ewords,noeor,nwords)
 idscr1 = iz(iestx)
 IF( idscr1  > 0) THEN
   GO TO   700
 ELSE
   GO TO  9000
 END IF
 
!  CHECK FOR MATCH OF ESTX(SCR1) ENTRY ID WITH OES1 ENTRY ID.
 
 700 IF( idoes1 - idscr1  < 0) THEN
   GO TO   670
 ELSE IF ( idoes1 - idscr1  == 0) THEN
   GO TO   710
 ELSE
   GO TO   680
 END IF
 
!  MATCH FOUND THUS BEGIN OES1M ENTRY CALCULATIONS
 
 710 mcsid = iz(iestx+1)
 loc = 710
 CALL tranem( mcsid, npts, z(iestx+npts+2), icomp, u(1), vec(1) )
 
!  FORM AND ADD ENTRY TO CORE. INVARIANTS WILL BE COMPUTED LATER.
 
 incr = 9
 IF (icmplx == 1) incr = 15
 icrq = noes1m + incr - jcore
 IF( noes1m+incr > jcore ) GO TO 9008
 iz(noes1m+1) = buf(1)
 iz(noes1m+2) = mcsid
 iz(noes1m+3) = icomp
 IF (icmplx == 1) GO TO 730
 
!  IF STRAINS DO MODIFICATION OF GAMMA
 
 IF( .NOT. strain ) GO TO 720
 rbuf(isig1+2) = rbuf(isig1+2) / 2.0
 rbuf(isig2+2) = rbuf(isig2+2) / 2.0
 720 CALL gmmats( u(1),3,3,0,  rbuf(isig1),3,1,0,  z(noes1m+4) )
 CALL gmmats( u(1),3,3,0,  rbuf(isig2),3,1,0,  z(noes1m+7) )
 
 noes1m = noes1m + 9
 GO TO 740
 
 730 IF (iformt /= 3) GO TO 732
 DO  mm1 = 3, 10, 7
   mm2 = mm1 + 4
   DO  lll = mm1, mm2, 2
     ztemp   = rbuf(lll)*COS(rbuf(lll+1)*degrad)
     rbuf(lll+1) = rbuf(lll)*SIN(rbuf(lll+1)*degrad)
     rbuf(lll) = ztemp
   END DO
 END DO
 732 zdum(1) = rbuf(3)
 zdum(2) = rbuf(5)
 zdum(3) = rbuf(7)
 CALL gmmats (u(1),3,3,0,zdum,3,1,0,z(noes1m+4))
 zdum(1) = rbuf(4)
 zdum(2) = rbuf(6)
 zdum(3) = rbuf(8)
 CALL gmmats (u(1),3,3,0,zdum,3,1,0,z(noes1m+7))
 zdum(1) = rbuf(10)
 zdum(2) = rbuf(12)
 zdum(3) = rbuf(14)
 CALL gmmats (u(1),3,3,0,zdum,3,1,0,z(noes1m+10))
 zdum(1) = rbuf(11)
 zdum(2) = rbuf(13)
 zdum(3) = rbuf(15)
 CALL gmmats (u(1),3,3,0,zdum,3,1,0,z(noes1m+13))
 
 IF (iformt /= 3) GO TO 738
 DO  mm1 = 4, 10, 6
   mm2 = mm1 + 2
   DO  lll = mm1, mm2
     ll1 = noes1m + lll
     ll2 = ll1 + 3
     ztemp   = SQRT (z(ll1)**2 + z(ll2)**2)
     IF (ztemp /= 0.0) GO TO 734
     z(ll2) = 0.0
     GO TO 736
     734 z(ll2) = ATAN2 (z(ll2), z(ll1))*raddeg
     IF (z(ll2) < -0.00005E0) z(ll2) = z(ll2) + 360.0
     736 z(ll1) = ztemp
   END DO
 END DO
 
 738 noes1m = noes1m + 15
 
 
!  IF THIS IS THE FIRST ELEMENT ENTRY TO BE FOUND
!  AND OES1G IS TO BE FORMED, THE ID-RECORD IS SAVED FOR USE BY
!  CURV3 OVERLAY.
 
 740 IF( .NOT. foes1g ) GO TO 790
 IF( any ) GO TO 750
 FILE = scr3
 loc = 740
 CALL OPEN(*9001,scr3,iz(ibuf4),wrtrew)
 CALL WRITE( scr3, idrec(1), 146, eor )
 CALL CLOSE( scr3, clsrew )
 
 FILE = scr2
 CALL OPEN(*9001,scr2,iz(ibuf4),wrtrew)
 device = MOD( buf(1), 10 )
 any = .true.
 
!  OUTPUT SPECIAL ESTXX (SCR2) ENTRY FOR USE BY CURV3.
 
 750 CALL WRITE( scr2, mcsid, 1, noeor )
 CALL WRITE( scr2, z(noes1m-5), 6, noeor )
 CALL WRITE( scr2, vec(1), 3, noeor )
 CALL WRITE( scr2, npts, 1, noeor )
 CALL WRITE( scr2, iz(iestx+2), npts4, noeor )
 790 CONTINUE
 GO TO 670
!*****
!  END OF ENTRY DATA POSSIBLE FOR THIS ELEMENT TYPE
!******
 
!  SKIP ANY UNUSED DATA IN ESTX (SCR1) DATA RECORD FOR THIS ELEMENT TYPE
 
 900 FILE = scr1
 loc = 900
 CALL fwdrec(*9002,scr1)
 GO TO 960
 
!  SKIP ANY UNUSED DATA IN OES1 DATA RECORD FOR THIS ELEMENT TYPE.
 
 950 FILE = oes1
 loc = 950
 CALL fwdrec(*9002,oes1)
 
!  IF ANY ENTRIES WERE FOUND AND COMPLETED AND PLACED IN CORE
!  THEY ARE SORTED ON -MCSID- AND OUTPUT. AS THEY ARE OUTPUT
!  THE INVARIANTS ARE COMPUTED.
 
 960 IF( noes1m < ioes1m ) GO TO 100
 
!  YES THERE ARE SOME ENTRIES
 
 loes1m = noes1m - ioes1m + 1
 CALL sort( 0, 0, incr, 2, iz(ioes1m), loes1m )
 
!  OUTPUT ID-RECORD, REDEFINE MAJOR-ID FOR OFP MODULE
 
!  RE-DEFINITION MISSING FOR NOW.
 
 idrec(3) = idrec(3) + 1000
 CALL WRITE( oes1m, idrec(1), 146, eor )
 mcb(1) = oes1m
 CALL wrttrl( mcb(1) )
 any1m = .true.
 
!  MOVE AXIS CODE AND COMPLETE INVARIANTS OF EACH ENTRY.
 
 kmcsid = iz(ioes1m+1)
 kount = 0
 
 DO  i = ioes1m, noes1m, incr
   buf(1) = iz(i)
   buf(2) = iz(i+1)
   rbuf(3) = z(i+3)
   IF (icmplx == 1) GO TO 963
   rbuf(4) = z(i+4)
   rbuf(5) = z(i+5)
   CALL curvps( rbuf(3), rbuf(6) )
   IF( .NOT. strain ) GO TO 961
   rbuf(5) = 2.0 * rbuf(5)
   rbuf(9) = 2.0 * rbuf(9)
   961 buf(10) = iz(i+2)
   rbuf(11) = z(i+6)
   rbuf(12) = z(i+7)
   rbuf(13) = z(i+8)
   CALL curvps ( rbuf(11), rbuf(14) )
   IF( .NOT. strain ) GO TO 962
   rbuf(13) = 2.0 * rbuf(13)
   rbuf(17) = 2.0 * rbuf(17)
   962 CALL WRITE( oes1m, buf(1), 17, noeor )
   GO TO 964
   963 rbuf( 4) = z(i+ 6)
   rbuf( 5) = z(i+ 4)
   rbuf( 6) = z(i+ 7)
   rbuf( 7) = z(i+ 5)
   rbuf( 8) = z(i+ 8)
   buf(  9) = iz(i+2)
   rbuf(10) = z(i+ 9)
   rbuf(11) = z(i+12)
   rbuf(12) = z(i+10)
   rbuf(13) = z(i+13)
   rbuf(14) = z(i+11)
   rbuf(15) = z(i+14)
   CALL WRITE (oes1m, buf(1), 15, noeor)
   
!  KEEP COUNT OF ELEMENTS IN EACH MCSID GROUP
   
   964 IF( iz(i+1) /= kmcsid ) GO TO 965
   kount = kount + 1
   IF( i+incr-1 < noes1m ) CYCLE
   
!  CHANGE IN -MCSID- OF OUTPUT ENTRIES OR LAST ENTRY.
!  ADD COUNT OF ELEMENTS OF CURRENT TYPE TO TOTAL COUNT
!  OF ELEMENTS OF THIS -MCSID-.
   
   965 loc = 965
   CALL bisloc(*9000,kmcsid,iz(imcsid),2,mcsids,jp)
   iz(imcsid+jp) = iz(imcsid+jp) + kount
   kount = 1
   kmcsid = iz(i+1)
 END DO
 CALL WRITE( oes1m, 0, 0, eor )
 GO TO 100
!*****
!  ALL PROCESSING OF ONE SUBCASE COMPLETE FOR OES1M.
!*****
 5000 CALL CLOSE( oes1m, cls )
 RETURN
!*****
!  ERROR CONDITION ENCOUNTERED
!*****
 9000 imsg = -logerr
 GO TO 5000
 9001 imsg = -1
 GO TO 5000
 9002 imsg = -2
 GO TO 5000
 9003 imsg = -3
 GO TO 5000
 9008 imsg = -8
 lcore = icrq
 GO TO 5000
END SUBROUTINE curv2
