SUBROUTINE curv3
!*****
!  THIS OVERLAY WILL FORM OES1G (IF REQUESTED BY DMAP PARAMETER = 0)
 
!  OES1G OUTPUTS FOR CURRENT SUBCASE WILL BE GROUPED ON THE BASIS OF
!  THE MCSID. THUS THERE WILL BE A PASS FOR EACH -MCSID- HAVING A NON-
!  ZERO COUNT IN TABLE(IMCSID-NMCSID).
 
!  TO CONSERVE CORE FOR SSPLIN UTILITY, THE SIGMAS FOR EACH -MCSID- PASS
!  WILL BE WRITTEN TO SCR4 AS ENTRIES ARE SELECTED FROM SCR2. DEPENDENT
!  POINTS, EXTERNAL-IDS, AND INDEPENDENT POINTS WILL BE PLACED IN CORE
!  AND THEN REDUCED DURING THE PROJECTION SURFACE DETERMINATION PHASE.
 
!     OPEN CORE MAP DURING -CURV3- EXECUTION.
!     =======================================
 
!     FROM-------+------------+
!     CURV1      I  Z(IELTYP) I  MASTER LIST OF ELEMENT TYPES  ON
!     EXECUTION  I    THRU    I  ESTX(SCR1)
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
!     AND DURING I  Z(IINDEP) I  INDEPENDENT POINT COORDINATES FOR ONE
!     CURV3      I    THRU    I  -MCSID- OF CURRENT SUBCASE.
!     EXECUTION  I  Z(NINDEP) I  TWO OR THREE WORD ENTRIES POSSIBLE.
!                +------------+
!                I  Z(IDEP)   I  DEPENDENT POINT COORDINATES FOR ONE
!                I    THRU    I  -MCSID- OF CURRENT SUBCASE.
!                I  Z(NDEP)   I  TWO OR FOUR WORD ENTRIES POSSIBLE.
!                +------------+
!                I  Z(IGMAT)  I  G MATRIX FROM SSPLIN UTILITY
!                I    THRU    I  (N-DEPENDENT-PTS BY N-INDEPENDENT-PTS)
!                I  Z(NGMAT)  I
!                +------------+
!                I  Z(ISIGMA) I  OES1M SIGMAS FOR ONE -MCSID- OF CURRENT
!                I    THRU    I  SUBCASE.  6X1 ENTRIES.
!                I  Z(NSIGMA) I
!                +------------+
!                I     .      I  AVAILABLE CORE.
!                I     .      I  (SSPLIN UTILITY USES Z(ISIGMA) THRU
!                I     .      I  Z(LCORE) FOR WORKING SPACE.)
!                I     .      I
!                I  Z(JCORE)  I
!                +------------+
!                I  Z(IBUF4)  I  GINO-BUFFER
!                I            I
!                +------------+
!                I  Z(IBUF3)  I  GINO-BUFFER
!                I            I
!                +------------+
!                I  Z(IBUF2)  I  GINO-BUFFER
!                I            I
!                +------------+
!                I  Z(IBUF1)  I  GINO-BUFFER
!                I            I
!                I  Z(LCORE)  I
!                +------------+
 
!  INPUTS - SCR2 CONTAINING ACTUAL ELEMENT ENTRIES USED TO FORM
!                OES1M FOR CURRENT OES1 SUBCASE. MAY BE MORE THAN
!                ONE -MCSID-. HAS THE SIX SIGMAS OF EACH ELEMENT
!                APPENDED TO EACH ELEMENT.
 
!                       -ELEMENT-ENTRY-
 
!                        MCSID = MATERIAL COORDINATE SYSTEM ID
!                        SIGMA1-X
!                        SIGMA1-Y
!                        SIGMA1-XY
!                        SIGMA2-X
!                        SIGMA2-Y
!                        SIGMA2-XY
!                        XC  *
!                        YC   * MEAN CENTER OF INDEPENDENT POINT
!                        ZC  *
!                        NPTS = NUMBER OF CONNECTED DEPENDENT GRIDS
!                        EXTERNAL GRID IDS (1 FOR EACH POINT)
!                        X,Y,Z  COMPONENTS OF EACH DEPENDENT GRID
 
!           SCR3 CONTAINING OFP TYPE -ID- RECORD TO USE AS A MODEL
!                FOR OES1G -ID- RECORD.
 
 
!           TABLE(IMCSID) THRU Z(NMCSID) CONTAINS PAIRS OF MCSID-S AND
!                COUNTS. (ONE PAIR FOR EACH UNIQUE MCSID OF CURRENT
!                SUBCASE.)
 
 
!           TABLE  Z(ICSTM) TO Z(NCSTM) CONTAINING TRANSFORMATIONS
 
!*****
 REAL :: z(1)     ,rbuf(100)
 
 INTEGER :: mcb(7)
 
 INTEGER :: cstms    ,scr1     ,scr2     ,scr3
 INTEGER :: scr4     ,oes1m    ,oes1g    ,oes1     ,scr5
 INTEGER :: cstm     ,est      ,sil      ,gpl
 INTEGER :: eltype   ,subcas   ,FILE     ,estwds
 INTEGER :: ewords   ,owords   ,depts    ,cstype
 INTEGER :: device   ,oldid    ,buf      ,sbuf
 INTEGER :: rd       ,rdrew    ,wrt      ,wrtrew
 INTEGER :: cls      ,clsrew   ,eor      ,sysbuf
 
 LOGICAL :: any      ,eofos1   ,first    ,anyout
 LOGICAL :: foes1g   ,strain   ,any1m    ,any1g
 
 COMMON/BLANK /     ip1      ,ip2
 
 COMMON/system/     sysbuf   ,ioutpt
 
 COMMON/names /     rd       ,rdrew    ,wrt      ,wrtrew ,clsrew   ,cls
 
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
 
!  BRING OES1G -ID- RECORD INTO CORE AND MODIFY AS NECESSARY.
 
 FILE = scr3
 loc = 50
 CALL OPEN(*9001,scr3,iz(ibuf1),rdrew)
 CALL READ(*9002,*9003,scr3,idrec(1),146,noeor,nwords)
 CALL CLOSE( scr3, clsrew )
 
 
 
 
!  OVERALL LOOP IS ON ENTRIES OF TABLE(IMCSID-NMCSID)
 
 jmcsid = imcsid
 
 100 mcsid = iz(jmcsid)
 indpts = iz(jmcsid+1)
 loc = 100
 IF( indpts  < 0) THEN
   GO TO  9000
 ELSE IF ( indpts  == 0) THEN
   GO TO   980
 END IF
 
!  COLLECT DATA REQUIRED FROM SCR2.
 
 110 FILE = scr2
 
!  CORE ALLOCATION FOR XC, YC, ZC OF EACH INDEPENDENT POINT.
 
 iindep = ncstm + 1
 nindep = ncstm + 3*indpts
 
!  CORE ALLOCATION FOR EXT-ID,X,Y,Z OF EACH UNIQUE DEPENDENT POINT.
!  (THE QUANTITY OF DEPENDENT POINTS IS NOT YET KNOWN.)
 
 idep = nindep + 1
 ndep = nindep
 loc = 110
 icrq = ndep - jcore
 IF( ndep > jcore ) GO TO 9008
 
 CALL OPEN(*9001,scr2,iz(ibuf1),rdrew)
 FILE = scr3
 CALL OPEN(*9001,scr3,iz(ibuf2),wrtrew)
 
 jindep = iindep
 FILE = scr2
 
!  FIND -INDPTS- NUMBER OF INDEPENDENT ELEMENT POINTS ENTRIES
!  FOR CURRENT -MCSID- PASS. (LOGIC ERROR IF CAN NOT FIND THIS MANY)
 
 DO  i = 1,indpts
   
!  READ ELEMENT INDEPENDENT PORTION OF ENTRY
   
   150 loc = 150
   CALL READ(*9002,*9003,scr2,buf(1),11,noeor,nwords)
   npts = buf(11)
   npts4 = 4*npts
   
!  CHECK MCSID OF ENTRY TO BE SAME AS ONE OF THIS PASS.
   
   IF( buf(1) == mcsid ) GO TO 170
   
!  NO IT IS NOT THUS SKIP BALANCE OF ENTRY.
   
   loc = 170
   CALL READ(*9002,*9003,scr2,0,-npts4,noeor,nwords)
   GO TO 150
   
!  YES, THIS ENTRY IS OF CURRENT PASS MCSID. ADD POINT DATA TO CORE.
!  FIRST OUTPUT SIGMAS TO SCR3
   
   170 CALL WRITE( scr3, buf(2), 6, noeor )
   z(jindep  ) = rbuf(8)
   z(jindep+1) = rbuf(9)
   z(jindep+2) = rbuf(10)
   jindep = jindep + 3
   
!  INDEPENDENT POINTS NOT YET IN CORE ARE ADDED.
   
   CALL READ(*9002,*9003,scr2,buf(1),npts4,noeor,nwords)
   k = npts
   DO  j = 1,npts
     
!  CHECK IF EXTERNAL ID IS IN TABLE YET.
     
     IF( ndep < idep ) GO TO 220
     DO  l = idep,ndep,4
       IF( buf(j) == iz(l) ) GO TO 290
     END DO
     
!  NOT YET IN THUS ADD IT TO TABLE
     
     220 icrq = ndep + 4 - jcore
     IF( ndep+4 > jcore ) GO TO 9008
     iz(ndep+1) = buf(j)
     z(ndep+2) = rbuf(k+1)
     z(ndep+3) = rbuf(k+2)
     z(ndep+4) = rbuf(k+3)
     ndep = ndep + 4
     
     290 k = k + 3
     
   END DO
   
 END DO
!*****
!  ALL DATA FOR CURRENT MCSID HAS BEEN COLLECTED FROM SCR2.
!*****
 CALL CLOSE( scr2, clsrew )
 CALL CLOSE( scr3, clsrew )
 
!  DEPENDENT COORDINATES ARE SORTED ON EXTERNAL-ID.
 
 CALL sort( 0, 0, 4, 1, z(idep), ndep-idep+1 )
!*****
!  CONVERSION OF INDEPENDENT AND DEPENDENT POINTS TO LOCAL
!  MATERIAL COORDINATE SYSTEM. FIRST GET CSTM DATA TO USE.
!*****
 loc = 400
 CALL bisloc(*9000,mcsid,iz(icstm),14,cstms,jp)
 ivmat = icstm + jp + 1
 itran = ivmat + 3
 ictype = iz(ivmat-1)
 
!  FOR EACH POINT
!                               T
!                   (R     )=( T ) ( R     - V     )
!                     LOCAL           BASIC   MCSID
 
!                    (3X1)   (3X3)   (3X1)   (3X1)
 
 DO  i = iindep,nindep,3
   vec(1) = z(i  ) - z(ivmat  )
   vec(2) = z(i+1) - z(ivmat+1)
   vec(3) = z(i+2) - z(ivmat+2)
   CALL gmmats( z(itran),3,3,1, vec(1),3,1,0, z(i) )
 END DO
 
 DO  i = idep,ndep,4
   vec(1) = z(i+1) - z(ivmat  )
   vec(2) = z(i+2) - z(ivmat+1)
   vec(3) = z(i+3) - z(ivmat+2)
   CALL gmmats( z(itran),3,3,1, vec(1),3,1,0, z(i+1) )
 END DO
!*****
!  CONVERSION OF INDEPENDENT POINT LOCAL COORDINATES TO MAPPING
!  COORDINATES. (IF MCSID IS A RECTANGULAR SYSTEM THEN NO CHANGE.)
!*****
 loc = 490
 IF( ictype < 1 .OR. ictype > 3 ) GO TO 9000
 SELECT CASE ( ictype )
   CASE (    1)
     GO TO  589
   CASE (    2)
     GO TO 510
   CASE (    3)
     GO TO 530
 END SELECT
 
!  CYLINDRICAL COORDINATES
 
 510 avgl = 0.0
 DO  i = iindep,nindep,3
   vec(1) = SQRT( z(i)**2 + z(i+1)**2 )
   avgl = avgl + vec(1)
   IF( vec(1) <= 0.0 ) GO TO 515
   z(i+1) = ATAN2( z(i+1), z(i) )
   GO TO 517
   515 z(i+1) = 0.0
   517 z(i) = vec(1)
 END DO
 avgl = avgl / FLOAT(indpts)
 GO TO 589
 
!  SPHERICAL COORDINATES
 
 530 avgl = 0.0
 DO  i = iindep,nindep,3
   xsqysq = z(i)**2 + z(i+1)**2
   fl = SQRT( xsqysq )
   vec(1) = SQRT( xsqysq + z(i+2)**2 )
   avgl = avgl + vec(1)
   IF( vec(1) > 0.0 ) GO TO 540
   vec(2) = 0.0
   GO TO 550
   540 vec(2) = ATAN2( fl, z(i+2) )
   550 IF( fl > 0.0 ) GO TO 560
   vec(3) = 0.0
   GO TO 570
   560 vec(3) = ATAN2( z(i+1), z(i) )
   570 z(i  ) = vec(1)
   z(i+1) = vec(2)
   z(i+2) = vec(3)
 END DO
 avgl = avgl / FLOAT(indpts)
!*****
!  CONVERSION OF DEPENDENT POINT LOCAL COORDINATES TO MAPPING
!  COORDINATES.
!  (IF MCSID IS RECTANGULAR SYSTEM THEN NO CHANGE.)
!*****
 589 SELECT CASE ( ictype )
   CASE (    1)
     GO TO 609
   CASE (    2)
     GO TO 590
   CASE (    3)
     GO TO 600
 END SELECT
 
!  CYLINDRICAL COORDINATES
 
 590 DO  i = idep,ndep,4
   vec(1) = SQRT( z(i+1)**2 + z(i+2)**2 )
   IF( vec(1) <= 0.0 ) GO TO 592
   z(i+2) = ATAN2( z(i+2), z(i+1) )
   GO TO 593
   592 z(i+2) = 0.0
   593 z(i+1) = vec(1)
 END DO
 GO TO 609
 
!  SPHERICAL COORDINATES
 
 600 DO  i = idep,ndep,4
   xsqysq = z(i+1)**2 + z(i+2)**2
   fl = SQRT( xsqysq )
   vec(1) = SQRT( xsqysq + z(i+3)**2 )
   IF( vec(1) > 0.0 ) GO TO 602
   vec(2) = 0.0
   GO TO 604
   602 vec(2) = ATAN2( fl, z(i+3) )
   604 IF( fl > 0.0 ) GO TO 605
   vec(3) = 0.0
   GO TO 606
   605 vec(3) = ATAN2( z(i+2), z(i+1) )
   606 z(i+1) = vec(1)
   z(i+2) = vec(2)
   z(i+3) = vec(3)
 END DO
 
!  SET MAXIMUM AND MIMIMUM X,Y,Z VALUES.
 
 609 DO  i = 1,3
   vmax(i) = z(iindep+i-1)
   vmin(i) = z(iindep+i-1)
 END DO
 
 DO  i = iindep,nindep,3
   DO  j = 1,3
     vmax(j) = AMAX1( z(i+j-1), vmax(j) )
     vmin(j) = AMIN1( z(i+j-1), vmin(j) )
   END DO
 END DO
 
!  SET THE X,Y,Z RANGES
 
 DO  i = 1,3
   vmax(i) = vmax(i) - vmin(i)
   vec(i) = vmax(i)
 END DO
 
 IF( ictype == 1 ) GO TO 680
 vmax(2) = avgl * vmax(2)
 IF( ictype == 2 ) GO TO 680
 vmax(3) = avgl * vmax(3)
 
!  DIRECTION YIELDING MINIMUM RANGE DETERMINES PROJECTION
 
 680 IF( vmax(1) < vmax(2) ) GO TO 700
 IF( vmax(2) < vmax(3) ) GO TO 690
 685 k1 = 1
 k2 = 2
 kctype = 3
 GO TO 710
 690 k1 = 1
 k2 = 3
 kctype = 2
 GO TO 710
 700 IF( vmax(3) < vmax(1) ) GO TO 685
 k1 = 2
 k2 = 3
 kctype = 1
 
 710 xrange = vec(k1)
 yrange = vec(k2)
 IF( xrange  == 0.0) THEN
   GO TO   711
 ELSE
   GO TO   712
 END IF
 711 xrange = 1.0
 712 IF( yrange  == 0.0) THEN
   GO TO   713
 ELSE
   GO TO   714
 END IF
 713 yrange = 1.0
 
!  COORDINATES -K1- AND -K2- WILL BE KEPT.
 
!  TABLE OF INDEPENDENT AND DEPENDENT POINTS ARE REDUCED TO
!  TABLES OF X,Y PAIRS. FIRST TO GAIN SOME CORE, EXTERNAL
!  IDS- ARE WRITTEN TO SCR4.
 
 714 FILE = scr4
 loc = 714
 CALL OPEN(*9001,scr4,iz(ibuf1),wrtrew)
 DO  i = idep,ndep,4
   CALL WRITE( scr4, iz(i), 1, noeor )
 END DO
 CALL CLOSE( scr4, clsrew )
 
!  REDUCE INDEPENDENT POINTS TO XY PAIRS, SCALE BY X AND Y RANGES
!  RESPECTIVELY, AND COMPRESS IN CORE.
 
 j = iindep
 DO  i = iindep,nindep,3
   z(j  ) = z(i+k1-1) / xrange
   z(j+1) = z(i+k2-1) / yrange
   j = j + 2
 END DO
 nindep = j - 1
 
!  REDUCE DEPENDENT POINTS LIST. (J IS STILL GOOD)
 
 DO  i=idep,ndep,4
   z(j  ) = z(i+k1) / xrange
   z(j+1) = z(i+k2) / yrange
   j = j + 2
 END DO
 idep = nindep + 1
 ndep = j - 1
 depts = (ndep - idep + 1) / 2
!*****
!  INDEPENDENT AND DEPENDENT POINT COORDINATE LISTS ARE NOW
!  COMPLETE.  CALL FOR INTERPOLATION.
!*****
 CALL curvit( z(iindep), indpts, z(idep), depts, scr5,  &
     z(ndep+1), iz(ndep+1), lcore-ndep-1, ip2, 15.0, mcsid, xrange, yrange )
 
!  BRING -OES1M- SIGMAS INTO CORE FOR CURRENT -MCSID- PASS.
 
 isigma = iindep + 1
 nsigma = iindep + 6*indpts
 jsigma = isigma - 7
 icrq = nsigma - ibuf3
 IF( nsigma >= ibuf3 ) GO TO 9008
 FILE = scr3
 loc = 800
 CALL OPEN(*9001,scr3,iz(ibuf1),rdrew)
 CALL READ(*9002,*810,scr3,iz(isigma),ibuf3-isigma,noeor,nwords)
 loc = 810
 GO TO 9000
 
 810 IF( nwords /= 6*indpts ) GO TO 9000
 CALL CLOSE( scr3, clsrew )
 
!    (SIGMAS                ) = (G)(SIGMAS                        )
!           DEPENDENT POINTS              OES1M INDEPENDENT POINTS
 
!  SINCE THE ORDER OF THE ROWS IN THE G MATRIX ARE IN SORTED EXTERNAL
!  GRID ORDER EACH OUTPUT LINE OF OES1G WILL BE HANDLED ON ITS
!  OWN. THIS ELIMINATES NECESSITY OF HOLDING ANOTHER SIGMA ARRAY
!  IN CORE.
 
 FILE = oes1g
 loc = 815
 CALL OPEN(*9001,oes1g,iz(ibuf1),wrt)
 
!  OUTPUT ID RECORD. PREVIOUSLY PREPARED.
 
 idrec(3) = idrec(3) + 2000
 CALL WRITE( oes1g, idrec(1), 146, eor )
 mcb(1) = oes1g
 CALL wrttrl( mcb(1) )
 any1g = .true.
 
!  OPEN SCR5 CONTAINING ROWS OF THE G-MATRIX.
 
 FILE = scr5
 CALL OPEN(*9001,scr5,iz(ibuf3),rdrew)
 CALL fwdrec(*9002,scr5)
 
!  OPEN SCR4 CONTAINING LIST OF EXTERNAL IDS )
 
 FILE = scr4
 CALL OPEN(*9001,scr4,iz(ibuf2),rdrew)
 
!  COMPUTE AND OUTPUT SIGMAS FOR THE DEPENDENT POINTS
 
 buf(2) = mcsid
 DO  i=1,depts
   
!  READ THE EXTERNAL ID
   
   FILE = scr4
   CALL READ(*9002,*9003,scr4,buf(1),1,noeor,nwords)
   FILE = scr5
   
!  INITIALIZE SIGMAS(DEPENDENT POINT) TO ZERO
   
   DO  j = 3,8
     rbuf(j) = 0.0
   END DO
   
   k = 0
   loc = 825
   
!  READ ACTIVE INDEX AND G-VALUE FROM SCRATCH 5
   
   825 CALL READ(*9002,*840,scr5,rbuf(11),2,noeor,nwords)
   k = k + 10
   idx = jsigma + 6*buf(11)
   DO  j = 1,6
     rbuf(j+2) = rbuf(j+2) + rbuf(12)*z(idx+j)
   END DO
   GO TO 825
   
!  IF THERE WERE ANY G-VALUES THEN NOW COMPLETE THE OUTPUT LINE.
   
   840 IF( k <= 0 ) CYCLE
   
   buf(10) = k + kctype
   
   rbuf(11) = rbuf(6)
   rbuf(12) = rbuf(7)
   rbuf(13) = rbuf(8)
   
!  COMPUTE INVARIANTS FOR EACH LINE
   
   CALL curvps( rbuf( 3), rbuf( 6) )
   CALL curvps( rbuf(11), rbuf(14) )
   IF( .NOT. strain ) GO TO 881
   rbuf(5) = 2.0 * rbuf(5)
   rbuf(9) = 2.0 * rbuf(9)
   rbuf(13) = 2.0 * rbuf(13)
   rbuf(17) = 2.0 * rbuf(17)
   
!  APPEND DEVICE CODE TO EXTERNAL ID AND OUTPUT LINE
   
   881 buf(1) = 10*buf(1) + device
   CALL WRITE( oes1g, buf(1), 17, noeor )
 END DO
 
 CALL WRITE( oes1g, 0, 0, eor )
 IF(eofos1 .AND. jmcsid+2 > nmcsid) CALL CLOSE(oes1g,clsrew)
 CALL CLOSE( oes1g, cls )
 CALL CLOSE( scr4, clsrew )
 CALL CLOSE( scr5, clsrew )
!*****
!  ALL INDEPENDENT POINTS OUTPUT TO OES1G FOR 1 ACTIVE MCSID OF
!  CURRENT SUBCASE. GO TO NEXT MCSID.
!*****
 980 jmcsid = jmcsid + 2
 IF( jmcsid <= nmcsid ) GO TO 100
!*****
!  ALL THROUGH FORMING OES1G FOR CURRENT SUBCASE.
!*****
 5000 RETURN
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
END SUBROUTINE curv3
