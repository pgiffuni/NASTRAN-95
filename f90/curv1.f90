SUBROUTINE curv1
    !*****
    !  INITIALIZATION OVERLAY. ALL LOGIC INDEPENDENT OF PROCESSING
    !  THE SUBCASE DATA ON OES1 IS HANDLED IN THIS INITIALIZATION
    !  ROUTINE OF THE -CURV- MODULE
 
    !     OPEN CORE MAP DURING -CURV1- EXECUTION.
    !     =======================================
    !       INITIAL               AFTER CURV1 RETURNS
    !     +-----------+             +----------------+
    !     I Z(IELTYP) I             I  Z(IELTYP)     I
    !     I    .      I             I     .          I
    !     I  ELEMENT  I             I  REDUCED       I
    !     I  TYPES    I             I  ELEMENT-TYPES I
    !     I  BEING    I             I  LIST          I
    !     I  PLACED   I             I     .          I
    !     I  IN SCR1  I             I  Z(NELTYP)     I
    !     I    .      I             +----------------+
    !     I Z(NELTYP) I             I  Z(IMCSID)     I
    !     +-----------+             I     .          I
    !     I Z(IMID)   I             I  MCSID LIST    I
    !     I    .      I             I  OF MCSIDS     I
    !     I MATID-    I             I  ACTUALLY      I
    !     I MCSID-    I             I  REFERENCED    I
    !     I FLAG-     I             I     .          I
    !     I ENTRIES   I             I  Z(NMCSID)     I
    !     I    .      I             +----------------+
    !     I Z(NMID)   I             I  Z(ICSTM)      I
    !     +-----------+             I     .          I
    !     I Z(ISIL)   I             I  CSTMS IN      I
    !     I    .      I             I  EXISTENCE     I
    !     I SILS IN   I             I  FOR MCSIDS    I
    !     I INTERNAL  I             I  IN ABOVE      I
    !     I SORT      I             I  TABLE         I
    !     I    .      I             I     .          I
    !     I Z(NSIL)   I             I  Z(NCSTM)      I
    !     +-----------+             +----------------+
    !     I Z(IEXT)   I             I     .          I
    !     I    .      I             I  AVAILABLE     I
    !     I EXTERNAL  I             I  CORE          I
    !     I IDS IN    I             I     .          I
    !     I INTERNAL  I             I     .          I
    !     I SORT      I             I     .          I
    !     I    .      I             I     .          I
    !     I Z(NEXT)   I             I     .          I
    !     +-----------+             I     .          I
    !     I   .       I             I     .          I
    !     I AVAILABLE I             I     .          I
    !     I CORE      I             I     .          I
    !     I   .       I             I     .          I
    !     I Z(JCORE)  I             I  Z(JCORE)      I
    !     +-----------+             +----------------+
    !     I Z(IBUF4)  I             I  Z(IBUF4)      I
    !     I Z(IBUF3)  I             I  Z(IBUF3)      I
    !     I Z(IBUF2)  I             I  Z(IBUF2)      I
    !     I Z(IBUF1)  I             I  Z(IBUF1)      I
    !     I GINO-BUFS I             I  GINO-BUFS     I
    !     I Z(LCORE)  I             I  Z(LCORE)      I
    !     +-----------+             +----------------+
 
    !*****
    REAL :: z(1)     ,rbuf(100)
 
    INTEGER :: cstms    ,scr1     ,scr2     ,scr3
    INTEGER :: scr4     ,oes1m    ,oes1g    ,oes1     ,scr5
    INTEGER :: cstm     ,est      ,sil      ,gpl
    INTEGER :: eltype   ,subcas   ,file     ,estwds
    INTEGER :: ewords   ,owords   ,depts    ,cstype
    INTEGER :: device   ,oldid    ,buf      ,sbuf
    INTEGER :: rd       ,rdrew    ,wrt      ,wrtrew
    INTEGER :: cls      ,clsrew   ,eor      ,sysbuf
 
    INTEGER :: mat(6)   ,elem(5,4)
 
    LOGICAL :: any      ,eofos1   ,first    ,anyout
    LOGICAL :: foes1g   ,strain   ,any1m    ,any1g
 
    COMMON/BLANK /     ip1      ,ip2      ,icmplx   ,zdum(3)
 
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
 
    DATA mat / 103,1,12,   203,2,17  /
 
    !        - - - - - - - - CURV-MODULE ELEMENTS DATA - - - - - - - -
 
    !                   ELEMENT   EST       CONNECT.  MATID     BGPDT
    !                   TYPE      WORDS     POINTS    INDEX     INDEX
    !                   =======   =======   =======   =======   =======
    !  TRIA1
    DATA elem /           6        ,27       ,3        ,6        ,15  &
        !  TRIA2  &
        ,17       ,21       ,3        ,6        ,9   &
        !  QUAD1  &
        ,19       ,32       ,4        ,7        ,16  &
        !  QUAD2  &
        ,18       ,26       ,4        ,7        ,10       /
 
    !  IF EITHER OF THESE PARAMS IS EXCEEDED RESET AND RE-DIMENSION
    !  SBUF OR BUF.
 
    lsbuf  = 10
    lbuf   = 100
    nelems = 4
    logerr = 37
    !*****
    !  INITIALIZATION OF CORE AND FLAGS
    !*****
    foes1g = .true.
    IF (ip1 > 0) foes1g = .false.
    any1m = .false.
    any1g = .false.
    lmcsid = 0
 
    lcore = korsz( iz(1) )
    DO  i = 1,lcore
        iz(i) = 0
    END DO
    ibuf1 = lcore - sysbuf
    ibuf2 = ibuf1 - sysbuf
    ibuf3 = ibuf2 - sysbuf
    ibuf4 = ibuf3 - sysbuf
 
    !  SET FILE NUMBERS EXPLICITYLY.  ALL OVERLAYS REFERENCE /CURVTB/
 
    oes1  = 101
    mpt   = 102
    cstm  = 103
    est   = 104
    sil   = 105
    gpl   = 106
    oes1m = 201
    oes1g = 202
    scr1  = 301
    scr2  = 302
    scr3  = 303
    scr4  = 304
    scr5 = 305
    jcore = ibuf4 - 1
    FILE = 0
    loc = 300
    icrq =-ibuf4
    IF( ibuf4  > 0) THEN
        GO TO   300
    ELSE
        GO TO  9008
    END IF
    !*****
    !  ALLOCATE TABLE OF ELEMENT TYPES PLACED ON ESTX(SCR1).  MAXIMUM
    !  SIZE NOW AND REDUCED LATER TO ACTUAL SIZE.
    !*****
300 ieltyp = 1
    jeltyp = ieltyp
    neltyp = nelems
    !*****
    !  CONSTRUCTION OF TABLE CONTAINING ENTRIES OF,
 
    !     MID   = MATERIAL-ID
    !     MCSID = MATERIAL-COORDINATE-SYSTEM-ID
    !     FLAG  = REFERENCE-FLAG
 
    !  ALL MAT1 AND MAT2 BULK DATA CARDS CONTAINING A NON-ZERO -MCSID-
    !  RESULT IN AN ENTRY BEING ADDED TO THIS TABLE. TABLE IS THEN SORTED
    !  ON -MID-.
    !*****
    imid = neltyp + 1
    nmid = imid - 1
 
    !  OPEN MPT USING -PRELOC- FUNCTION.
 
    FILE = mpt
    loc = 400
    CALL preloc(*9001,iz(ibuf1),mpt)
 
    !  PASS MAT1 AND MAT2 DATA IF ANY.
 
    DO  i = 1,6,3
        iwords = mat(i+2)
        IF( iwords > lbuf ) GO TO 9000
        CALL locate(*480,iz(ibuf1),mat(i),idum)
410     CALL READ(*9002,*480,mpt,buf(1),iwords,noeor,nwords)
        IF( buf(iwords) <= 0 ) GO TO 410
        icrq = nmid + 3 - jcore
        IF( nmid+3 > jcore ) GO TO 9008
        iz(nmid+1) = buf(1)
        iz(nmid+2) = buf(iwords)
        iz(nmid+3) = 0
        nmid = nmid + 3
        GO TO 410
   
    !  EOR HIT READING MAT1 OR MAT2 CARDS
480 CONTINUE
    END DO
 
    !  TABLE COMPLETE, THUS NOW SORT IT. IF TABLE IS EMPTY WE ARE THROUGH
 
    CALL CLOSE( mpt, clsrew )
    lmid = nmid - imid + 1
    nmids = lmid / 3
    loc = 570
    IF (lmid < 0) THEN
        GO TO  9000
    ELSE IF (lmid == 0) THEN
        GO TO   950
    END IF
570 CALL sort( 0, 0, 3, 1, iz(imid), lmid )
    !*****
    !  LOAD LIST OF SILS INTO CORE, FOLLOWED BY LIST OF EXTERNAL IDS.
    !  THIS IS REQUIRED ONLY IF OES1G IS TO BE FORMED.
    !*****
    IF( .NOT. foes1g ) GO TO 630
    FILE = sil
    loc = 580
    isil = nmid + 1
    CALL gopen( sil, iz(ibuf1), 0 )
    CALL READ(*9002,*580,sil,iz(isil),jcore-isil,noeor,lsil)
    icrq = jcore - isil
    GO TO 9008
 
580 nsil = isil + lsil - 1
    CALL CLOSE( sil, clsrew )
 
    FILE = gpl
    loc = 590
    iext = nsil + 1
    CALL gopen( gpl, iz(ibuf1), 0 )
    CALL READ(*9002,*590,gpl,iz(iext),jcore-iext,noeor,lext)
    icrq = jcore - iext
    GO TO 9008
 
590 next = iext + lext - 1
    CALL CLOSE( gpl, clsrew )
    IF( lsil /= lext ) GO TO 9000
    !*****
    !  EST IS NOW READ. ANY ELEMENTS IN THE EST WHOSE MATERIAL ID REFERENCES
    !  A MAT1 OR MAT2 ENTRY WHICH CONTAINS A NON-ZERO MATERIAL-COORDINATE-
    !  SYSTEM-ID, WILL BE PLACED IN AN ABBREVIATED EST ON SCRATCH1.
 
    !  FORMAT OF EACH ELEMENT TYPE RECORD.
 
    !         ELEMENT TYPE NUMBER
    !         NUMBER OF WORDS PER EACH OF THE FOLLOWING ENTRIES.
    !         NUMBER OF POINTS PER THIS ELEMENT TYPE.
 
    !        * ELEMENT-ID
    !       *  MCSID = MATERIAL-COORDINATE-SYSTEM-ID
    ! ENTRY*
    !       *  EXTERNAL-GRID-IDS THIS ELEMENT CONNECTS (1 OR MORE)
    !        * X,Y,Z BASIC COORDINATE SETS OF EACH CONNECTED POINT(1 OR MORE
 
    !           ( ABOVE ELEMENT ENTRY REPEATS FOR EACH ELEMENT
    !             REFERENCING A MAT1 OR MAT2 CARD HAVING A NON-ZERO MCSID.)
 
    !*****
630 loc = 630
    FILE = scr1
    CALL OPEN(*9001,scr1,iz(ibuf2),wrtrew)
    FILE = est
    CALL gopen(  est, iz(ibuf1), 0 )
 
    oldid = -99999998
 
    !  READ ELEMENT TYPE OF NEXT EST RECORD AND DETERMINE IF IT IS
    !  AMONG ELEMENT TYPES TO BE EVEN CONSIDERED.
 
645 loc = 645
    CALL READ(*800,*9003,est,eltype,1,noeor,nwords)
    DO  i = 1,nelems
        IF( eltype == elem(1,i) ) GO TO 670
    END DO
    CALL fwdrec(*9002,est)
    GO TO 645
 
    !  OK THIS  ELEMENT TYPE RECORD IS TO BE CONSIDERED.
 
670 estwds = elem(2,i)
    loc = 670
    IF( estwds > lbuf ) GO TO 9000
    any    = .false.
    npts   = elem(3,i)
    imatid = elem(4,i)
    ixyz1 = elem(5,i)
    ixyz2 = ixyz1 + 4*npts - 1
    k1 = 2 + npts
    loc = 680
    IF( k1 > lsbuf ) GO TO 9000
 
    !  READ AN ELEMENT ENTRY AND CHECK TO DETERMINE IF IT IS TO BE USED.
 
690 loc = 690
    CALL READ(*9002,*780,est,buf(1),estwds,noeor,nwords)
    matid = buf(imatid)
    IF( matid == oldid ) GO TO 730
    CALL bisloc(*690,matid,iz(imid),3,nmids,jp)
    mcsid = iz(imid+jp)
    oldid = matid
    iz(imid+jp+1) = 7
 
    !  DEVELOP AND OUTPUT ABBREVIATED ENTRY TO SCRATCH1.
    !  (INITIALIZE RECORD WITH THREE-WORD-HEADER ENTRY.)
 
730 IF( any ) GO TO 733
    sbuf(1) = eltype
    sbuf(2) = 4*npts + 2
    sbuf(3) = npts
    CALL WRITE( scr1, sbuf(1), 3, noeor )
    iz(jeltyp) = eltype
    jeltyp = jeltyp + 1
    any = .true.
 
733 sbuf(1) = buf(1)
    sbuf(2) = mcsid
 
    !  CONVERT SILS TO EXTERNAL-IDS IF OES1G IS TO BE BUILT
 
    IF( foes1g ) GO TO 740
    DO  i=3,k1
        sbuf(i) = 0
    END DO
    GO TO 760
 
740 jsil = 2
    loc = 740
    DO  i = 3,k1
        CALL bisloc(*9000,buf(jsil),iz(isil),1,lsil,jp)
        sbuf(i) = iz(iext+jp-1)
        jsil = jsil + 1
    END DO
 
    !  OUTPUT THIS PORTION OF ENTRY AND THEN XYZ COMPONENTS OF CONNECTED
    !  POINTS
 
760 CALL WRITE( scr1, sbuf(1), npts+2, noeor )
 
    DO  i = ixyz1,ixyz2,4
        CALL WRITE( scr1, buf(i+1), 3, noeor )
    END DO
 
    !  GO FOR NEXT ELEMENT OF THIS TYPE
 
    GO TO 690
 
    !  END OF ENTRIES FUR CURRENT ELEMENT TYPE.
 
780 loc = 780
    IF( nwords /= 0 ) GO TO 9000
    IF( any ) CALL WRITE( scr1, 0, 0, eor )
    GO TO 645
 
    !  END OF ALL ELEMENT TYPES IN EST
 
800 CALL CLOSE(  est, clsrew )
    CALL CLOSE( scr1, clsrew )
    !*****
    !  REDUCTION OF MATERIAL-ID AND COORDINATE-SYSTEM-ID TO THOSE
    !  ACTUALLY REFERENCED BY ELEMENTS BEING CONSIDERED.
    !*****
    neltyp = jeltyp - 1
 
    !  RESORT MID-MCSID TABLE ON MCSID.
 
    CALL sort( 0, 0, 3, 2, iz(imid), lmid )
    imcsid = neltyp + 1
    nmcsid = neltyp
    loc = 820
    oldid = 0
    DO  i = imid,nmid,3
        IF( iz(i+2)  < 0) THEN
            GO TO  9000
        ELSE IF ( iz(i+2)  == 0) THEN
            GO TO   840
        END IF
   
        !  ELIMINATE DUPLICATE MCSIDS.
   
820     IF( iz(i+1) == oldid ) CYCLE
        oldid = iz(i+1)
        nmcsid = nmcsid + 2
        iz(nmcsid-1) = iz(i+1)
        iz(nmcsid  ) = 0
840 CONTINUE
    END DO
    lmcsid = nmcsid - imcsid + 1
    mcsids = lmcsid / 2
 
    !  IF TABLE IS NOW EMPTY THERE IS NOTHING MORE TO DO
 
    loc = 860
    IF (lmcsid < 0) THEN
        GO TO  9000
    ELSE IF (lmcsid == 0) THEN
        GO TO   950
    END IF
    !*****
    !  COORDINATE SYSTEMS WHICH MAY BE REFERENCED ARE AT THIS TIME
    !  PULLED INTO CORE FROM THE -CSTM- DATA BLOCK. (SILS AND EXTERNAL-IDS
    !  TABLES IF IN CORE ARE NO LONGER REQUIRED.)
    !*****
860 icstm = nmcsid + 1
    ncstm = nmcsid
    FILE = cstm
    CALL gopen( cstm, iz(ibuf1), 0 )
870 icrq = ncstm + 14 - jcore
    IF( ncstm+14 > jcore ) GO TO 9008
880 CALL READ(*9002,*900,cstm,iz(ncstm+1),14,noeor,nwords)
    kid = iz(ncstm+1)
    CALL bisloc (*880, kid, iz(imcsid), 2, mcsids, jp)
    ncstm = ncstm + 14
    GO TO 870
 
    !  END OF COORDINATE SYSTEM DATA
 
900 CALL CLOSE( cstm, clsrew )
    lcstm = ncstm - icstm + 1
    cstms = lcstm / 14
    CALL sort( 0, 0, 14, 1, z(icstm), lcstm )
    CALL pretrs( iz(icstm), lcstm )
    !*****
    !  INITIALIZE INPUT AND OUTPUT FILE POSITIONS.
    !*****
950 CALL gopen (oes1, iz(ibuf1), 0)
 
    !  CHECK FOR STRAIN OPTION
 
    FILE = oes1
    loc = 910
    CALL READ(*9002,*9003,oes1,idrec(1),2,0,flag)
    i = idrec(2)
    IF (i /= 5.AND.i /= 21.AND.i /= 1005) GO TO 9000
    strain = .false.
    IF (i == 21) strain = .true.
    icmplx = 0
    IF (i == 1005) icmplx = 1
    CALL bckrec (oes1)
 
    CALL CLOSE( oes1, cls )
    eofos1 = .false.
 
    CALL gopen( oes1m, iz(ibuf1), 1 )
    CALL CLOSE( oes1m, cls )
 
    IF( .NOT. foes1g ) GO TO 5000
    CALL gopen( oes1g, iz(ibuf1), 1 )
    CALL CLOSE( oes1g, cls )
    !*****
    !  END OF INITIALIZATION
    !*****
5000 RETURN
    !*****
    !  ERROR CONDITION ENCOUNTERED.
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

END SUBROUTINE curv1
