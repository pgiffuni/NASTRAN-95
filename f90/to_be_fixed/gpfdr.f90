SUBROUTINE gpfdr
     
!     GRID-POINT-FORCE-DATA-RECOVERY (MODULE)
 
!     THIS MODULE FORMULATES OFP TYPE OUTPUT DATA BLOCKS OF ELEMENT-
!     STRAIN ENERGYS AND GRID-POINT FORCE BALANCES.
 
!     DMAP CALLING SEQUENCES.
 
!     SOLUTION 1 -
!     GPFDR  CASECC,UGV,KMAT,KDICT,ECT,EQEXIN,GPECT,PG,QG/ONRGY1,OGPF1/
!            *STATICS* $
!     SOLUTION 3 -
!     GPFDR  CASECC,PHIG,KMAT,KDICT,ECT,EQEXIN,GPECT,LAMA,/ONRGY1,OGPF1/
!            *REIG* $
 
!     COMMENT FROM G.CHAN/UNISYS, 1/88 -
!     FOR MACHINES OF 32 OR 36 BIT WORDS, THE STRAIN ENERGY COMPUTATION
!     (OTHER COMPUTATIONS TOO) MUST BE DONE IN DOUBLE PRECISION. SINCE
!     THE K-MATRIX NORMALLY IN 10**7, AND THE DISPLACEMENT VECTOR IN
!     10**-2 OR 10**-3 RANGE, SINGLE PRECISION COMPUTATION GIVES BAD
!     RESULT.
 
 LOGICAL :: dicout   ,engout   ,enflag   ,anygp    ,diagm    ,any   ,  &
     DOUBLE   ,silin    ,enfile   ,gpfile   ,eorst4   ,axic  , axif
 INTEGER :: z        ,casecc   ,scrt1    ,eor      ,sysbuf   ,title ,  &
     names(2) ,typout   ,ug       ,scrt2    ,core     ,subr(2)  &
     ,        gsize    ,ect      ,scrt3    ,symflg   ,extgp    ,subtit,  &
     buf(100) ,trl(7)   ,onrgy1   ,gpset    ,points   ,scrt4 ,  &
     pg       ,qg       ,ugpgqg   ,oload(2) ,ospcf(2) ,iii(2),  &
     isum(10) ,scale(2) ,kvec(10) ,clseof   ,recidx(3),outpt ,  &
     FILE     ,mcb(7)   ,eqexin   ,ogpf1    ,elnset   ,branch,  &
     rd       ,app      ,gpect    ,set      ,gpdvis   ,subcas,  &
     rdrew    ,ectwds   ,grdpts   ,comps(32),eldvis   ,buf1  ,  &
     wrt      ,eltype   ,grid1    ,exelid   ,dicloc   ,buf2  ,  &
     wrtrew   ,elem     ,NAME(2)  ,gridl    ,idrec(10),buf3  ,  &
     cls      ,phead(3) ,recid(3) ,gpsil    ,comp     ,buf4  ,  &
     clsrew   ,estid    ,oldcod   ,out(10)  ,oldid    ,buf5  ,  &
     pivot    ,extid    ,ptr      ,entrys   ,total    ,buf6  , method(20)
 REAL :: rz(1)    ,rbuf(5)  ,rout(10) ,vec(6)   ,ridrec(146)     ,  &
     rsum(10) ,fvec(10)
 DOUBLE PRECISION :: diii     ,elengy   ,toteng   ,dz(1)
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON  /xmssg /   ufm      ,uwm      ,uim      ,sfm      ,swm
 COMMON  /system/   sysbuf   ,outpt
 COMMON  /names /   rd       ,rdrew    ,wrt      ,wrtrew  ,clsrew ,  &
     cls      ,clseof
 COMMON  /gpta1 /   nelems   ,last     ,incr     ,elem(1)
 COMMON  /unpakx/   typout   ,irow     ,nrow     ,incrx
 COMMON  /zntpkx/   a(4)     ,irowx    ,ieol
 COMMON  /zzzzzz/   z(1)
 COMMON  /BLANK /   app(2)
 EQUIVALENCE        (z(1),rz(1),dz(1)) ,(buf(1),rbuf(1)),  &
     (out(1),rout(1))   ,(name1,names(1)),  &
     (name2,names(2))   ,(idrec(1),ridrec(1)),  &
     (diii,iii(1))      ,(isum(1),rsum(1)), (kvec(1),fvec(1))
 DATA     enoeor,   eor / 0,1/   ,  lbuf/100/, subr/4HGPFD,4HR    /
 DATA     casecc,   ug, kmat,kdict,ect,eqexin,gpect,pg, qg        /  &
     101   ,   102,103, 104,  105,106,   107,  108,109       /
 DATA     onrgy1,   ogpf1,scrt1,scrt2,scrt3,scrt4  ,lama          /  &
     201   ,   202,  301,  302,  303,  304,    108           /
 DATA     meths /   10/, oload/4HAPP-,4HLOAD/, ospcf/4HF-of,4H-spc/
 DATA     scale /   5, 0/, isum  / 0,0,4H*tot,4HALS*,0,0,0,0,0,0  /
 DATA     method/   4HSTAT,4HICS , 4HREIG,4HEN  , 4HDS0 ,4H       ,  &
     4HDS1 ,4H    , 4HFREQ,4H    , 4HTRAN,4HSNT    ,  &
     4HBKL0,4H    , 4HBKL1,4H    , 4HCEIG,4HEN     , 4HPLA ,4H    /
 
!     CASE CONTROL POINTERS
 
 DATA     title ,   subtit, label        / 39, 71,103             /
 DATA     isym  ,   igp,ieln,ilsym,isubc / 16,167,170,200,1       /
 
!     DETERMINE APPROACH
 
 n = 2*meths - 1
 DO  i = 1,n,2
   IF (app(1) == method(i)) GO TO 40
 END DO
 WRITE  (outpt,30) uwm,app
 30 FORMAT (a25,' 2342, UNRECOGNIZED APPROACH PARAMETER ',2A4,  &
     ' IN GPFDR INSTRUCTION.')
 i = 19
 nerror = 0
 GO TO 1810
 
 40 branch = (i+1)/2
 
!     INITIALIZATION AND BUFFER ALLOCATION.
 
 core = korsz(z)
 buf1 = core - sysbuf - 2
 buf2 = buf1 - sysbuf - 2
 buf3 = buf2 - sysbuf - 2
 buf4 = buf3 - sysbuf - 2
 buf5 = buf4 - sysbuf - 2
 buf6 = buf5 - sysbuf - 2
 core = buf6 - 1
 
!     READ IN FREQUENCIES IF APPROACH IS REIGEN
 
 IF (branch /= 2) GO TO 70
 mode = 0
 CALL OPEN (*70,lama,z(buf1),rdrew)
 CALL fwdrec (*60,lama)
 CALL fwdrec (*60,lama)
 lfeq = core
 50 CALL READ (*60,*60,lama,buf,7,0,iwords)
 rz(core) = rbuf(5)
 core = core - 1
 GO TO 50
 60 CALL CLOSE (lama,clsrew)
 
!     GPTA1 DUMMY ELEMENT SETUP CALL.
 
 70 CALL delset
 nerror = 1
 IF (core > 0.0) THEN
   GO TO    80
 ELSE
   GO TO  1800
 END IF
 
!     OPEN CASE CONTROL
 
 80 FILE = casecc
 nerror = 2
 CALL OPEN (*1760,casecc,z(buf1),rdrew)
 CALL fwdrec (*1770,casecc)
 
!     OPEN VECTOR FILE.
 
 FILE = ug
 CALL OPEN (*1760,ug,z(buf2),rdrew)
 CALL fwdrec (*1770,ug)
 trl(1) = ug
 CALL rdtrl (trl)
 gsize = trl(3)
 
!     PREPARE OUTPUT BLOCKS FOR ANY OUTPUTS POSSIBLE
 
 enfile = .false.
 CALL OPEN (*90,onrgy1,z(buf3),wrtrew)
 enfile = .true.
 CALL fname (onrgy1,NAME)
 CALL WRITE (onrgy1,NAME,2,eor)
 CALL CLOSE (onrgy1,clseof)
 mcb(1) = onrgy1
 CALL rdtrl (mcb)
 mcb(2) = 0
 CALL wrttrl (mcb)
 
 90 gpfile = .false.
 nerror = 4
 CALL OPEN (*100,ogpf1,z(buf3),wrtrew)
 gpfile = .true.
 CALL fname (ogpf1,NAME)
 CALL WRITE (ogpf1,NAME,2,eor)
 CALL CLOSE (ogpf1,clseof)
 
 100 movepq = 1
 silin  = .false.
 trl(1) = eqexin
 CALL rdtrl (trl)
 points = trl(2)
 isilex = 1
 nsilex = 2*points
 nerror = 5
 IF (nsilex > core) GO TO 1800
 iccz = nsilex
 icc  = iccz + 1
 GO TO 120
 
!     OPEN CASECC AND UGV WITH NO REWIND
 
 110 FILE   = casecc
 nerror = 8
 CALL OPEN (*1760,casecc,z(buf1),rd)
 FILE = ug
 CALL OPEN (*1760,ug,z(buf2),rd)
 
!     READ NEXT CASE CONTROL RECORD.
 
 120 CALL READ (*1750,*130,casecc,z(iccz+1),core-iccz,eor,iwords)
 nerror = 7
 GO TO 1800
 
 130 ncc    = iccz + iwords
 itemp  = iccz + isubc
 subcas = z(itemp)
 
!     SYMMETRY-REPCASE, GP-FORCE REQUEST, AND EL-ENERGY REQUEST CHECKS
 
 itemp  = iccz + isym
 symflg = z(itemp)
 
!     SET REQUEST PARAMETERS FOR GP-FORCE AND EL-ENERGY.
 
 itemp  = iccz + igp
 gpset  = z(itemp)
 IF (.NOT.gpfile) gpset = 0
 gpdvis = z(itemp+1)
 itemp  = iccz + ieln
 elnset = z(itemp)
 IF (.NOT. enfile) elnset = 0
 eldvis = z(itemp+1)
 IF (gpset <= 0 .AND. elnset <= 0) GO TO 170
 
!     POINTERS TO SET LIST DOMAINS
 
 itemp = iccz + ilsym
 lsym  = z(itemp)
 itemp = itemp + lsym + 1
 140 set   = z(itemp)
 iset  = itemp + 2
 lset  = z(itemp+1)
 
!     CHECK IF THIS SET IS THE ONE FOR GP-FORCE
 
 IF (set /= gpset) GO TO 150
 igplst = iset
 lgplst = lset
 
!     CHECK IF THIS SET IS THE ONE FOR EL-ENERGY
 
 150 IF (set /= elnset) GO TO 160
 iellst = iset
 lellst = lset
 
 160 itemp = iset + lset
 IF (itemp < ncc) GO TO 140
 
!     IS THIS A REPCASE.  IF SO BACK-RECORD UG (REP-CASE OK ONLY FOR
!                                               STATICS)
 
 170 IF (symflg < 0.0) THEN
   GO TO   180
 ELSE
   GO TO   190
 END IF
 
!     NEGATIVE SYMFLG IMPLIES A REP-CASE.
 
 180 IF (app(1) /= method(1)) GO TO 120
 
!     REP-CASE AND STATICS APPROACH THUS POSITION BACK ONE
!     VECTOR ON UG UNLESS THERE IS NO REQUEST FOR GP-FORCE OR
!     EL-ENERGY TO BEGIN WITH.
 
 IF (gpset == 0 .AND. elnset == 0) GO TO 120
 CALL bckrec (ug)
 movepq = movepq - 1
 GO TO 210
 
!     NOT A REP-CASE BUT STILL IF THERE IS NO REQUEST FOR
!     GP-FORCE OR EL-ENERGY POSITION OVER VECTORS ASSOCIATED
!     WITH THIS CASE.
 
 190 IF (gpset /= 0 .OR. elnset /= 0) GO TO 210
 IF (symflg == 0.0) THEN
   GO TO   200
 ELSE
   GO TO   120
 END IF
 
!     NOT A SYMMETRY CASE (WHICH WOULD USE VECTORS ALREADY READ, THUS
!     SKIP A VECTOR ASSOCIATED WITH THIS CASE.
 
 200 nerror = 8
!IBMD 6/93 CALL FWDREC (*1770,UG)
!IBMNB 6/93
! MAJOR LOOP OF MODULE TERMINATES WITH ENDING OF CASE CONTROL OR
! END OF EIGENVECTORS COMPUTED.  IF MODES CARD IS USED AND SPECIFIES
! MORE MODES THAN WERE COMPUTED, THEN THE FOLLOWING WILL TERMINATE
! THE LOOP.  (SEE DEMO T03011A WHICH COMPUTED 4 EIGENVALUES BUT HAD
! A MODES CARD SPECIFYING 5 MODES)
 CALL fwdrec (*1750,ug)
!IBMNE
 movepq = movepq + 1
 GO TO 120
 
!  BRING VECTOR INTO CORE, BRANCH IF SYMMETRY CASE.
 
 210 ivec   = ncc + 1
 ivecz  = ncc
 nvec   = ivecz + gsize
 nerror = 9
 IF (nvec > core) GO TO 1800
 ASSIGN 320 TO iretrn
 ugpgqg = ug
 220 IF (symflg > 0.0) THEN
   GO TO   260
 END IF
 
 230 irow   = 1
 nrow   = gsize
 incrx  = 1
 typout = 1
 CALL unpack (*240,ugpgqg,rz(ivec))
 GO TO 310
 
!     NULL VECTOR (SET VECTOR SPACE TO ZERO)
 
 240 DO  i = ivec,nvec
   rz(i) = 0.0
 END DO
 GO TO 310
 
!     SYMMETRY SEQUENCE.  SUM VECTORS OF SEQUENCE APPLYING COEFFICIENTS.
 
 260 itemp = iccz + ilsym
 lsym  = z(itemp)
 
!     BACK UP OVER THE VECTORS OF THE SEQUENCE
 
 DO  i = 1,lsym
   CALL bckrec (ugpgqg)
 END DO
 
 DO  i = ivec,nvec
   rz(i) = 0.0
 END DO
 
 DO  i = 1,lsym
   itemp = itemp + 1
   coef  = rz(itemp)
   
!     SUM IN COEF*VECTOR(I)
   
   CALL intpk (*300,ugpgqg,0,1,0)
   290 CALL zntpki
   j = ivecz + irowx
   rz(j) = rz(j) + coef*a(1)
   IF (ieol == 0) THEN
     GO TO   290
   ELSE
     GO TO   300
   END IF
 END DO
 310 GO TO iretrn, (320,1460)
 
!     AT THIS POINT VECTOR IS IN CORE ALONG WITH THE CASE CONTROL RECORD
 
!     NOW START ECT PASS.  IN THIS PASS GP-FORCES REQUESTED WILL BE
!     WRITTEN TO PMAT (A SCRATCH SET ACTUALLY=SCRT1), AND BY THE GINO
!     DIRECT-ACCESS METHOD.  ALSO EL-ENERGY OUTPUTS WILL BE FORMED FOR
!     ANY REQUESTED ELEMENTS.
 
!     NOTE.  THE ASSEMBLY OF GP-FORCES FOR OUTPUT IS ACCOMPLISHED AFTER
!     ALL GP-FORCES REQUESTED HAVE BEEN WRITTEN TO PMAT.
 
 320 CALL CLOSE (casecc,cls)
 CALL CLOSE (ug,cls)
 IF (silin) GO TO 370
 
!     GET SECOND RECORD OF EQEXIN INTO CORE AND TRANSFER CODES FROM
!     SILS TO EXTERNALS AND THEN INSURE SORT ON SILS.
 
 nerror = 6
 FILE   = eqexin
 CALL OPEN (*1760,eqexin,z(buf1),rdrew)
 CALL fwdrec (*1770,eqexin)
 CALL fwdrec (*1770,eqexin)
 CALL READ (*1770,*350,eqexin,z(isilex),core-isilex,noeor,iwords)
 330 WRITE  (outpt,340) swm,eqexin
 340 FORMAT (a27,' 2343.  DATA BLOCK',i5,' IS EITHER NOT -EQEXIN- OR ',  &
     'POSSIBLY INCORRECT.')
 GO TO 1810
 
 350 IF (iwords /= 2*points) GO TO 330
 CALL CLOSE (eqexin,clsrew)
 DO  i = isilex,nsilex,2
   z(i  ) = 10*z(i) + MOD(z(i+1),10)
   z(i+1) = z(i+1)/10
 END DO
 silin  = .true.
 CALL sort (0,0,2,2,z(isilex),nsilex-isilex+1)
 
!     SET UP OFP ID RECORD WITH TITLE, SUBTITLE, AND LABEL.
 
 370 itit = iccz + title
 isub = iccz + subtit
 ilab = iccz + label
 DO  i = 1,32
   idrec(i+ 50) = z(itit)
   idrec(i+ 82) = z(isub)
   idrec(i+114) = z(ilab)
   itit = itit + 1
   isub = isub + 1
   ilab = ilab + 1
 END DO
 DO  i = 1,50
   idrec(i) = 0
 END DO
 FILE   = ect
 nerror = 10
 CALL OPEN (*1760,ect,z(buf4),rdrew)
 FILE = kmat
 CALL OPEN (*1760,kmat,z(buf5),rdrew)
 
!     DETERMINE PRECISION OF KMAT DATA
 
 mcb(1) = kmat
 CALL rdtrl (mcb)
 DOUBLE = .false.
 IF (mcb(2) == 2) DOUBLE = .true.
 FILE = kdict
 CALL OPEN (*1760,kdict,z(buf6),rdrew)
 CALL fwdrec (*1770,kdict)
 
!     PMAT WILL BE ON SCRATCH1
!     PDICT WILL BE ON SCRATCH2
 
 FILE   = scrt1
 nerror = 11
 CALL OPEN (*1760,scrt1,z(buf1),wrtrew)
 FILE = scrt2
 CALL OPEN (*1760,scrt2,z(buf2),wrtrew)
 
!     REQUESTED OUTPUT ELEMENT ENERGIES WILL BE TEMPORARILY WRITTEN ON
!     SCRT3 WHILE THE TOTAL ENERGY IS SUMMED.
 
 FILE   = scrt3
 IF (elnset /= 0) CALL OPEN (*1760,scrt3,z(buf3),wrtrew)
 nextgp = 1
 lastid = 0
 oldcod = 0
 toteng = 0.0D0
 estid  = 0
 axic   = .false.
 axif   = .false.
 
!     ECT PASS OF ALL ELEMENT TYPES PRESENT.
 
!     DETERMINE NEXT ELEMENT TYPE TO FIND ON ECT AND THEN FIND ITS
!     TYPE IN ECT.
 
 400 FILE = kdict
 nerror = 12
 CALL READ (*990,*1780,kdict,recid,3,noeor,iwords)
 kt = recid(1)
 
!           CCONAX       CTRIAAX       CTRAPAX
 IF (kt == 35 .OR. kt == 70 .OR. kt == 71) axic = .true.
!         CFLUID2/3/4  AND CFMASS
 IF (kt >= 43 .AND. kt <= 46) axif = .true.
!         CAXIF2/3/4 AND CSLOT3/4
 IF (kt >= 47 .AND. kt <= 51) axif = .true.
 
 FILE = ect
 CALL fwdrec (*1770,ect)
 410 CALL READ (*1770,*1780,ect,recidx,3,noeor,iwords)
!     2147483647 = 2**31-1
 IF (recidx(1) == 2147483647) GO TO 1770
 DO  i = 1,last,incr
   IF (elem(i+3) /= recidx(1)) CYCLE
   eltype = (i/incr) + 1
   ectwds = elem(i+5)
   IF (ectwds <= lbuf) GO TO 430
   WRITE  (outpt,420) swm,elem(i),elem(i+1)
   420 FORMAT (a27,' 2344. GPFDR FINDS ELEMENT = ',2A4,' HAS AN ECT ',  &
       'ENTRY LENGTH TOO LONG FOR A PROGRAM LOCAL ARRAY.')
   GO TO 1810
   
   430 grdpts= elem(i+ 9)
   grid1 = elem(i+12)
   name1 = elem(i   )
   name2 = elem(i+ 1)
   GO TO 470
 END DO
 
!     UNRECOGNIZED ELEMENT DATA ON ECT.
 
 WRITE  (outpt,450) swm,recidx
 450 FORMAT (a27,' 2345.  GPFDR FINDS AND IS IGNORING UNDEFINED ECT ',  &
     'DATA WITH LOCATE NUMBERS = ',3I8)
 FILE = ect
 
!     PASS THIS ECT RECORD BUT KEEP ESTID COUNTER IN SYNC.
 
 460 CALL READ (*1770,*410,ect,buf,ectwds,noeor,iwords)
 estid = estid + 1
 GO TO 460
 
 470 IF (eltype /= recid(1)) GO TO 460
 FILE  = kdict
 ldict = recid(2)
 IF (recid(3) == grdpts) GO TO 500
 480 WRITE  (outpt,490) swm,eltype,kdict
 490 FORMAT (a27,' 2346.  GPFDR FINDS DATA FOR EL-TYPE =',i9,  &
     ' IN DATA BLOCK',i9, /5X,  &
     'NOT TO BE IN AGREEMENT WITH THAT WHICH IS EXPECTED.')
 GO TO 1810
 
 500 ikdic  = nvec + 1
 nkdic  = nvec + ldict
 dicout = .false.
 engout = .false.
 
!     ALLOCATE A P-DICTIONARY FOR THE ELEMENTS GP-FORCE VECTOR
!     CONTRIBUTION.  CONTENTS = ESTID, EXT-EL.-ID, GINO-LOCS (GRDPTS)
 
 ipdic = nkdic + 1
 npdic = ipdic + grdpts + 1
 lpdic = grdpts + 2
 nerror= 13
 IF (npdic > core) GO TO 1800
 iloc1 = nkdic - grdpts
 phead(1) = eltype
 phead(2) = lpdic
 phead(3) = grdpts
 
!     LOOP IS NOW MADE ON THE ELEMENT ENTRIES OF THIS ELEMENT TYPE.
 
 nexten = 1
 
!     READ NEXT ELEMENT DICTIONARY FROM KDICT OF CURRENT ELEMENT TYPE
!     AND FIND ECT ENTRY WITH SAME ESTID.
 
 510 FILE = kdict
 CALL READ (*1770,*980,kdict,z(ikdic),ldict,noeor,iwords)
 FILE = ect
 nerror = 14
 520 CALL READ (*1770,*1780,ect,buf,ectwds,noeor,iwords)
 estid = estid + 1
 IF (z(ikdic)-estid < 0.0) THEN
   GO TO   480
 ELSE IF (z(ikdic)-estid == 0.0) THEN
   GO TO   530
 ELSE
   GO TO   520
 END IF
 
!     DECODE THE CODE WORD INTO A LIST OF INTEGERS
 
 530 IF (z(ikdic+3) == oldcod) GO TO 540
 oldcod = z(ikdic+3)
 CALL decode (oldcod,comps,ncomps)
 ncomp2 = ncomps
 IF (DOUBLE) ncomp2 = ncomps + ncomps
 
!     DETERMINE ACTIVE CONNECTIONS
 
 540 nsize  = z(ikdic+2)
 ngrids = nsize / ncomp2
 IF (ngrids <= grdpts) GO TO 560
 WRITE  (outpt,550) uwm,buf(1)
 550 FORMAT (a25,' 2347.  GPFDR FINDS TOO MANY ACTIVE CONNECTING GRID',  &
     ' POINTS FOR ELEMENT ID =',i9)
 GO TO 1810
 
!     ELEMENT ONLY DISPLACEMENT AND LOAD SPACE.
 
 560 iuge = npdic + 1
 IF (DOUBLE ) iuge = iuge/2 + 1
 nuge = iuge + nsize - 1
 ipge = nuge + 1
 npge = nuge + nsize
 IF (npge > core) GO TO 1800
 
!     ECT ENTRY AND K-DICTIONARY ENTRY NOW AT HAND.
 
!     SET FLAG IF EL-ENERGY IS TO BE OUTPUT FOR THIS ELEMENT.
 
 exelid    = buf(1)
 z(ipdic)  = estid
 z(ipdic+1)= exelid
 enflag    = .false.
 IF (axic) exelid = MOD(exelid,10000  )
 IF (axif) exelid = MOD(exelid,1000000)
 IF (elnset < 0.0) THEN
   GO TO   580
 ELSE IF (elnset == 0.0) THEN
   GO TO   590
 END IF
 
!     FIND THIS EXTERNAL ELEMENT ID IN THE REQUESTED SET LIST FOR
!     ELEMENT ENERGY OUTPUTS.
 
 570 CALL setfnd (*590,z(iellst),lellst,exelid,nexten)
 580 enflag = .true.
 590 gridl  = grid1 + grdpts - 1
 
!     REORDER ECT CONNECTION LIST ACCORDING TO SIL SEQUENCE.
 
 j = grid1 - 1
 600 j = j + 1
 IF (j >= gridl) GO TO 620
 gpsil = isilex + 2*buf(j) - 1
 lsil  = z(gpsil)
 i = j
 610 i = i + 1
 IF (i > gridl) GO TO 600
 gpsil = isilex + 2*buf(i) - 1
 isil  = z(gpsil)
 IF (isil > lsil) GO TO 610
 lsil   = buf(j)
 buf(j) = buf(i)
 buf(i) = lsil
 lsil   = isil
 GO TO 610
 
!     NOW SET INTERNAL GRID POINT ID-S IN THE ECT ENTRY NEGATIVE IF THEY
!     ARE TO HAVE THEIR GP-FORCE BALANCE OUTPUT.
 
 620 anygp = .false.
 IF (gpset == 0) GO TO 670
 DO  i = grid1,gridl
   IF (buf(i) > 0.0) THEN
     GO TO   630
   ELSE
     GO TO   660
   END IF
   630 IF (gpset  < 0.0) THEN
     GO TO   650
   ELSE IF (gpset  == 0.0) THEN
     GO TO   660
   END IF
   640 idx = isilex + 2*buf(i)
   id  = z(idx-2)/10
   IF (axic) id = MOD(id,1000000)
   IF (axif) id = MOD(id,500000 )
   IF (id < lastid) nextgp = 1
   lastid = id
   CALL setfnd (*660,z(igplst),lgplst,id,nextgp)
   650 buf(i) = -buf(i)
   anygp  = .true.
 END DO
 
!     IF NO GRID POINTS OF THIS ELEMENT WERE FLAGGED AND THERE IS
!     NO POTENTIAL OF ANY ELEMENT ENERGY OUTPUTS THEN SKIP THIS ELEMENT
!     AT THIS POINT.
 
 670 IF (.NOT.anygp .AND. elnset == 0) GO TO 510
 
!     BUILD A NON-EXPANDED ELEMENT DISPLACEMENT VECTOR AT THIS TIME.
 
 j = iuge
 DO  i = grid1,gridl
   IF (buf(i) < 0.0) THEN
     GO TO   680
   ELSE IF (buf(i) == 0.0) THEN
     GO TO   720
   ELSE
     GO TO   690
   END IF
   680 gpsil = isilex - 2*buf(i) - 1
   GO TO 700
   690 gpsil = isilex + 2*buf(i) - 1
   700 isil  = z(gpsil)
   DO  k = 1,ncomps
     lsil  = isil + comps(k)
     dz(j) = DBLE(rz(ivecz+lsil))
     j = j + 1
   END DO
 END DO
 
 IF (j-1 == nuge) GO TO 740
 WRITE  (outpt,730) swm,buf(1)
 730 FORMAT (a27,' 2348.  GPFDR DOES NOT UNDERSTAND THE MATRIX-',  &
     'DICTIONARY ENTRY FOR ELEMENT ID =',i9)
 GO TO 1810
 
!     TOTAL ELEMENT FORCE VECTOR IS NOW COMPUTED.
 
 740 DO  i = ipge,npge
   dz(i) = 0.0D0
 END DO
 
 jsize = nsize
 ikmat = npge  + 1
 IF (.NOT.DOUBLE) GO TO 760
 jsize = jsize + nsize
 ikmat = npge*2+ 1
 760 nkmat = ikmat + jsize - 1
 IF (nkmat > core) GO TO 1800
 diagm = .false.
 IF (z(ikdic+1) == 2) diagm = .true.
 
!     LOOP THROUGH ALL PARTITIONS ON KMAT FOR THIS ELEMENT.
 
 jpge = ipge
 DO  i = 1,grdpts
   itemp = iloc1 + i
   IF (z(itemp) == 0.0) THEN
     GO TO   870
   END IF
   770 CALL filpos (kmat,z(itemp))
   IF (diagm) GO TO 830
   
!     FULL MATRIX.  READ COLUMNS OF ROW-STORED VERETICAL PARTITION.
   
   nerror = 16
   DO  k = 1,ncomps
     CALL READ (*1770,*1780,kmat,z(ikmat),jsize,noeor,iwords)
     jkmat = ikmat
     IF (DOUBLE) GO TO 790
     DO  j = iuge,nuge
       dz(jpge) = dz(jpge) + dz(j)*DBLE(rz(jkmat))
       jkmat = jkmat + 1
     END DO
     GO TO 810
     
     790 DO  j = iuge,nuge
       iii(1) = z(jkmat  )
       iii(2) = z(jkmat+1)
       dz(jpge) = dz(jpge) + dz(j)*diii
       jkmat = jkmat + 2
     END DO
     
     810 jpge = jpge + 1
   END DO
   CYCLE
   
!     DIAGONAL MATRIX.  THUS ONLY DIAGONAL TERMS OF PARTITION CAN
!     BE READ.
   
   830 nerror = 17
   CALL READ (*1770,*1780,kmat,z(ikmat),ncomp2,noeor,iwords)
   IF (DOUBLE) GO TO 850
   
   DO  j = 1,ncomps
     dz(jpge) = dz(iuge+j-1)*DBLE(rz(ikmat+j-1))
     jpge = jpge + 1
   END DO
   CYCLE
   
   850 jkmat = ikmat
   DO  j = 1,ncomps
     iii(1) = z(jkmat)
     iii(2) = z(jkmat+1)
     dz(jpge) = dz(iuge+j-1)*diii
     jkmat = jkmat + 2
     jpge  = jpge + 1
   END DO
   
 END DO
 
!     ENERGY COMPUTATION IS NOW MADE IF NECESSARY.
 
!       U   =  0.5(PG ) X (UG )
!        T           E         E
 
 
 IF (elnset == 0.0) THEN
   GO TO   900
 END IF
 880 jpge  = ipge
 elengy= 0.0D0
 DO  i = iuge,nuge
   elengy= elengy + dz(i)*dz(jpge)
   jpge  = jpge + 1
 END DO
 
!     NOTE, TOTAL ENERGY WILL BE DIVIDED BY 2.0 LATER.
 
 toteng = toteng + elengy
 
!     WRITE THIS ELEMENTS ENERGY ON SCRT3 FOR LATER OUTPUT IF REQUESTED.
 
 IF (.NOT. enflag) GO TO 900
 out (1) = buf(1)
 rout(2) = SNGL(elengy)*0.50
 IF (.NOT.engout) CALL WRITE (scrt3,names,2,noeor)
 CALL WRITE (scrt3,out,2,noeor)
 engout  = .true.
 
!     GRID POINT FORCE BALANCE OUTPUTS FOR REQUESTED GIRD POINTS.
 
 900 IF (.NOT. anygp) GO TO 970
 
!     EXPAND TO 6X1 FROM PGE EACH GRID POINT FORCE TO BE OUTPUT.
 
!     FORCES COMPUTED FOR COMPONENTS OTHER THAN 1 THRU 6 ARE NOT
!     NOW OUTPUT FROM MODULE GPFDR...  FUTURE ADDITIONAL CAPABLILITY.
!     OFP MODS NEEDED AT THAT TIME.
 
 jpge   = ipge
 dicloc = ipdic + 2
 DO  i = dicloc,npdic
   z(i) = 0
 END DO
 DO  i = grid1,gridl
   IF (buf(i) < 0.0) THEN
     GO TO   930
   ELSE IF (buf(i) == 0.0) THEN
     GO TO   960
   END IF
   
!     THIS GRID POINT NOT IN GP-FORCE BALANCE REQUEST LIST.
   
   920 jpge = jpge + ncomps
   dicloc = dicloc + 1
   CYCLE
   
!     OK THIS GRID POINT GETS OUTPUT.
   
   930 DO   j = 1,6
     vec(j) = 0.0
   END DO
   DO  j = 1,ncomps
     comp = comps(j)
     IF (comp <= 5) vec(comp+1) =-SNGL(dz(jpge))
     jpge = jpge + 1
   END DO
   
   CALL WRITE  (scrt1,vec,6,eor)
   CALL savpos (scrt1,z(dicloc))
   dicloc = dicloc + 1
 END DO
 
!     OUTPUT THE DICTIONARY
 
 IF (.NOT.dicout) CALL WRITE (scrt2,phead,3,noeor)
 CALL WRITE (scrt2,z(ipdic),lpdic,noeor)
 dicout = .true.
 
!     GO FOR NEXT ELEMENT OF CURRENT TYPE.
 
 970 GO TO 510
 
!     END OF ELEMENT ENTRIES OF CURRENT ELEMENT TYPE.
!     COMPLETE RECORDS IN PDIC, AND SCRT3=EL-ENERGY.
 
 980 IF (dicout) CALL WRITE (scrt2,0,0,eor)
 IF (engout) CALL WRITE (scrt3,0,0,eor)
 
!     GO FOR NEXT ELEMENT TYPE
 
 GO TO 400
 
!     END OF ALL ELEMENT DATA ON ECT (WRAP UP PHASE I OF GPFDR).
 
 990 CALL CLOSE (kmat,clsrew)
 CALL CLOSE (kdict,clsrew)
 CALL CLOSE (ect,clsrew)
 CALL CLOSE (scrt1,clsrew)
 CALL CLOSE (scrt2,clsrew)
 CALL CLOSE (scrt3,clsrew)
 
!     PREPARE AND WRITE THE ELEMENT ENERGY OUTPUTS NOW RESIDENT ON SCRT3
 
 IF (elnset == 0) GO TO 1050
 
!     OFP ID RECORD DATA
!     DEVICE, OFP-TYPE, TOTAL ENERGY, SUBCASE, ELEMENT NAME, WORDS
!     PER ENTRY.
 
 idrec( 1) = 10*branch + eldvis
 idrec( 2) = 18
 ridrec(3) = SNGL(toteng)*0.50
 idrec( 4) = subcas
 idrec(10) = 3
 
!     IF APPROACH IS REIG, PUT MODE NO. AND FREQ. INTO IDREC, 8 AND 9
!     WORDS
 
 IF (branch /= 2) GO TO 1000
 ridrec(9) = rz(lfeq-mode)
 mode      = mode + 1
 idrec( 8) = mode
 
 1000 nerror = 22
 FILE = onrgy1
 CALL OPEN (*1760,onrgy1,z(buf2),wrt)
 FILE = scrt3
 CALL OPEN (*1760,scrt3,z(buf3),rdrew)
 
!     TOTENG FACTOR FOR MULTIPLICATION TO GET DECIMAL PERCENTAGE BELOW
 
 IF (toteng /= 0.0D0) toteng = 200.0D0/toteng
 
!     READ ELEMENT NAME INTO IDREC RECORD.
 
 jtype = 0
 1010 CALL READ  (*1040,*1780,scrt3,idrec(6),2,noeor,iwords)
 CALL WRITE (onrgy1,idrec,146,eor)
 1020 CALL READ  (*1770,*1030,scrt3,buf,2,noeor,iwords)
 jtype  = jtype + 1
 buf(1) = 10*buf(1) + eldvis
 rbuf(3)= rbuf(2)*SNGL(toteng)
 CALL WRITE (onrgy1,buf,3,noeor)
 GO TO 1020
 
 1030 CALL WRITE (onrgy1,0,0,eor)
 GO TO 1010
 
 1040 CALL CLOSE (onrgy1,clseof)
 mcb(1) = onrgy1
 CALL rdtrl (mcb)
 mcb(2) = mcb(2) + jtype
 CALL wrttrl (mcb)
 CALL CLOSE (scrt3,clsrew)
 idrec(3) = 0
 idrec(6) = 0
 idrec(7) = 0
 
!     A GRID-POINT-FORCE-BALANCE-OUTPUT-MAP IS NOW CONSTRUCTED. (GPFBOM)
 
!     CONTENTS...  1 LOGICAL RECORD FOR EACH GRID POINT TO BE OUTPUT
!     ===========
 
!     REPEATING 4       * EXTERNAL-ELEMENT-ID
!     WORD ENTRIES     *  ELEMENT NAME FIRST 4H
!     OF THE CON-      *  ELEMENT NAME LAST  4H
!     NECTED ELEMENTS   * GINO-LOC OF THE 6X1 FORCE VECTOR CONTRIBUTION
 
!     FOR EACH RECORD WRITTEN ABOVE, A 3-WORD ENTRY IS WRITTEN TO A
!     COMPANION DICTIONARY FILE GIVING,
 
!                      *  1-THE EXTERNAL GRID POINT ID
!     REPEATING ENTRY *   2-THE GINO-LOC TO THE ABOVE RECORD
!                      *  3-THE NUMBER OF ENTRIES IN THE RECORD
 
 
!     ALLOCATE A TABLE WITH AN ENTRY FOR EACH ELEMENT TYPE.
!     POSSIBLE IN IT.  EACH ENTRY TO HAVE 3 WORDS.
 
!     ENTRY I =      1= PTR TO DICTIONARY DATA FOR ELEMENT TYPE-I
!     *********      2= LENGTH OF DICTIONARY DATA
!                    3= NUMBER OF ENTRIES
 
 1050 IF (gpset == 0) GO TO 110
 idtab = ncc + 1
 ndtab = idtab + nelems*3 - 1
 jdicts= ndtab + 1
 IF (jdicts > core) GO TO 1800
 DO  i = idtab,ndtab
   z(i) = 0
 END DO
 
!     READ IN DICTIONARIES OF PMAT VECTORS.  (SCRT2)
 
 FILE = scrt2
 CALL OPEN (*1760,scrt2,z(buf2),rdrew)
 
!     READ AN ELEMENT TYPE HEADER (FIRST 3-WORDS OF EACH RECORD)
 
 1070 CALL READ (*1090,*1780,scrt2,buf,3,noeor,iwords)
 itype = buf(1)
 ldict = buf(2)
 grdpts= buf(3)
 k     = incr*itype - incr
 j     = idtab + 3*itype - 3
 z(j)  = jdicts
 
!     BLAST READ IN THE DICTIONARIES OF THIS TYPE.
 
 CALL READ (*1770,*1080,scrt2,z(jdicts),core-jdicts,noeor,iwords)
 nerror = 18
 GO TO 1800
 
 1080 z(j+1) = iwords
 z(j+2) = iwords/ldict
 jdicts = jdicts + iwords
 nerror = 19
 IF (core-jdicts > 0.0) THEN
   GO TO  1070
 ELSE
   GO TO  1800
 END IF
 
 1090 CALL CLOSE (scrt2,clsrew)
 
!     DICTIONARIES ALL IN CORE.  SCRT2 IS AVAILABLE FOR USE AS THE
!     -GPFBOM-.
 
 ndicts = jdicts - 1
 
!     PASS THE -GPECT- AND BUILD THE -GPFBOM- (ON SCRT2) AND ITS
!     COMPANION DICTIONARY FILE (ON SCRT3).
 
 FILE = scrt2
 CALL OPEN (*1760,scrt2,z(buf2),wrtrew)
 FILE = scrt3
 CALL OPEN (*1760,scrt3,z(buf3),wrtrew)
 
 FILE = gpect
 CALL OPEN (*1760,gpect,z(buf4),rdrew)
 CALL fwdrec (*1770,gpect)
 oldid = 0
 next  = 1
 
!     READ PIVOT HEADER DATA FROM -GPECT- RECORD.
 
 1100 CALL READ (*1240,*1780,gpect,buf,2,noeor,iwords)
 pivot = buf(1)
 
!     CONVERT SIL TO EX-ID
 
 CALL bisloc (*1200,pivot,z(isilex+1),2,points,j)
 j = isilex + j - 1
 extid = z(j)/10
 idext = extid
 IF (axic) idext = MOD(extid,1000000)
 IF (axif) idext = MOD(extid,500000 )
 nentry = 0
 
!     CHECK FOR OUTPUT REQUEST THIS EX-ID
 
 IF (gpset < 0.0) THEN
   GO TO  1120
 ELSE IF (gpset == 0.0) THEN
   GO TO  1220
 END IF
 1110 IF (idext < oldid) next = 1
 oldid = idext
 CALL setfnd (*1220,z(igplst),lgplst,idext,next)
 
!     YES GP-FORCE BALANCE FOR PIVOT IS TO BE OUTPUT.
 
!     PROCESS ALL ELEMENTS CONNECTING THIS PIVOT.
 
 1120 CALL READ (*1770,*1230,gpect,length,1,noeor,iwords)
 length = IABS(length)
 IF (length <= lbuf) GO TO 1140
 WRITE  (outpt,1130) swm,pivot,gpect
 1130 FORMAT (a27,' 2349.  GPFDR FINDS AN ELEMENT ENTRY CONNECTING ',  &
     'PIVOT SIL =',i9,' ON DATA BLOCK',i5, /5X,  &
     'TOO LARGE FOR A LOCAL ARRAY. ENTRY IS BEING IGNORED.')
 CALL READ (*1770,*1780,gpect,0,-length,noeor,iwords)
 GO TO 1120
 
!     LOCATE ELEMENT FORCE DICTIONARY FOR THIS ELEMENT ENTRY.
 
 1140 CALL READ (*1770,*1780,gpect,buf,length,noeor,iwords)
 ktype = buf(2)*3 - 3 + idtab
 ptr   = z(ktype)
 ldicts= z(ktype+2)
 IF (ldicts == 0) GO TO 1180
 n = z(ktype+1)
 CALL bisloc (*1180,buf(1),z(ptr),n/ldicts,ldicts,j)
 j = ptr + j
 out(1) = z(j)
 
!     FOUND DICTIONARY.  DETERMINE GINO-LOC TO USE.
 
 DO  i = 3,length
   j = j + 1
   IF (buf(i) == pivot .AND. z(j) > 0) GO TO 1170
 END DO
 WRITE  (outpt,1160) swm,pivot,out(1),gpect
 1160 FORMAT (a27,' 2350.  GPFDR CANNOT FIND PIVOT SIL =',i10, /5X,  &
     'AMONG THE SILS OF ELEMENT ID =',i9,  &
     ' AS READ FROM DATA BLOCK',i5,',  ENTRY THUS IGNORED.')
 GO TO 1120
 
 1170 k = buf(2)*incr - incr
 out(2) = elem(k+1)
 out(3) = elem(k+2)
 out(4) = z(j)
 
!     GINO-LOC IN P-DICTIONARY NO LONGER NEEDED, THUS SET IT NEGATIVE
!     TO AVOID RE-USE IN CASE WHERE AN ELEMENT CONNECTS SAME GRID MORE
!     THAN ONCE.
 
 z(j) = -z(j)
 
!     OUTPUT THE 4-WORD ENTRY TO -GPFBOM-
 
 CALL WRITE (scrt2,out,4,noeor)
 
!     INCREMENT COUNTS
 
 nentry = nentry + 1
 
!     GET THE NEXT ELEMENT ENTRY.
 
 GO TO 1120
 
!     HERE WHEN PMAT DICTIONARY MISSING FOR AN ELEMENT
!     CONNECTED TO A GRID POINT TO HAVE GP-FORCE BALANCE OUTPUT.
 
 1180 kkk = buf(2)*incr - incr
 WRITE  (outpt,1190) uim,elem(kkk+1),elem(kkk+2),extid
 1190 FORMAT (a29,' 2351. A FORCE CONTRIBUTION  DUE TO ELEMENT TYPE = ',  &
     2A4,', ON POINT ID =',i10, /5X,  &
     'WILL NOT APPEAR IN THE GRID-POINT-FORCE-BALANCE SUMMARY.')
 GO TO 1120
 
!     SIL NOT FOUND IN LIST OF SILS, OR NOT REQUESTED.
 
 1200 WRITE  (outpt,1210) swm,pivot,gpect
 1210 FORMAT (a27,' 2352.  GPFDR IS NOT ABLE TO FIND PIVOT SIL =',i10,  &
     ' AS READ FROM DATA BLOCK',i5, /5X,'IN TABLE OF SILS.')
 
 1220 CALL fwdrec (*1770,gpect)
 GO TO 1100
 
!     HERE WHEN END OF RECORD ON GPECT.
!     COMPLETE THE RECORD ON -GPFBOM- AND WRITE DICTIONARY ENTRY FOR THE
!     COMPLETED RECORD.
 
 1230 CALL WRITE (scrt2,0,0,eor)
 buf(1) = extid
 buf(3) = nentry
 CALL savpos (scrt2,buf(2))
 CALL WRITE (scrt3,buf,3,noeor)
 
!     GO FOR NEXT PIVOT SIL
 
 GO TO 1100
 
!     HERE WHEN END OF FILE ON -GPECT-.
 
 1240 CALL CLOSE (gpect,clsrew)
 CALL CLOSE (scrt2,clsrew)
 CALL CLOSE (scrt3,clsrew)
 
!     SO AS TO OUTPUT THE FORCE BALANCES IN EXTERNAL GRID POINT ORDER
!     THE FOLLOWING STEPS ARE NOW PERFORMED ON THE DICTIONARY ENTRIES OF
!     THE -GPFBOM- COMPANION FILE (SCRT3).
 
!     1) ALL OF THE COMPANION FILE DICTIONARIES ARE READ INTO CORE.
!     2) THEY ARE SORTED ON THE EXTERNAL IDS.
!     3) THEY ARE PARTITIONED INTO GROUPS BASED ON A CONSIDERATION OF
!        THE NEED FOR 12 WORDS OF CORE FOR EACH ENTRY OF EACH -GPFBOM-
!        RECORD REPRESENTED BY THE GROUP IN THE FINAL OUTPUT PASS.
!     4) EACH ENTRYS 3-RD WORD (THE NUMBER OF ENTRIES IN THE RECORD) IS
!        REPLACED WITH THE INTEGER POSITION OF THE ENTRY IN THE GROUP.
!     5) EACH GROUP IS SORTED ON GINO-LOC AND WRITTEN BACK
!        TO THE COMPANION FILE AS A LOGICAL RECORD.  (THIS INSURES THAT
!        NO MORE THAN ONE PASS OF THE -GPFBOM- IS MADE PER GROUP WHEN
!        CONSTRUCTING TABLE-1 AND TABLE-2 IN THE FINAL OUTPUT PASS.)
 
 FILE = scrt3
 nerror = 20
 CALL OPEN (*1760,scrt3,z(buf3),rdrew)
 
!     BLAST-READ 3-WORD -GPFBOM- DICTIONARY ENTRIES INTO CORE.
 
 idicts = ncc + 1
 CALL READ (*1770,*1250,scrt3,z(idicts),core-idicts,noeor,iwords)
 GO TO 1800
 
 1250 ndicts = idicts + iwords - 1
 CALL CLOSE (scrt3,clsrew)
 nerror = 21
 CALL OPEN (*1760,scrt3,z(buf3),wrtrew)
 
!     SORT ENTRIES ON EXTERNAL ID
 
 CALL sort (0,0,3,1,z(idicts),iwords)
 
!     DETERMINE A -GPFBOM- GROUP OF RECORDS FOR OUTPUT.  EACH -GPFBOM-
!     RECORDS ENTRY WILL REQUIRE 12 WORDS OF CORE IN THE FINAL OUTPUT
!     PROCEEDURES.
 
 entrys = (core-ncc)/12
 1260 j = idicts
 total = 0
 1270 IF (total+z(j+2) > entrys) GO TO 1280
 total = total + z(j+2)
 j = j + 3
 IF (j < ndicts) GO TO 1270
 
!     GROUP RANGE HAS BEEN FOUND.  REPLACE EACH ENTRYS -GPFBOM- ENTRY
!     COUNT WITH THE OUTPUT ORDER OF THE EXTERNAL ID ENTRY HERE.
 
 1280 jdicts = j - 1
 k  = 1
 DO  i = idicts,jdicts,3
   jk = z(i+2)
   z(i+2) = k
   k  = k + jk
 END DO
 
!     SORT THIS GROUP OF 3-WORD ENTRIES ON THE GINO-LOCS.
 
 length = jdicts - idicts + 1
 CALL sort (0,0,3,2,z(idicts),length)
 
!     OUTPUT AS A LOGICAL RECORD.
 
 CALL WRITE (scrt3,z(idicts),length,eor)
 
!     PROCESS NEXT GROUP IF THERE ARE MORE.
 
 idicts = jdicts + 1
 IF (idicts < ndicts) GO TO 1260
 
!     ALL GROUPS HAVE BEEN DETERMINED, SEQUENCED, SORTED ON GINO-LOCS,
!     AND OUTPUT.
 
 CALL CLOSE (scrt3,clsrew)
 
!     PREPARE GRID-POINT-FORCE-BALANCE ENTRIES WITH RESPECT TO APPLIED-
!     LOAD AND SINGLE-POINT-CONSTRAINT FORCES.
 
!     LINE ENTRIES WILL BE WRITTEN TO SCRT4 FROM THE VECTOR IN CORE
!     FOR EACH OF PG AND QG CONTAINING,
 
!     EXTERNAL GP ID, 0, 4H----, 4H----, T1, T2, T3, R1, R2, R3,
 
!     ONLY FOR THOSE POINTS WHICH MAY BE OUTPUT IN THE GRID-POINT FORCE
!     BALANCE.
 
!     (NULL ENTRIES ARE NOT OUTPUT)
 
!     AFTER ALL ENTRIES FOR PG AND QG DESIRED HAVE BEEN WRITTEN TO
!     SCRT4 THEY ARE BROUGHT BACK INTO CORE, SORTED ON EXTERNAL GP ID
!     AND RE-OUTPUT TO SCRT4.
 
 FILE = scrt4
 CALL OPEN (*1760,scrt4,z(buf1),wrtrew)
 
!     PROCESS PG.
 
 ugpgqg = pg
 buf(2) = 0
 buf(3) = oload(1)
 buf(4) = oload(2)
 lastid = 0
 nextgp = 1
 ASSIGN 1300 TO icont
 GO TO 1400
 
!     PROCESS QG
 
 1300 ugpgqg = qg
 buf(3) = ospcf(1)
 buf(4) = ospcf(2)
 lastid = 0
 nextgp = 1
 ASSIGN 1310 TO icont
 GO TO 1400
 
!     SORT SCRT4 ENTRIES ON EXTERNAL GP ID
 
 1310 CALL WRITE (scrt4,0,0,eor)
 CALL CLOSE (scrt4,clsrew)
 movepq = 0
 CALL OPEN (*1760,scrt4,z(buf1),rdrew)
 CALL READ (*1770,*1330,scrt4,z(icc),buf1-icc,noeor,iwords)
 WRITE  (outpt,1320) uwm,subcas
 1320 FORMAT (a25,' 2353.  INSUFFICIENT CORE TO HOLD ALL NON-ZERO APP-',  &
     'LOAD AND F-OF-SPC OUTPUT LINE ENTRIES OF', /5X,  &
     'GRID-POINT-FORCE-BALANCE REQUESTS. SOME POINTS REQUESTED',  &
     ' FOR OUTPUT WILL BE MISSING THEIR APP-LOAD OR F-OF-SPC',  &
     /5X,'CONTRIBUTION IN THE PRINTED BALANCE.')
 iwords = buf1 - icc - MOD(buf1-icc,10)
 1330 CALL sort  (0,0,10,1,z(icc),iwords)
 CALL CLOSE (scrt4,clsrew)
 CALL OPEN  (*1760,scrt4,z(buf1),wrtrew)
 CALL WRITE (scrt4,z(icc),iwords,eor)
 CALL CLOSE (scrt4,clsrew)
 GO TO 1560
 
!     INTERNAL ROUTINE TO GET A VECTOR IN CORE (PG OR QG) AND WRITE
!     SELECTED NON-ZERO ENTRIES TO SCRT4 FOR INCLUSION LATER IN THE
!     GRID-POINT-FORCE-BALANCE.
 
 1400 CALL OPEN (*1550,ugpgqg,z(buf2),rd)
 IF (movepq < 0) THEN
   GO TO  1410
 ELSE IF (movepq == 0) THEN
   GO TO  1450
 ELSE
   GO TO  1430
 END IF
 
!     BACK POSITION DATA BLOCK
 
 1410 j = IABS(movepq)
 DO  i = 1,j
   CALL bckrec (ugpgqg)
 END DO
 GO TO 1450
 
!     FORWARD POSITION DATA BLOCK
 
 1430 FILE = ugpgqg
 DO  i = 1,movepq
   CALL fwdrec (*1770,ugpgqg)
 END DO
 
!     GET VECTOR INTO CORE.
 
 1450 ASSIGN 1460 TO iretrn
 GO TO 220
 
!     OUTPUT NON-ZERO ENTRIES REQUESTED
 
 1460 CALL CLOSE (ugpgqg,cls)
 DO  i = isilex,nsilex,2
   icode = MOD(z(i),10)
   i1 = ivecz + z(i+1)
   i2 = i1 + scale(icode)
   DO  j = i1,i2
     IF (rz(j) == 0.0) THEN
       GO TO  1470
     ELSE
       GO TO  1480
     END IF
   END DO
   CYCLE
   
!     NON-ZERO ENTRY.  CHECK FOR OUTPUT.
   
   1480 buf(1) = z(i)/10
   ibuf1  = buf(1)
   IF (axic) ibuf1 = MOD(ibuf1,1000000)
   IF (axif) ibuf1 = MOD(ibuf1,500000 )
   IF (ibuf1 < lastid) nextgp = 1
   lastid = ibuf1
   IF (gpset < 0.0) THEN
     GO TO  1500
   ELSE IF (gpset == 0.0) THEN
     GO TO  1550
   END IF
   1490 CALL setfnd (*1540,z(igplst),lgplst,ibuf1,nextgp)
   1500 l = 5
   DO  j = i1,i2
     buf(l) = z(j)
     l = l + 1
   END DO
   IF (l >= 11) GO TO 1530
   DO  j = l,10
     rbuf(l) = 0.0
   END DO
   1530 buf(1) = buf(1)*10 + gpdvis
   CALL WRITE (scrt4,buf,10,noeor)
 END DO
 1550 GO TO icont, (1300,1310)
 1560 CONTINUE
 
!     FINAL OUTPUT PHASE FOR CURRENT CASE CONTROL.
 
!     THE -GPFBOM- COMPANION FILE IS PROCESSED RECORD BY RECORD.
 
!     FOR EACH RECORD THEN,
 
!     1) A 3-WORD ENTRY IS READ GIVING 1) EXTERNAL GP-ID
!                                      2) GINO-LOC OF -GPFBOM- RECORD
!                                      3) OUTPUT ORDER WITHIN THE GROUP.
 
!     2) -GPFBOM- IS POSITIONED USING THE GINO-LOC.
 
!     3) A POINTER IS DETERMINED INTO TABLE-2 OF WHERE OUTPUTS BELONG
!        =10*ORDER - 10  (A ZERO POINTER) + TABLE BASE (A ZERO POINTER)
 
!     4) ENTRIES ARE READ FROM -GPFBOM- CONTAINING,
 
!                                      1) EXTERNAL ELEMENT ID
!                                      2) ELEMENT NAME FIRST 4H
!                                      3) ELEMENT NAME LAST  4H
!                                      4) GINO LOC TO 6X1 FORCE VECTOR
 
!        UNTIL AN EOR IS ENCOUNTERED.
 
!        FOR EACH ENTRY READ A 2-WORD ENTRY IS ADDED TO TABLE-1
!        CONSISTING OF                 1) GINO-LOC TO THE 6X1 VECTOR
!                                      2) PTR INTO TABLE-2
 
!        AND A 10-WORD ENTRY IS ADDED TO TABLE-2 AT Z(PTR)
!        CONSISTING OF                 1) EXTERNAL GP-ID
!                                      2) EXTERNAL ELEMENT-ID
!                                      3) NAME FIRST 4H
!                                      4) NAME LAST  4H
!                                      5 THRU 10)   NOT SET YET.
 
!     5) WHEN ALL ENTRIES OF THE -GPFBOM- RECORDS OF THE GROUP
!        (AS SPECIFIED BY ONE RECORD ON THE COMPANINON FILE) ARE IN CORE
!        TABLE-1 IS SORTED ON GINO LOCS.
!        THIS WILL PREVENT HAVING TO MAKE MORE THAN ONE PASS
!        OF THE PMAT DATA PER GROUP.
 
!     6) A SERIAL PASS OF TABLE-1 IS MADE AND EACH 6X1 VECTOR IS
!        READ DIRECTLY INTO Z(PTR+4) OF TABLE-2.
 
!     7) OUTPUT TO THE FORCE BALANCE DATA BLOCK IS MADE WITH THE
!        STANDARD OFP METHOD OF HEADER RECORD, AND REPEATING ENTRY DATA
!        RECORD.  A HEADER RECORD WILL BE OUTPUT EACH TIME THE GRID
!        POINT CHANGES.
 
 
!     ALLOCATE TABLE-1 AND TABLE-2
 
 itab1 = ncc + 1
 ntab1 = ncc + 2*entrys
 itab2 = ntab1 + 1
 
!     OPEN -GPFBOM- (SCRT2) AND ITS COMPANION DICTIONARY FILE (SCRT3).
 
 FILE   = scrt2
 nerror = 23
 CALL OPEN (*1760,scrt2,z(buf2),rdrew)
 FILE = scrt3
 CALL OPEN (*1760,scrt3,z(buf3),rdrew)
 
!     OPEN THE OUTPUT FILE FOR GP-FORCES.
 
 FILE = ogpf1
 CALL OPEN (*1760,ogpf1,z(buf4),wrt)
 lines    = 0
 idrec(1) = 10*branch + gpdvis
 idrec(2) = 19
 idrec(4) = subcas
 idrec(10)= 10
 
!     OPEN THE PMAT 6X1 FORCE VECTORS FILE.
 
 FILE = scrt1
 CALL OPEN (*1760,scrt1,z(buf1),rdrew)
 
!     INITIALIZE INPUT OF APP-LOAD AND F-OF-SPC LINE ENTRIES FROM SCRT4.
 
 FILE = scrt4
 CALL OPEN (*1760,scrt4,z(buf5),rdrew)
 CALL READ (*1770,*1570,scrt4,kvec,10,noeor,iwords)
 eorst4 = .false.
 GO TO 1580
 1570 eorst4 = .true.
 1580 CONTINUE
 
!     PROCESS ONE GROUP OF -GPFBOM- RECORDS AS SPECIFIED BY THE 3-WORD
!     ENTRIES OF ONE RECORD ON SCRT3.
 
 any   = .false.
 oldid = 0
 CALL WRITE (ogpf1,idrec,146,eor)
 1590 iptr1 = itab1 - 1
 jtab1 = itab1 - 1
 jtab2 = itab2 - 1
 FILE  = scrt2
 1600 CALL READ (*1740,*1620,scrt3,buf,3,noeor,iwords)
 extgp = buf(1)
 loc   = buf(2)
 iptr2 = itab2 + 10*buf(3) - 11
 
!     POSITION -GPFBOM- TO RECORD OF 4-WORD ENTRIES FOR THIS EXTERNAL GP
 
 CALL filpos (scrt2,loc)
 nerror = 24
 
!     READ AND DISTRIBUTE THE DATA OF THE 4-WORD ENTRIES.
 
 1610 CALL READ (*1770,*1600,scrt2,buf,4,noeor,iwords)
 z(iptr1+1) = buf(4)
 z(iptr1+2) = iptr2
 z(iptr2+1) = extgp
 z(iptr2+2) = buf(1)
 z(iptr2+3) = buf(2)
 z(iptr2+4) = buf(3)
 iptr1 = iptr1 + 2
 iptr2 = iptr2 + 10
 jtab1 = jtab1 + 2
 jtab2 = jtab2 + 10
 GO TO 1610
 
!     HERE ON END OF A GROUP.  SORT TABLE-1 ON GINO LOCS.
!     AND FILL TABLE-2 WITH 6X1 FORCE VECTORS.
 
 1620 CALL sort (0,0,2,1,z(itab1),jtab1-itab1+1)
 
 nerror= 25
 FILE  = scrt1
 DO  i = itab1,jtab1,2
   CALL filpos (scrt1,z(i))
   ptr = z(i+1)
   CALL READ (*1770,*1780,scrt1,z(ptr+5),6,noeor,iwords)
 END DO
 
!     OUTPUT DATA.  START NEW SUM WHEN ENCOUNTERING A NEW GP-ID.
!     APPLIED-LOADS AND FORCES-OF-SPC WILL INITIALIZE SUM, IF THEY EXIST
!     FOR GRID POINT IN QUESTION,  OHTERWISE SUM IS INITIALIZED TO ZERO.
 
 DO  i = itab2,jtab2,10
   
!     IS THIS SAME GRID POINT ID AS CURRENTLY BEING SUMMED.  IF SO,
!     CONTINUE OUTPUT OF LINE ENTRY AND SUM IN.  OTHERWISE OUTPUT
!     SUM LINE, AND NEW ID-S APPLIED-LOAD AND F-OF-SPC ENTRY.
   
   1640 IF (z(i) == oldid) GO TO 1710
   
!     CHANGE IN GRID POINT ID.
   
   isum(1) = oldid*10 + gpdvis
   IF (any) CALL WRITE (ogpf1,isum,10,noeor)
   IF (any) lines = lines + 1
   any = .false.
   
!     OUTPUT ALL LINE ENTRIES OF APP-LOADS AND F-OF-SPC UNTIL
!     MATCH ON NEW ID IS FOUND OR CURRENT FVEC IS NOT YET NEEDED.
   
   IF (eorst4) GO TO 1690
   IF (kvec(1)/10 > z(i)) GO TO 1690
   DO  j = 5,10
     rsum(j) = fvec(j)
   END DO
   oldid = kvec(1)/10
   CALL WRITE (ogpf1,kvec,10,noeor)
   lines = lines + 1
   any   = .true.
   
!     SUM IN ANY MORE FROM SCRT4 OF CURRENT ID, OUTPUT LINE ENTRIES.
   
   1660 CALL READ (*1770,*1680,scrt4,kvec,10,noeor,iwords)
   IF (kvec(1)/10 /= oldid) GO TO 1640
   CALL WRITE (ogpf1,kvec,10,noeor)
   lines = lines + 1
   DO  j = 5,10
     rsum(j) = rsum(j) + fvec(j)
   END DO
   GO TO 1660
   
   1680 eorst4 = .true.
   GO TO 1640
   
!     NO APP-LOAD OR F-OF-SPC ENTRIES LEFT OR CURRENT ONE NOT NEEDED YET
   
   1690 DO  j = 5,10
     rsum(j) = 0.0
   END DO
   any  = .true.
   oldid= z(i)
   
   1710 z(i) = 10*z(i) + gpdvis
   CALL WRITE (ogpf1,z(i),10,noeor)
   lines = lines + 1
   DO  j = 5,10
     rsum(j) = rsum(j) + rz(i+j-1)
   END DO
   
 END DO
 
 isum(1) = oldid*10 + gpdvis
 IF (any) CALL WRITE (ogpf1,isum,10,noeor)
 IF (any) lines = lines + 1
 any = .false.
 
!     GO FOR NEXT GROUP FROM THE -GPFBOM-.
 
 GO TO 1590
 
!     HERE ON EOF ON -GPFBOM- COMPANION FILE.  THUS AT CONCLUSION OF
!     OUTPUT PHASE FOR GP-FORCE BALANCE ONE SUBCASE, OR ONE TIME STEP OF
!     ONE SUBCASE.
 
 1740 CALL CLOSE (scrt1,clsrew)
 CALL CLOSE (scrt2,clsrew)
 CALL CLOSE (scrt3,clsrew)
 CALL CLOSE (scrt4,clsrew)
 mcb(1) = ogpf1
 CALL rdtrl (mcb)
 mcb(2) = mcb(2) + lines
 CALL wrttrl (mcb)
 CALL CLOSE (ogpf1,clseof)
 GO TO 110
 
!     NORMAL COMPLETION.
 
 1750 CALL CLOSE (casecc,clsrew)
 CALL CLOSE (ug,clsrew)
 RETURN
 
!     HERE ON ERROR CONDITIONS.
 
 1760 mm = 1
 GO TO 1790
 1770 mm = 2
 GO TO 1790
 1780 mm = 3
 1790 CALL mesage (mm,FILE,subr)
 GO TO 1810
 1800 CALL mesage (8,0,subr)
 1810 WRITE  (outpt,1820) swm,nerror
 1820 FORMAT (a27,' 2354.' ,/5X,'GPFDR MODULE IS UNABLE TO CONTINUE ',  &
     'AND HAS BEEN TERMINATED DUE TO ERROR MESSAGE PRINTED ',  &
     'ABOVE OR BELOW THIS MESSAGE.', /5X,'THIS ERROR OCCURRED ',  &
     'IN GPFDR CODE WHERE THE VARIABLE -NERROR- WAS SET =',i5)
 DO  i = 100,300,100
   DO  j = 1,9
     CALL CLOSE (i+j,clsrew)
   END DO
 END DO
 RETURN
END SUBROUTINE gpfdr
