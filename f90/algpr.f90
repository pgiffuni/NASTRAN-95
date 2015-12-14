SUBROUTINE algpr (ierr)
     
 
 INTEGER, INTENT(OUT)                     :: ierr
 LOGICAL :: debug
 INTEGER :: sysbuf,NAME(2),edt,eqexin,cstm,ugv,FILE,corwds,  &
     pgeom,buf1,buf2,scr1,scr2,ret2,typout,bgpdt,  &
     itrl(7),stream(3),apress,atemp,strml,algdb,  &
     idata(24),kptsa(10),ifangs(10),rd,rdrew,wrt,  &
     wrtrew,clsrew,norew,LEN(3),ifill(3),algdd
 REAL :: rfill(3),z(1),ta(9),rdata(6),xsta(21,10),  &
     rsta(21,10),r(21,10),b1(21),b2(21),rle(21),  &
     tc(21),te(21),cord(21),delx(21),dely(21),zed(21),  &
     phi(2,21),zr(21),pp(21),qq(21),cord2(21),  &
     fchord(21),jz(21),xb(21,10),yb(21,10),zb(21,10),  &
     dispt(3),dispr(3),dispt1(21,10),dispt2(21,10),  &
     dispt3(21,10),blafor(21,10),dispr1(21,10), dispr2(21,10),dispr3(21,10)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /BLANK / apress,atemp,strml,pgeom,iprtk,ifail,SIGN,zorign,  &
     fxcoor,fycoor,fzcoor
 COMMON /system/ sysbuf,nout
 COMMON /names / rd,rdrew,wrt,wrtrew,clsrew,norew
 COMMON /zzzzzz/ iz(1)
 COMMON /condas/ pi,twopi,radeg
 COMMON /unpakx/ typout,ir1,ir2,incr
 EQUIVALENCE     (iz(1),z(1)),(idata(1),rdata(1)), (ifill(1),rfill(1))
 DATA    NAME  / 4HALGP,4HR   /
 DATA    stream/ 3292, 92,292 /
 DATA    LEN   / 18, 24, 6    /
 DATA    iblk  , izero,rzero  / 4H    , 0, 0.0          /
 DATA    edt   , eqexin,ugv,algdd,cstm,bgpdt, scr1,scr2 /  &
     102   , 103   ,104,105  ,106 ,107  , 301 ,302  /
 
 
!     PERFORM GENERAL INITIALIZATION
 
 debug =.false.
 CALL sswtch (20,j)
 IF (j == 1) debug =.true.
 buf1 = korsz(iz) - sysbuf
 buf2 = buf1 - sysbuf
 left = corwds(iz(1),iz(buf2-1))
 m8   =-8
 IF (left <= 0) CALL mesage (m8,0,NAME)
 ir1  = 1
 incr = 1
 typout = 1
 ierr = 0
 
 ifill(1) = iblk
 ifill(2) = izero
 rfill(3) = rzero
 
!     CREATE ALGDB WITH CORRECT LENGTH RECORDS -
!     BCD(18 WORDS), INTEGER(24 WORDS), REAL(6 WORDS)
 
 CALL gopen (algdd,iz(buf1),rdrew)
 CALL gopen (scr2,iz(buf2),wrtrew)
 itrl(1) = algdd
 CALL rdtrl (itrl)
 itrl(1) = scr2
 CALL wrttrl (itrl)
 1 CALL READ (*7,*2,algdd,idata,99,1,nwar)
 2 CALL algpb (idata(1),ntype)
 length = LEN(ntype)
 
!     REMOVE NUMERIC ZEROS FROM BCD STRING
 
 IF (ntype /= 1) GO TO 4
 3 IF (idata(nwar) /= 0) GO TO 4
 nwar = nwar - 1
 IF (nwar >      0) GO TO 3
 4 IF (nwar >= length) GO TO 6
 nwar1 = nwar + 1
 DO  i = nwar1,length
   idata(i) = ifill(ntype)
 END DO
 6 CALL WRITE (scr2,idata,length,1)
 GO TO 1
 7 CALL CLOSE (algdd,clsrew)
 CALL CLOSE (scr2,clsrew)
 algdb = scr2
 
!     IF UGV IS NOT IN FIST (PURGED) THEN THERE WILL BE NO DATA
!     MODIFICATION
 
 itrl(1) = ugv
 CALL rdtrl (itrl)
 IF (itrl(1) < 0) GO TO 997
 
!     READ EQEXIN INTO CORE
 
 FILE = eqexin
 CALL gopen (eqexin,iz(buf1),rdrew)
 CALL READ (*901,*10,eqexin,iz(1),left,1,neqex)
 CALL mesage (m8,0,NAME)
 10 CALL fread (eqexin,iz(neqex+1),neqex,1)
 CALL CLOSE (eqexin,clsrew)
 kn = neqex/2
 IF (debug) CALL bug1 ('EQEX    ',10,iz(1),neqex)
 IF (debug) CALL bug1 ('EQEX    ',10,iz(neqex+1),neqex)
 
!     READ CSTM INTO CORE (CSTM MAY BE PURGED)
 
 FILE  = cstm
 icstm = 2*neqex + 1
 ncstm = 0
 CALL OPEN (*30,cstm,z(buf1),rdrew)
 CALL fwdrec (*901,cstm)
 CALL READ (*901,*20,cstm,iz(icstm),buf1-icstm,1,ncstm)
 CALL mesage (m8,0,NAME)
 20 CALL CLOSE (cstm,clsrew)
 IF (debug) CALL bug1 ('CSTM    ',20,iz(icstm),ncstm)
 
!     SET-UP FOR CALLS TO TRANSS
 
 CALL pretrs (iz(icstm),ncstm)
 
!     UNPACK UGV DISPLACEMENT VECTOR (SUBCASE 2) INTO CORE
 
 30 ivec = icstm + ncstm
 FILE = ugv
 itrl(1) = FILE
 CALL rdtrl (itrl)
 
!     CHECK FOR VALID UGV VECTOR
!     THIS ROUTINE WILL ONLY PROCESS A REAL S.P. RECT. VECTOR
!     OF SIZE G X 2
!     (EXPANDED TO INCLUDE REAL D.P. RECT. VECTOR, G X 2,
!     BY G.CHAN/UNISYS)
 
 nvects = itrl(2)
 kform  = itrl(4)
 ktype  = itrl(5)
 IF (nvects /= 2 .OR. kform /= 2) GO TO 902
 ivecn = ivec + ktype*itrl(3) - 1
 IF (ivecn >= buf1) CALL mesage (m8,0,NAME)
 
!     OPEN UGV AND SKIP FIRST COLUMN (SUBCASE 1)
 
 CALL gopen (ugv,iz(buf1),rdrew)
 CALL fwdrec (*901,ugv)
 ir2 = itrl(3)
 CALL unpack (*40,ugv,iz(ivec))
 GO TO 60
 
!     NULL COLUMN
 
 40 DO  i = ivec,ivecn
   z(i) = 0.0
 END DO
 60 CALL CLOSE (ugv,clsrew)
 IF (debug) CALL bug1 ('UGV     ',60,iz(ivec),ir2)
 
!     LOCATE STREAML1 CARDS ON EDT AND STORE IN CORE
 
 FILE   = edt
 ichord = ivecn + 1
 CALL preloc (*903,iz(buf1),edt)
 CALL locate (*904,iz(buf1),stream,idx)
 CALL READ (*901,*70,edt,iz(ichord),buf1-ichord,1,nchord)
 CALL mesage (m8,0,NAME)
 70 CALL CLOSE (edt,clsrew)
 IF (debug) CALL bug1 ('CHOR    ',70,iz(ichord),nchord)
 lchord = ichord + nchord -1
 
!     READ THE BGPDT INTO CORE
 
 ibgpdt = lchord + 1
 FILE   = bgpdt
 CALL gopen (bgpdt,iz(buf1),rdrew)
 CALL READ (*901,*80,bgpdt,iz(ibgpdt),buf1-ibgpdt,1,nbgpdt)
 CALL mesage (m8,0,NAME)
 80 CALL CLOSE (bgpdt,clsrew)
 IF (debug) CALL bug1 ('BGPD    ',80,iz(ibgpdt),nbgpdt)
 
!     FOR EACH STREAML1 CARD -
!     (1) FIND BLADE NODES
!     (2) FIND EQUIVALENT INTERNAL NUMBERS OF THESE NODES
!     (3) LOCATE CORRESPONDING COMPONENTS OF DISPLACEMENT AND
!         CONVERT THEN TO BASIC VIA CSTM
!     (4) LOCATE BASIC GRID POINT DATA FOR BLADE NODES
 
 ic  = ichord + 1
 icc = ichord
 jchord = 1
 nnodes = 0
 100 istatn = 0
 110 id = iz(ic)
 IF (id /= -1) GO TO 120
 icc = ic + 1
 ic  = ic + 2
 nnodes = nnodes + istatn
 jchord = jchord + 1
 IF (ic >= lchord) GO TO 150
 GO TO 100
 120 istatn = istatn + 1
 GO TO 1005
 
!     STORE BASIC GRID POINT COORDINATES FROM BGPDT
 
 130 xb(jchord,istatn) = z(icid+1)
 yb(jchord,istatn) = z(icid+2)
 zb(jchord,istatn) = z(icid+3)
 dispt1(jchord,istatn) = dispt(1)
 dispt2(jchord,istatn) = dispt(2)
 dispt3(jchord,istatn) = dispt(3)
 dispr1(jchord,istatn) = dispr(1)
 dispr2(jchord,istatn) = dispr(2)
 dispr3(jchord,istatn) = dispr(3)
 IF (debug) CALL bug1 ('NODE    ',id,z(icid+1),3)
 IF (debug) CALL bug1 ('NODE    ',id,dispt,3)
 IF (debug) CALL bug1 ('NODE    ',id,dispr,3)
 ic = ic + 1
 GO TO 110
 150 CONTINUE
 jchord = jchord - 1
 IF (jchord > 21) GO TO 906
 
!     MODIFY AERODYNAMIC INPUT  (OPEN ALGDB DATA BLOCK)
 
 FILE = algdb
 CALL gopen (algdb,iz(buf1),rdrew)
 CALL fwdrec (*907,algdb)
 CALL READ (*901,*908,algdb,idata,2,1,nwar)
 naero = idata(2)
 CALL skprec (algdb,1)
 CALL fread (algdb,idata,17,1)
 nlines = idata(1)
 nstns  = idata(2)
 nspec  = idata(4)
 ipunch = idata(8)
 isecn  = idata(9)
 ifcord = idata(10)
 isplit = idata(13)
 irle   = idata(15)
 irte   = idata(16)
 nsign  = idata(17)
 CALL skprec (algdb,1)
 DO  isk = 1,nstns
   CALL fread (algdb,idata,2,1)
   kptsa(isk) = idata(1)
   ifangs(isk)= idata(2)
   CALL skprec (algdb,idata(1))
   DO  inl = 1,nlines
     CALL fread (algdb,rdata,2,1)
     blafor(inl,isk) = rdata(2)
   END DO
 END DO
 DO  isk = 1,nspec
   CALL fread (algdb,rdata,6,1)
   zr(isk)  = rdata(1)
   jz(isk)  = rdata(1) + 0.4
   b1(isk)  = rdata(2)
   b2(isk)  = rdata(3)
   pp(isk)  = rdata(4)
   qq(isk)  = rdata(5)
   rle(isk) = rdata(6)
   CALL fread (algdb,rdata,6,1)
   tc(isk)  = rdata(1)
   te(isk)  = rdata(2)
   zed(isk) = rdata(3)
   cord(isk)= rdata(4)
   delx(isk)= rdata(5)
   dely(isk)= rdata(6)
   IF (isecn == 1 .OR. isecn == 3) CALL skprec (algdb,1)
 END DO
 CALL CLOSE (algdb,clsrew)
 
!     NUMBER OF BLADE STATIONS
 
 nblstn = irte - irle + 1
 IF (nlines /= jchord) GO TO 909
 IF (nnodes /= nlines*nblstn) GO TO 909
 
!     COMPUTE FCORD AND PHI
 
 DO  k = 1,nspec
   j    = jz(k)
   temp = (xb(j,nblstn)-xb(j,1))**2 + (zb(j,nblstn)-zb(j,1))**2
   IF (ifcord == 1) temp = temp   + (yb(j,nblstn)-yb(j,1))**2
   fchord(k) = cord(k)/SQRT(temp)
   phi(1,k)  = ATAN((zb(j,2)-zb(j,1))/(xb(j,2)-xb(j,1)))
   phi(2,k)  = ATAN((zb(j,nblstn)-zb(j,nblstn-1))/  &
       (xb(j,nblstn)-xb(j,nblstn-1)))
 END DO
!     COMPUTE NEW COORDINATES
!     GENERATE XSTA, RSTA AND R , SET KPTS = NLINES
 DO  i = 1,nlines
   DO  j = 1,nblstn
     xb(i,j) = xb(i,j) + SIGN*dispt1(i,j)*fxcoor
     yb(i,j) = yb(i,j) + SIGN*dispt2(i,j)*fycoor
     zb(i,j) = zb(i,j) + SIGN*dispt3(i,j)*fzcoor
     xsta(i,j) = xb(i,j)
     rsta(i,j) = zb(i,j) + zorign
     r(i,j)  = rsta(i,j)
   END DO
 END DO
 
!     COMPUTE CORD2
 
 DO  k = 1,nspec
   j    = jz(k)
   temp = (xb(j,nblstn)-xb(j,1))**2 + (zb(j,nblstn)-zb(j,1))**2
   IF (ifcord == 1) temp = temp   + (yb(j,nblstn)-yb(j,1))**2
   cord2(k) = fchord(k)*SQRT(temp)
 END DO
 
!     MODIFY B1, B2, RLE, TC, TE, CORD, DELX AND DELY
 
 i1 = (nblstn+1)/2
 i2 = i1
 IF (i1*2 /= nblstn+1) i2 = i2 + 1
 DO  k = 1,nspec
   j  = jz(k)
   b1(k) = b1(k) - nsign*SIGN*radeg*(dispr3(j,1)*COS(phi(1,k)) -  &
       dispr1(j,1)*SIN(phi(1,k)))
   b2(k) = b2(k) - nsign*SIGN*radeg*(dispr3(j,nblstn)*COS(phi(2,k)) -  &
       dispr1(j,nblstn)*SIN(phi(2,k)))
   temp   = cord(k)/cord2(k)
   rle(k) = rle(k) *temp
   tc(k)  = tc(k)  *temp
   te(k)  = te(k)  *temp
   cord(k) = cord2(k)
   delx(k) = delx(k) + 0.5*SIGN*fxcoor*(dispt1(j,i1)+dispt1(j,i2))
   dely(k) = dely(k) + 0.5*SIGN*fycoor*(dispt2(j,i1)+dispt2(j,i2))
 END DO
 
!     GENERATE NEW ALGDB DATA BLOCK
 
 CALL gopen (algdb,iz(buf1),rdrew)
 CALL gopen (scr1,iz(buf2),wrtrew)
 itrl(1) = algdb
 CALL rdtrl (itrl)
 
!     MODIFY THE NUMBER OF CARDS IN ALGDB
 
 ncdsx = 0
 DO  kpt = irle,irte
   ncdsx = ncdsx + nlines - kptsa(kpt)
 END DO
 itrl(2) = itrl(2) + ncdsx
 itrl(1) = scr1
 CALL wrttrl (itrl)
 ASSIGN 322 TO ret2
 nrec = 5
 GO TO 1300
 
!     COPY DATA FOR STATIONS 1 THRU (IRLE-1)
 
 322 IF (irle == 1) GO TO 335
 nles = irle - 1
 nrec = nles + nles*nlines
 DO  ikp = 1,nles
   nrec = nrec + kptsa(ikp)
 END DO
 ASSIGN 326 TO ret2
 GO TO 1300
 
!     SKIP OVER EXISTING RECORDS FOR STATIONS IRLE THRU IRTE
 
 326 nrec = nblstn + nblstn*nlines
 DO  ikp = irle,irte
   nrec = nrec + kptsa(ikp)
 END DO
 CALL skprec (algdb,nrec)
 
!     CREATE NEW DATA RECORDS FOR STATIONS IRLE THRU IRTE
 
 ksta = 0
 DO  jsta = irle,irte
   ksta = ksta + 1
   idata(1) = nlines
   idata(2) = ifangs(jsta)
   CALL WRITE (scr1,idata,2,1)
   IF (debug) CALL bug1 ('ALGPR   ',329,idata,2)
   DO  i = 1,nlines
     rdata(1) = xsta(i,ksta)
     rdata(2) = rsta(i,ksta)
     IF (debug) CALL bug1 ('ALGPR   ',330,rdata,2)
     CALL WRITE(scr1,rdata,2,1)
   END DO
   DO  i = 1,nlines
     rdata(1) = r(i,ksta)
     rdata(2) = blafor(i,ksta)
     IF (debug) CALL bug1 ('ALGPR   ',332,rdata,2)
     CALL WRITE (scr1,rdata,2,1)
   END DO
 END DO
 335 CONTINUE
 
!     COPY DATA FOR STATIONS (IRTE+1) THRU NSTNS
 
 IF (irte == nstns) GO TO 338
 irte1 = irte  + 1
 irte2 = nstns - irte
 nrec  = irte2 + irte2*nlines
 DO  ikp = irte1,nstns
   nrec  = nrec  + kptsa(ikp)
 END DO
 ASSIGN 338 TO ret2
 GO TO 1300
 338 CONTINUE
 
!     MODIFY THE NEXT NSPEC RECORDS
 
 DO  i = 1,nspec
   CALL skprec(algdb,2)
   rdata(1) = zr(i)
   rdata(2) = b1(i)
   rdata(3) = b2(i)
   rdata(4) = pp(i)
   rdata(5) = qq(i)
   rdata(6) = rle(i)
   CALL WRITE (scr1,rdata,6,1)
   IF (debug) CALL bug1 ('ALGPR   ',338,rdata,6)
   rdata(1) = tc(i)
   rdata(2) = te(i)
   rdata(3) = zed(i)
   rdata(4) = cord(i)
   rdata(5) = delx(i)
   rdata(6) = dely(i)
   CALL WRITE (scr1,rdata,6,1)
   IF (debug) CALL bug1 ('ALGPR   ',339,rdata,6)
   IF (isecn /= 1 .AND. isecn /= 3) CYCLE
   CALL fread (algdb,rdata,2,1)
   CALL WRITE (scr1,rdata,2,1)
   IF (debug) CALL bug1 ('ALGPR   ',340,rdata,2)
 END DO
 
!     COPY REST OF ANALYTIC DATA
 
 IF (isplit < 1) GO TO 344
 nrec = nspec
 DO  i = 1,nstns
   IF (ifangs(i) == 2) nrec = nrec + nlines
 END DO
 ASSIGN 344 TO ret2
 GO TO 1300
 344 CONTINUE
 IF (naero /= 1 .AND. ipunch /= 1) GO TO 352
 nrec  = 1
 ASSIGN 346 TO ret2
 GO TO 1300
 346 nrad  = idata(1)
 ndpts = idata(2)
 ndatr = idata(3)
 ASSIGN 347 TO ret2
 nrec  = 2
 GO TO 1300
 347 nb = nblstn - 1
 i  = 1
 348 nrec = 1
 ASSIGN 349 TO ret2
 GO TO  1300
 349 nrec = idata(1)
 ASSIGN 350 TO ret2
 GO TO  1300
 350 i = i + 1
 IF (i <= nb) GO TO 348
 nrec = nrad*(ndpts+1) + ndatr
 ASSIGN 352 TO ret2
 GO TO  1300
 
!     PROCESS AERODYNAMIC INPUT
 
 352 IF (naero == 0) GO TO 366
 ASSIGN 354 TO ret2
 nrec  = 3
 GO TO 1300
 354 nstns = idata(1)
 ncase = idata(6)
 nmany = idata(16)
 nle   = idata(19)
 nte   = idata(20)
 nsign = idata(21)
 IF (nstns == 0) nstns = 11
 IF (ncase == 0) ncase = 1
 nrec  = ncase + 3
 IF (nmany > 0) nrec = ncase + 4
 ASSIGN 356 TO ret2
 GO TO  1300
 356 CONTINUE
 
!     COPY DATA FOR STATIONS 1 THRU (NLE-1)
 
 IF (nle == 1) GO TO 361
 nle1 = nle - 1
 i    = 1
 357 nrec = 1
 ASSIGN 358 TO ret2
 GO TO  1300
 358 nrec = idata(1)
 ASSIGN 360 TO ret2
 GO TO  1300
 360 i = i + 1
 IF (i <= nle1) GO TO 357
 361 jsta = 0
 
!     MODIFY DATA FOR STATIONS NLE THRU NTE
 
 DO  i = nle,nte
   jsta = jsta + 1
   CALL fread (algdb,nspec,1,1)
   CALL skprec (algdb,nspec)
   CALL WRITE (scr1,nlines,1,1)
   IF (debug) CALL bug1 ('ALGPR   ',361,nlins,1)
   DO  nl = 1,nlines
     rdata(1) = xsta(nl,jsta)
     rdata(2) = rsta(nl,jsta)
     IF (debug) CALL bug1 ('ALGPR   ',362,rdata,2)
     CALL WRITE (scr1,rdata,2,1)
   END DO
 END DO
 
!     COPY REST OF DATA
 
 ASSIGN 366 TO ret2
 nrec = 65000
 GO TO  1300
 
!     CLOSE ALGDB AND SCR1
 
 366 CALL CLOSE (algdb,clsrew)
 CALL CLOSE (scr1,clsrew)
 
!     PUNCH NEW ALGDB TABLE INTO DTI CARDS IF PGEOM=3.
 
 IF (pgeom == 3) CALL algap (algdd,scr1)
 GO TO 999
 
 
!     INTERNAL BINARY SEARCH ROUTINE
 
!     SEARCH EQEXIN FOR INTERNAL NUMBER AND SIL NUMBER OF EXTERNAL NODE
 
 1005 klo = 1
 khi = kn
 1010 k   = (klo + khi + 1) / 2
 1020 IF (id - iz(2*k-1) < 0) THEN
   GO TO  1030
 ELSE IF (id - iz(2*k-1) == 0) THEN
   GO TO  1090
 ELSE
   GO TO  1040
 END IF
 1030 khi = k
 GO TO 1050
 1040 klo = k
 1050 IF (khi - klo - 1 < 0) THEN
   GO TO   905
 ELSE IF (khi - klo - 1 == 0) THEN
   GO TO  1060
 ELSE
   GO TO  1010
 END IF
 1060 IF (k == klo) GO TO 1070
 k   = klo
 GO TO 1080
 1070 k   = khi
 1080 klo = khi
 GO TO 1020
 1090 intn = iz(2*k)
 isil = iz(2*k+neqex)/10
 kode = iz(2*k+neqex) - 10*isil
 IF (debug) CALL bug1('ISTL    ',1090,isil,1)
 IF (debug) CALL bug1('KODE    ',1090,kode,1)
 
!     LOCATE COORDINATE SYSTEM ID FOR THIS NODE IN THE BGPDT
 
 icid = 4*(intn-1) + ibgpdt
 
!     SET-UP COORDINATE SYSTEM TRANSFORMATION FOR DISPLACEMENTS.
 
 IF (iz(icid) > 0) CALL transs (iz(icid),ta)
 
!     COMPUTE POINTER INTO UGV
!     JVEC = IVEC + KTYPE *(ISIL-1)
 
 jvec = ivec + typout*(isil-1)
 
!     PICK-UP DISPLACEMENTS
 
 IF (kode == 1) GO TO 1092
 
!     SCALAR POINT
 
 dispt(1) = z(jvec)
 dispt(2) = 0.0
 dispt(3) = 0.0
 dispr(1) = 0.0
 dispr(2) = 0.0
 dispr(3) = 0.0
 GO TO 1100
 
!     GRID POINT
 
 1092 IF (iz(icid) > 0) GO TO 1094
 
!     DISPLACEMENTS ALREADY IN BASIC SYSTEM
 
 dispt(1) = z(jvec  )
 dispt(2) = z(jvec+1)
 dispt(3) = z(jvec+2)
 dispr(1) = z(jvec+3)
 dispr(2) = z(jvec+4)
 dispr(3) = z(jvec+5)
 GO TO 1100
 
!     DISPLACEMENTS MUST BE TRANSFORMED TO BASIC
 
 1094 CALL gmmats (ta,3,3,0,z(jvec  ),3,1,0,dispt)
 CALL gmmats (ta,3,3,0,z(jvec+3),3,1,0,dispr)
 1100 CONTINUE
 GO TO 130
 1300 DO  icopy = 1,nrec
   CALL READ (*1306,*1302,algdb,idata,99,1,nwar)
   1302 CALL WRITE (scr1,idata,nwar,1)
   IF (debug) CALL bug1 ('ALGPR   ',1302,idata,nwar)
 END DO
 IF (nrec < 65000) GO TO 1306
 WRITE  (nout,1305)
 1305 FORMAT (/,' *** NO. OF RECORDS EXCEEDS HARDWARE LIMIT/ALGPR')
 CALL mesage (-37,0,0)
 1306 GO TO ret2, (322,326,338,344,346,347,349,350,352,354,356,358, 360,366)
 
 901 CALL mesage (-2,FILE,NAME)
 GO TO 998
 902 WRITE (nout,2001) ufm
 GO TO 998
 903 WRITE (nout,2002) ufm
 GO TO 998
 904 WRITE (nout,2003) ufm
 GO TO 998
 905 WRITE (nout,2004) ufm,iz(icc),id
 GO TO 998
 906 WRITE (nout,2005) uwm
 GO TO 999
 907 WRITE (nout,2006) ufm
 GO TO 998
 908 CALL mesage (-3,FILE,NAME)
 GO TO 998
 909 WRITE (nout,2007) ufm
 GO TO 998
 997 ierr = 1
 GO TO 999
 998 ierr = -1
 999 RETURN
 
 2001 FORMAT (a23,' - ALG MODULE - UGV DATA BLOCK IS NOT A REAL S.P. ',  &
     'RECTANGULAR MATRIX OF ORDER G BY 2.')
 2002 FORMAT (a23,' - ALG MODULE - EDT DATA BLOCK MAY NOT BE PURGED.')
 2003 FORMAT (a23,' - ALG MODULE - STREAML1 BULK DATA CARD MISSING ',  &
     'FROM BULK DATA DECK.')
 2004 FORMAT (a23,' - ALG MODULE - STREAML1 BULK DATA CARD (SLN NO. =',  &
     i3,') REFERENCES UNDEFINED NODE NO.',i8)
 2005 FORMAT (a25,' - ALG MODULE - MORE THAN 21 STREAML1 CARDS READ. ',  &
     'FIRST 21 WILL BE USED.')
 2006 FORMAT (a23,' - ALG MODULE - ALGDB DATA BLOCK (FILE 105) DOES ',  &
     'NOT HAVE ENOUGH RECORDS.')
 2007 FORMAT (a23,' - ALG MODULE - INPUT IN ALGDB DATA BLOCK (FILE 105',  &
     ') INCONSISTENT WITH DATA ON STREAML1 BULK DATA CARDS.',  &
     /39X,'CHECK THE NUMBER OF COMPUTING STATIONS AND THE ',  &
     'NUMBER OF STREAMSURFACES ON THE BLADE.')
END SUBROUTINE algpr
