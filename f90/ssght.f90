SUBROUTINE ssght
     
!     THIS IS THE STATIC-SOLUTION-GENERATOR FOR HEAT TRANSFER.
 
!     DMAP CALLING SEQUENCE.
 
!     SSGHT  USET,SIL,GPTT,GM,EST,MPT,DIT,PF,PS,KFF,KFS,KSF,KSS,RFN,RSN,
!            LFILE,UFILE/UGV,QG,RULV/V,N,NLK/V,N,NLR/C,Y,EPS0/C,Y,TABS/
!            C,Y,MAXITR/C,Y,IRES/V,N,MPCF1/V,N,SINGLE $
 
 LOGICAL :: nogo,noqg,rulvec,diagon,linear,loop1,nlrad
 INTEGER :: buf(10),sysbuf,outpt,tset,rd,rdrew,wrt,wrtrew,  &
     clsrew,cls,precis,core,single,eor,umcb,bmcb,xmcb,  &
     pkin,pkout,pkirow,pknrow,pkincr,eol,buf1,buf2,  &
     mcb(7),rulmcb(7),gsize,FILE,fsize,ssize,flag,word,  &
     mcb2(7),NAME(2),uset,gptt,gm,est,dit,ufile,pf,ps,  &
     rfn,rsn,ugv,qg,rulv,subr(2),ditx,z,scrt1,scrt2,  &
     scrt3,treqst,tstart,tend,tloop,scrt4,telaps, alibi(5,5)
 REAL :: rbuf(10),rz(1)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /fbsx  / jlmcb(7),jumcb(7),jbmcb(7),jxmcb(7),jzzz,jprec, jsign
 COMMON /gfbsx / lmcb(7),umcb(7),bmcb(7),xmcb(7),lz,iprec,ISIGN
 COMMON /packx / pkin,pkout,pkirow,pknrow,pkincr
 COMMON /zntpkx/ ai(4),irow,eol
 COMMON /zblpkx/ ao(4),jrow
 COMMON /system/ ksystm(65)
 COMMON /stime / treqst
 COMMON /hmatdd/ ihmat,nhmat,mptx,ditx
 COMMON /names / rd,rdrew,wrt,wrtrew,clsrew,cls
 COMMON /zzzzzz/ z(1)
 COMMON /BLANK / nlk,nlr,eps0,tabs,maxitr,ires,mpcf1,single
 EQUIVALENCE     (ksystm(1),sysbuf),(ksystm(2),outpt),  &
     (ksystm(10),tset),(ksystm(55),iprec1),  &
     (buf(1),rbuf(1)),(z(1),rz(1)),(idfalt,defalt)
 DATA    diagon/ .false. /
 DATA    subr  / 4HSSGH,4HT     /, eor,noeor /1,0/
 DATA    uset  , gptt,gm,est,mptfil / 101,103,104,105,106/
!     DATA    SIL   / 102     /
 DATA    dit   , pf,ps,kff,kfs,ksf,kss/107,108,109,110,111,112,113/
 DATA    rfn   , rsn,lfile,ufile  / 114,115,116,117/
 DATA    ugv   , qg,rulv / 201,202,203 /
 DATA    scrt1 , scrt2,scrt3,scrt4/ 301,302,303,304/
 DATA    alibi / 4H nor,4HMAL ,4HCONV,4HERGE,4HNCE ,  &
     4H MAX,4HIMUM,4H ite,4HRATI,4HONS , 4H div,4HERGI,4HNG s,4HOLUT,4HION ,  &
     4H ins,4HUFFI,4HCIEN,4HT ti,4HME  , 4H MAX,4HIMUM,4H con,4HVERG,4HENCE/
 
!     OPEN CORE
 
!          PRE-ITERATIONS              DURING-ITERATIONS
!         +--------------+
!         I              I  Z(IUNI)
!         I    I         I
!         I  (U ) VECTOR I
!         I    N         I
!         I              I  Z(NUNI)
!         +--------------+
!         I              I  Z(IHMAT)
!         I  HMAT CORE   I
!         I  BLOCK IF    I
!         I  REQUIRED    I
!         I              I  Z(NHMAT)
!         +--------------+
!         I              I  Z(IUN)
!         I (U ) PARTIT. I
!         I   N  VECTOR  I
!         I      FOR F+S I
!         I              I  Z(NUN)
!         +--------------+  -  -  -  -  -  -  -  -  -  -  -  -  -
!         I              I  Z(IEQIV)        MSIZE
!         I              I           -  -  -  *  -  -  EQUIV TABLE
!         I   E          I                   *I*       WILL BE PLACED
!         I(UN )EQIV.TBL I                  * I *      ON SCRATCH FILE
!         I   G          I                    I        DURING ITERATIONS
!         I              I  Z(NEQIV)          I
!         +--------------+  -  -  -   +----------------+
!         I              I  Z(IUM)    I                I  Z(ISN)
!         I (U )PARTIT.  I            I (S ) DIAGONAL  I
!         I   M VECTOR   I            I   N            I
!         I              I  Z(NUM)    I                I
!         +--------------+            I                I  Z(NSN)
!         I              I  Z(IUME)   +----------------+
!         I   E          I            I                I  Z(IDELU)
!         I (U ) TABLE   I            I (DELTA U ) VEC.I
!         I   M          I  Z(NUME)   I         N      I
!         +--------------+            I                I  Z(NDELU)
!                                     +----------------+
!                                     I                I  Z(IDELP)
!                                     I (DELTA P ) VEC.I
!                                     I         N      I
!                                     I                I  Z(NDELP)
!                                     +----------------+
!                                            .
!                                      CORE  . UNUSED
!                           Z(CORE)          .            Z(CORE)
!         - - - - - - - - - - - - - - +----------------+
!                                     I  BUFFER 2      I  Z(BUF2)
!                                     I                I
!                                     +----------------+
!                                     I  BUFFER 1      I  Z(BUF1)
!                                     I                I
!                                     +----------------+
 
!     CORE SIZE AND BUFFERS
 
 lcore = korsz(z)
 buf1  = lcore - sysbuf - 2
 buf2  = buf1  - sysbuf - 2
 core  = buf2  - 1
 IF (core < 100) CALL mesage (-8,0,subr)
 precis = 1
 
!     SET MISC. FLAGS.
 
 eps010 = 10.0*eps0
 epsold = eps0 + 1.0
 nlrad  = .true.
 IF (nlr == -1) nlrad = .false.
 CALL sswtch (18,k)
 IF (k == 1) diagon = .true.
 linear = .true.
 IF (nlk == +1) linear = .false.
 
!     READ TRAILER OF USET TO GET GSIZE.
 
 mcb(1) = uset
 CALL rdtrl (mcb)
 FILE = uset
 IF (mcb(1) <= 0) GO TO 1290
 gsize = mcb(3)
 
!     READ GM TRAILER TO DETERMINE COUNT OF UM POINTS
 
 mcb(1) = gm
 CALL rdtrl (mcb)
 IF (mcb(1) <= 0) GO TO 30
 msize = mcb(3)
 GO TO 40
 30 msize = 0
 40 nsize = gsize - msize
 
!     CORE ALLOCATION.
 
 iuni  = 1
 nuni  = nsize
 iuniz = iuni - 1
 ihmat = nuni + 1
 nhmat = nuni
 IF (.NOT.linear) nhmat = core
 mptx  = mptfil
 ditx  = dit
 IF (.NOT.linear) CALL prehma (z)
 iun   = nhmat + 1
 nun   = nhmat + nsize
 ieqiv = nun + 1
 neqiv = nun + gsize
 
!     EQUIVALENCE TABLE WILL BE PUT ON SCRATCH DURING ITERATIONS.
 
 isn   = nun + msize + 1
 nsn   = isn + nsize - 1
 IF (.NOT.nlrad) nsn = isn - 1
 isnz  = isn - 1
 idelu = nsn + 1
 ndelu = nsn + nsize
 ideluz= idelu - 1
 idelp = ndelu + 1
 ndelp = ndelu + nsize
 IF (ndelp > core) CALL mesage (-8,0,subr)
 idelpz= idelp - 1
 ium   = neqiv + 1
 num   = neqiv + msize
 iumz  = ium - 1
 iume  = num + 1
 nume  = num + msize
 iumez = iume - 1
 IF (nume > core) CALL mesage (-8,0,subr)
 
!     CONSTRUCTION OF (U ) AND (U ) TABLES.
!                       M        N
 
 FILE   = uset
 CALL gopen (uset,z(buf1),rdrew)
 mpoint = iumz
 npoint = iun - 1
 isil   = 0
 fsize  = 0
 CALL fread (uset,z(ieqiv),gsize,0)
 CALL CLOSE (uset,clsrew)
 DO  i = ieqiv,neqiv
   word = z(i)
   isil = isil + 1
   
!     CHECK FOR M-POINT
   
   IF (MOD(word,2) <= 0) GO TO 90
   mpoint = mpoint + 1
   z(mpoint) = isil
   CYCLE
   
!     ASSUME N-POINT
   
   90 npoint = npoint + 1
   
!     CHECK FOR F OR S POINT
   
   IF (MOD(word/2,2) > 0.0) THEN
     GO TO   110
   END IF
   
!     OK N-POINT IS AN F-POINT.
   
   100 z(npoint) = -isil
   fsize = fsize + 1
   CYCLE
   
!     OK N-POINT IS ASSUMED AN S-POINT
   
   110 z(npoint) = +isil
 END DO
 ssize = nsize - fsize
 
!     U  AND U  ARE COMPLETE.
!      M      N
 
 IF (isil == gsize .AND. npoint == nun .AND. mpoint == num) GO TO 140
 WRITE  (outpt,130) sfm
 130 FORMAT (a25,' 3081, INCONSISTENT USET DATA DETECTED.')
 CALL mesage (-61,0,subr)
 
!             E
!     BUILD (U ) EQUIVALENCE (U ) POINTS FOR (U ).
!             M                N               M
 
 140 IF (nume < iume) GO TO 250
 DO  i = iume,nume
   z(i) = 0
 END DO
 CALL gopen (gm,z(buf1),rdrew)
 DO  i = 1,nsize
   
!     OPERATE ON A COLUMN OF GM.
   
   CALL intpk (*210,gm,0,precis,0)
   160 CALL zntpki
   
!                             E
!     ROW POSITION -IROW- IN U  GETS COLUMN NUMBER.
!                             M
   
   ipos = iumez + irow
   IF (z(ipos) > 0.0) THEN
     GO TO   190
   END IF
   170 z(ipos) = i
   180 IF (eol > 0.0) THEN
     GO TO   210
   ELSE
     GO TO   160
   END IF
   
!     ERROR
   
   190 WRITE  (outpt,200) uwm,irow,i
   200 FORMAT (a25,' 3082, M =',i10,'  N =',i10)
   GO TO 180
   210 CONTINUE
 END DO
 CALL CLOSE (gm,clsrew)
 
!     INSURE ALL UME SLOTS FILLED
 
 nogo = .false.
 DO  i = iume,nume
   IF (z(i) > 0.0) THEN
     GO TO   240
   END IF
   220 m = i - iumez
   isil = iumz + m
   WRITE  (outpt,230) ufm,m,z(isil)
   230 FORMAT (a23,' 3083, UM POSITION =',i10,', SIL =',i10)
   nogo = .true.
   240 CONTINUE
 END DO
 IF (nogo) CALL mesage (-61,0,subr)
 
!                        E
!     CONSTRUCTION OF (UN ) EQUIVALENCE TABLE.
!                        G
 
 250 mpoint = ium
 mpt  = z(mpoint)
 IF (mpoint > num ) mpt = 1000000
 mpte = iume
 mval = z(mpte)
 nval = 1
 k = ieqiv - 1
 DO  i = 1,gsize
   k = k + 1
   IF (i /= mpt) GO TO 260
   
!     M-POINT NEXT
   
   z(k) = mval
   mpte = mpte + 1
   mval = z(mpte)
   mpoint = mpoint + 1
   mpt  = z(mpoint)
   IF (mpoint > num ) mpt = 1000000
   CYCLE
   
!     N-POINT NEXT
   
   260 z(k) = nval
   nval = nval + 1
 END DO
 
!     SET UP RULV IF RESIDUAL LOAD MATRIX IS TO BE FORMED.
 
 rulvec = .false.
 IF (ires <= 0) GO TO 290
 CALL makmcb (rulmcb,rulv,fsize,2,precis)
 CALL gopen (rulv,z(buf1),wrtrew)
 CALL CLOSE (rulv,cls)
 rulvec = .true.
 
!     GRID POINT TEMPERATURE DATA IS EXPANDED INTO CORE NOW.  ONLY
 
!      1
!     U  IS FORMED.
!      N
 
 
 290 IF (tset > 0.0) THEN
   GO TO   310
 END IF
 300 k = 0
 GO TO 320
 310 k = 1
 320 DO  i = iuni,nuni
   z(i) = k
 END DO
 IF (tset > 0.0) THEN
   GO TO   340
 ELSE
   GO TO   510
 END IF
 
!     POSITION GPTT TO GRID TEMPERATURE DATA SECTION.
 
 340 FILE = gptt
 CALL OPEN (*1290,gptt,z(buf1),rdrew)
 CALL fread (gptt,buf,-2,0)
 NUMBER = 0
 350 CALL READ (*1300,*360,gptt,buf,3,noeor,flag)
 NUMBER = MAX0(NUMBER,buf(3))
 GO TO 350
 360 CALL skprec (gptt,NUMBER)
 
!     NOW AT GRID TEMP SECTION HEADER.
 
 CALL fread (gptt,buf,-2,0)
 400 CALL READ (*1300,*1330,gptt,buf,3,noeor,flag)
 IF (tset /= buf(1)) GO TO 400
 
!     BUF(1)=SET-ID,  BUF(2)=-1 OR DEFAULT TEMP,  BUF(3)=GPTT RECORD.
 
 defalt = rbuf(2)
 IF (buf(3) <= 0) GO TO 470
 CALL skprec (gptt,buf(3))
 
!     TEMP PAIRS IN INTERNAL-ID AND TEMPERATURE.
 
 iunat = iun
 isil  = IABS(z(iunat))
 iuniat= iuni
 
!     READ A TEMPERATURE PAIR.
 
 430 CALL READ (*1300,*470,gptt,buf,2,noeor,flag)
 440 IF (buf(1)-isil < 0.0) THEN
   GO TO   430
 ELSE IF (buf(1)-isil == 0.0) THEN
   GO TO   450
 ELSE
   GO TO   460
 END IF
 450 z(iuniat) = buf(2)
 460 iunat = iunat + 1
 isil  = IABS(z(iunat))
 iuniat= iuniat + 1
 IF (iuniat <= nuni) GO TO 440
 470 CALL CLOSE (gptt,clsrew)
 
!     CHECK FOR INTEGER 1-S WHICH GET THE DEFAULT TEMP.
 
 nogo = .false.
 DO  i = iuni,nuni
   IF (z(i) /= 1) CYCLE
   IF (idfalt /= -1) GO TO 490
   nogo = .true.
   k = iun + i - iuni
   isil = IABS(z(k))
   WRITE  (outpt,480) ufm,isil
   480 FORMAT (a23,' 3084, THERE IS NO TEMPERATURE DATA FOR SIL NUMBER', i10)
   CYCLE
   490 rz(i) = defalt
 END DO
 IF (nogo) CALL mesage (-61,0,subr)
 510 CONTINUE
 
!               1                  1
!     COMPUTE (P ) = (P ) - (K  )(U ) AND SAVE ON SCRATCH-4.
!               F      F      FS   S
 
 k = idelpz + fsize
 DO  i = idelp,k
   z(i) = 0
 END DO
 CALL OPEN (*540,pf,z(buf1),rdrew)
 CALL fwdrec (*540,pf)
 CALL intpk (*540,pf,0,precis,0)
 530 CALL zntpki
 k = idelpz + irow
 rz(k) = ai(1)
 IF (eol > 0.0) THEN
   GO TO   540
 ELSE
   GO TO   530
 END IF
 540 CALL CLOSE (pf,clsrew)
 
!                         1
!     SUBTRACT OFF (K  )(U )
!                    FS   S
 
 iat = iun - 1
 CALL OPEN (*590,kfs,z(buf1),rdrew)
 CALL fwdrec (*590,kfs)
 DO  i = 1,ssize
   
!     FIND NEXT US POINT TEMPERATURE DATA.
   
   550 iat = iat + 1
   IF (z(iat) > 0.0) THEN
     GO TO   560
   ELSE
     GO TO   550
   END IF
   560 k = iuniz + iat - iun + 1
   CALL intpk (*580,kfs,0,precis,0)
   value = rz(k)
   570 CALL zntpki
   k = idelpz + irow
   rz(k) = rz(k) - ai(1)*value
   IF (eol > 0.0) THEN
     GO TO   580
   ELSE
     GO TO   570
   END IF
   580 CONTINUE
 END DO
 590 CALL CLOSE (kfs,clsrew)
 
!                1
!     PACK OUT (P ) ON SCRATCH-4
!                F
 
 CALL gopen (scrt4,z(buf1),wrtrew)
 CALL makmcb (mcb,scrt4,fsize,2,precis)
 pkin   = precis
 pkout  = precis
 pkirow = 1
 pknrow = fsize
 pkincr = 1
 CALL pack (z(idelp),scrt4,mcb)
 CALL CLOSE (scrt4,clsrew)
 CALL wrttrl (mcb)
 
!     ELEMENT INITIAL PROCESSING PHASE.
 
 CALL gopen  (scrt1,z(buf2),wrtrew)
 IF (linear) GO TO 600
 CALL gopen  (est,z(buf1),rdrew)
 CALL ssght1 (est,scrt1,z(ieqiv))
 CALL CLOSE  (est,clsrew)
 
!        E
!     (UN ) EQUIVALENCE TABLE IS NOW APPENDED TO -SCRT1-.
!        G
 
 600 CALL WRITE (scrt1,0,0,1)
 CALL WRITE (scrt1,z(ieqiv),gsize,1)
 CALL CLOSE (scrt1,clsrew)
 
!                      1            3
!     FORM (S ) = 4( (U ) + (TABS) )  DIAGONAL MATRIX.
!            N         N
 
 IF (.NOT.nlrad) GO TO 630
 j = iuniz
 DO  i = isn,nsn
   j = j + 1
   rz(i) = 4.0*(rz(j) + tabs)**3
 END DO
 
!     SET PARTITIONING TABLE IN TERMS OF WHERE ELEMENTS ARE TO MOVE TO
!     WHEN GOING FROM N-SET TO F+S SETS.
 
 630 is = fsize
 IF = 0
 DO  i = iun,nun
   IF (z(i) > 0.0) THEN
     GO TO   650
   END IF
   
!     F-POINTER
   
   640 IF   = IF + 1
   z(i) = IF
   CYCLE
   
!     S-POINTER
   
   650 is = is + 1
   z(i)  = is
 END DO
 loop  = 0
 loop1 = .true.
 pfmag = 0.0
 
!     == ITERATION SECTION ==
 
!     ITERATIVE LOOPING
 
 670 loop = loop + 1
 
!     TIME LEFT AT START OF LOOP
 
 CALL tmtogo (tstart)
 DO  i = idelp,ndelp
   z(i) = 0
 END DO
 IF (loop1 .OR. linear) GO TO 690
 CALL gopen (scrt1,z(buf1),rdrew)
 CALL ssght2 (scrt1,z(idelp),z(iuni))
 CALL CLOSE (scrt1,clsrew)
 
!     PARTITION DELTA-P VECTOR INTO DELTA-F AND DELTA-S VECTORS.
 
 CALL ssghtp (z(iun),z(idelp),nsize)
 
!                     I
!     GENERATION OF (N ) WILL BE PERFORMED IN CORE SPACE OF (DELTA-P)
!                     F
!       I                          I        4         I
!     (N ) = (DELTA-P ) + (R  )( (U  + TABS)  - (S )(U ) )
!       F            F      FN     N               N  N
 
 690 IF (.NOT.nlrad) GO TO 730
 CALL OPEN (*720,rfn,z(buf2),rdrew)
 CALL fwdrec (*720,rfn)
 DO  i = 1,nsize
   
!     OPERATE ON A COLUMN OF RFN
   
   CALL intpk (*710,rfn,0,precis,0)
   
!     COMPUTE CONSTANT FOR COLUMN
   
   k  = iuniz + i
   un = rz(k)
   k  = isnz + i
   sn = rz(k)
   value = (un + tabs)**4 - sn*un
   
!     UNPACK NON-ZERO TERMS OF COLUMN.
   
   700 CALL zntpki
   k = idelpz + irow
   rz(k) = rz(k) + ai(1)*value
   IF (eol > 0.0) THEN
     GO TO   710
   ELSE
     GO TO   700
   END IF
   710 CONTINUE
 END DO
 720 CALL CLOSE (rfn,clsrew)
 
!          I      1      I
!     (PBAR ) = (P ) - (N )
!          F      F      F
!                    I
!     FIRST NEGATE (N ) SITTING IN DELTA-P CORE SPACE,
!                    F
!                                     1
!     THEN ADD IN NON-ZERO TERMS OF (P )
!                                     F
 
 730 k = idelpz + fsize
 DO  i = idelp,k
   rz(i) = -rz(i)
 END DO
 
!            1
!     OPEN (P ) FOR UNPACKING OF ONE COLUMN.
!            F
 
 CALL OPEN (*760,scrt4,z(buf2),rdrew)
 CALL fwdrec (*760,scrt4)
 CALL intpk (*760,scrt4,0,precis,0)
 750 CALL zntpki
 k = idelpz + irow
 rz(k) = rz(k) + ai(1)
 IF (loop1) pfmag = pfmag + ai(1)*ai(1)
 IF (eol > 0.0) THEN
   GO TO   760
 ELSE
   GO TO   750
 END IF
 760 CALL CLOSE (scrt4,clsrew)
 
!          I
!     (PBAR ) IS NOW PACKED OUT TO SCRATCH-2.
!          F
 
 IF (.NOT.loop1) GO TO 790
 pfmag = SQRT(pfmag)
 IF (pfmag > 0.0) THEN
   GO TO   790
 END IF
 770 WRITE  (outpt,780) ufm
 780 FORMAT (a23,' 3085, THE PF LOAD VECTOR IS EITHER PURGED OR NULL.')
 CALL mesage (-61,0,subr)
 790 CALL makmcb (mcb2,scrt2,fsize,2,2)
 CALL gopen (scrt2,z(buf2),wrtrew)
 pkin   = precis
 pkout  = iprec1
 pkirow = 1
 pknrow = fsize
 pkincr = 1
 CALL pack (z(idelp),scrt2,mcb2)
 CALL CLOSE (scrt2,clsrew)
 CALL wrttrl (mcb2)
 
!                       I           I
!     (DELTA-P ) = (PBAR ) - (K  )(U )
!             F         F      FF   F
!          I
!     (PBAR ) IS SITING IN CORE CURRENTLY.  (IT WILL BE GONE TOMORROW.)
!          F
!                       I       I        I
!     FIRST PARTITION (U ) TO (U ) AND (U )
!                       N       F        S
 
 CALL ssghtp (z(iun),z(iuni),nsize)
 CALL OPEN (*820,kff,z(buf1),rdrew)
 CALL fwdrec (*820,kff)
 DO  i = 1,fsize
   
!     OPERATE ON ONE COLUMN OF KFF
   
   CALL intpk (*810,kff,0,precis,0)
   
!                                 I
!     LOCATE COLUMN MULTIPLIER = U
!                                 FI
   
   k = iuniz + i
   value = rz(k)
   800 CALL zntpki
!                                                            I
!     SUBTRACT THIS ELEMENT*VALUE FROM IROW POSITION OF (PBAR )
!                                                            F
   k = idelpz + irow
   rz(k) = rz(k) - ai(1)*value
   IF (eol > 0.0) THEN
     GO TO   810
   ELSE
     GO TO   800
   END IF
   810 CONTINUE
 END DO
 820 CALL CLOSE (kff,clsrew)
 
!     COMPUTE EPSILON
!                    P
 
 k   = idelpz + fsize
 sum = 0.0
 DO  i = idelp,k
   sum = sum + rz(i)**2
 END DO
 sum = SQRT(sum)
 epsubp = sum/pfmag
 IF (loop1 .AND. diagon) WRITE (outpt,840) epsubp
 840 FORMAT ('1D I A G   1 8   O U T P U T   F R O M   S S G H T', //,  &
     ' ITERATION    EPSILON-P',9X,'LAMBDA-1',10X,'EPSILON-T',  &
     /1X,60(1H=), /,6H     1,1P,e19.6)
 
!                                                   I
!     IF -RULV- IS BEING FORMED, THEN WRITE (DELTA-P ) OUT ON -RULV-.
!                                                   F
 
 IF (.NOT. rulvec) GO TO 850
 CALL OPEN (*850,rulv,z(buf1),wrt)
 pkin   = precis
 pkout  = precis
 pkirow = 1
 pknrow = fsize
 pkincr = 1
 CALL pack (z(idelp),rulv,rulmcb)
 CALL CLOSE (rulv,cls)
 
!                     I+1
!     NOW SOLVE FOR (U   ) IN,
!                     F
!                           I+1         I
!                   (L)(U)(U   ) = (PBAR )
!                           F           F
 
 
 850 ISIGN   =+1
 iprec   = 2
 lmcb(1) = lfile
 CALL rdtrl (lmcb)
 umcb(1) = ufile
 CALL rdtrl (umcb)
 bmcb(1) = scrt2
 CALL rdtrl (bmcb)
 CALL makmcb (xmcb,scrt3,fsize,2,2)
 
!     INSURE EVEN BOUNDARY (ARRAY WILL BE USED AS DOUBLE PRECISION)
 
 jdelp = ndelp + 1 + MOD(ndelp+1,2) + 1
 lz = lcore - jdelp
!WKBI 3/94
 jzzz = lz
 DO  ijk = 1,31
   jlmcb(ijk) = lmcb(ijk)
 END DO
 IF (umcb(1) > 0) CALL gfbs (z(jdelp),z(jdelp))
 IF (umcb(1) <= 0) CALL  fbs (z(jdelp),z(jdelp))
 IF (umcb(1) > 0) CALL wrttrl( xmcb)
 IF (umcb(1) <= 0) CALL wrttrl(jxmcb)
 
!       I+1
!     (U   ) IS NOW MOVED FROM SCRATCH-3 INTO CORE IN (DELTA-P ) SPACE.
!       F                                                     N
 
 CALL gopen (scrt3,z(buf1),rdrew)
 k = idelpz + fsize
 DO  i = idelp,k
   z(i) = 0
 END DO
 CALL intpk (*880,scrt3,0,precis,0)
 870 CALL zntpki
 k = idelpz + irow
 rz(k) = ai(1)
 IF (eol > 0.0) THEN
   GO TO   880
 ELSE
   GO TO   870
 END IF
 880 CALL CLOSE (scrt3,clsrew)
 IF (loop1) GO TO 985
 
!                      I+1       I
!     ALPHA = SUM OF (U   ) (PBAR )          IROW = 1,FSIZE
!                      F         F
!                       IROW      IROW
 
!                           I+1      I
!     BETA = SUM OF (DELTA-U   )(PBAR )      IROW = 1,FSIZE
!                           F        F
!                            IROW     IROW
 
!                      I+1    I       I
!     GAMMA = SUM OF (U    - U  )(PBAR )     IROW = 1,FSIZE
!                      F      F       F
!                       IROW   IROW    IROW
 
!     WHERE I = ITERATION GREATER THAN 1.
 
 CALL gopen (scrt2,z(buf1),rdrew)
 alpha = 0.0
 beta  = 0.0
 gamma = 0.0
 CALL intpk (*900,scrt2,0,precis,0)
 
!     ONLY NON-ZERO TERMS OF (PBAR ) NEED BE CONSIDERED.
!                                 F
 890 CALL zntpki
 kufip1 = idelpz + irow
 kdelu  = ideluz + irow
 kufi   = iuniz  + irow
 alpha  = alpha  + rz(kufip1)*ai(1)
 beta   = beta   + rz(kdelu) *ai(1)
 gamma  = gamma  + (rz(kufip1) - rz(kufi))*ai(1)
 IF (eol > 0.0) THEN
   GO TO   900
 ELSE
   GO TO   890
 END IF
 900 CALL CLOSE (scrt2,clsrew)
 
!     CONVERGENCE TESTS ARE MADE HERE.
 
!     WHEN ENTERING EXIT MODE,
!         -IEXIT-        -REASON-
!            1           NORMAL CONVERGENCE
!            2           NO CONVERGENCE AT MAXIMUM ITERATIONS
!            3           NO CONVERGENCE UNSTABLE ITERATION
!            4           NO CONVERGENCE INSUFFICIENT TIME
!            5           MAXIMUM CONVERGENCE, BUT EPSHT NOT SATISFIED
 
 IF (gamma == 0.0) THEN
   GO TO   920
 END IF
 910 flamda = ABS(beta/gamma)
 GO TO 930
 920 flamda = 100.0
 epst   = 0.0
 GO TO 970
 930 IF (alpha == 0.0) THEN
   GO TO   960
 END IF
 940 IF (flamda-1.0 == 0.0) THEN
   GO TO   960
 END IF
 950 epst = ABS(gamma/((flamda - 1.0)*alpha))
 GO TO 970
 960 epst = 100.0
 970 CALL tmtogo (kleft)
 telaps = treqst - kleft
 tau    = 1.0 - FLOAT(tloop+telaps)/(.8*FLOAT(treqst))
 IF (diagon) WRITE (outpt,980) loop,epsubp,flamda,epst
 980 FORMAT (i6,1P,e19.6,1P,e18.6,1P,e18.6)
 iexit  = 1
 IF (epst < eps0 .AND. flamda > 1.0 .AND. epsubp < eps010) GO TO 1060
 
!     TEST FOR TWO SUCCESSIVE CASES PASSING TEST
 
 IF (epst < eps0 .AND. epsold < eps0) GO TO 1060
 epsold = epst
 iexit  = 2
 IF (loop >= maxitr) GO TO 1060
 iexit  = 3
 IF (flamda <= 1.0 .AND. loop >= 4) GO TO 1060
 iexit  = 5
 IF (gamma == 0.) GO TO 1060
 iexit  = 4
 IF (tau < 0.0) THEN
   GO TO  1060
 ELSE
   GO TO   990
 END IF
 
!                   I
!     COMPUTE (DELTA ) TO BE USED ON NEXT LOOP
!                   U
 
 985 iexit = 2
 IF (loop >= maxitr) GO TO 1051
 990 k  = idelpz + fsize
 kdelu = ideluz
 ki = iuniz
 kip1 = idelpz
 DO  i = idelp,k
   kdelu = kdelu + 1
   ki = ki + 1
   kip1 = kip1 + 1
   rz(kdelu) = rz(kip1) - rz(ki)
 END DO
 
!                       I+1
!     MOVE (U ) UNDER (U   ) BOTH TO BE IN (DELTA-P ) CORE.
!            S          F                          N
 
 ASSIGN 1050 TO iretrn
 1010 k1 = iuni + fsize
 k2 = idelpz + fsize
 IF (ssize <= 0) GO TO 1030
 DO  i = k1,nuni
   k2 = k2 + 1
   rz(k2) = rz(i)
 END DO
 
!             I+1                       I+1
!     MERGE (U   ) AND (U ) BACK INTO (U   ) FORM.
!             F          S              N
 
 1030 kuni = iuniz
 DO  i = iun,nun
   kuni = kuni + 1
   jpos = idelpz + z(i)
   z(kuni) = z(jpos)
 END DO
 GO TO iretrn, (1050,1180)
 
!     READY NOW FOR ANOTHER LOOP.
 
 1050 CALL tmtogo (tend)
 tloop = tstart - tend
 loop1 = .false.
 GO TO 670
 
!     == END ITERATION SECTION ==
 
!     ITERATION HALTED, NOW IN EXIT MODE.
!     IF QG FILE IS PRESENT, FORCES OF CONSTRAINT ARE PARTIALLY COMPUTED
!     QS WILL BE FORMED IN THE CORE SPACE USED UP TO NOW FOR (DELTA-U).
 
!                           I
!     (Q ) = -(P ) + (K  )(U ) + (K  )(U ) + (DELTA-P ) + (PRODUCT )
!       S       S      SF   F      SS   S            S            S
 
!                                 I         4        I
!     WHERE (PRODUCT ) = (R  )( (U   + TABS)  - (S  U ) )
!                   S      SN     NJ              NJ N
 
!                                      J = 1,NSIZE
 
!     LOAD (DELTA-P ) INTO QS FORMATION CORE SPACE.
!                  S
 
 1051 WRITE  (outpt,1052) uwm
 1052 FORMAT (a25,' 3132, SSGHT RECOVERING FROM SEVERE USER CONVERGENCE'  &
     ,       ' CRITERIA.')
 1060 WRITE  (outpt,1070) uim,iexit,(alibi(j,iexit),j=1,5)
 1070 FORMAT (a29,' 3086, ENTERING SSGHT EXIT MODE BY REASON NUMBER ',  &
     i2,2H (,5A4,1H) )
 noqg = .true.
 CALL OPEN (*1170,qg,z(buf2),wrtrew)
 noqg = .false.
 CALL fname (qg,NAME)
 CALL WRITE (qg,NAME,2,eor)
 iqs  = idelu
 nqs  = ideluz + ssize
 iqsz = ideluz
 k    = idelpz + fsize
 DO  i = iqs,nqs
   k    = k + 1
   z(i) = z(k)
 END DO
 
!     SUBTRACT OFF NON-ZERO TERMS OF PS VECTOR.
 
 CALL OPEN (*1100,ps,z(buf1),rdrew)
 CALL fwdrec (*1100,ps)
 CALL intpk (*1100,ps,0,precis,0)
 1090 CALL zntpki
 k = iqsz + irow
 rz(k) = rz(k) - ai(1)
 IF (eol > 0.0) THEN
   GO TO  1100
 ELSE
   GO TO  1090
 END IF
 1100 CALL CLOSE (ps,clsrew)
 
!                   I
!     ADD IN (K  )(U )
!              SF   F
 
 CALL OPEN (*1130,ksf,z(buf1),rdrew)
 CALL fwdrec (*1130,ksf)
 DO  i = 1,fsize
   CALL intpk (*1120,ksf,0,precis,0)
   k = idelpz + i
   value = rz(k)
   1110 CALL zntpki
   k = iqsz + irow
   rz(k) = rz(k) + ai(1)*value
   IF (eol > 0.0) THEN
     GO TO  1120
   ELSE
     GO TO  1110
   END IF
   1120 CONTINUE
 END DO
 1130 CALL CLOSE (ksf,clsrew)
 
!     ADD IN (K  )(U )
!              SS   S
 
 IF (ssize == 0) GO TO 1160
 CALL OPEN (*1160,kss,z(buf1),rdrew)
 CALL fwdrec (*1160,kss)
 iusz = iuniz + fsize
 DO  i = 1,ssize
   CALL intpk (*1150,kss,0,precis,0)
   k = iusz + i
   value = rz(k)
   1140 CALL zntpki
   k = iqsz + irow
   rz(k) = rz(k) + ai(1)*value
   IF (eol > 0.0) THEN
     GO TO  1150
   ELSE
     GO TO  1140
   END IF
   1150 CONTINUE
 END DO
 1160 CALL CLOSE (kss,clsrew)
 
!                                     I
!     TO COMPUTE ADDITIONAL PRODUCT (U ) IS NOW FORMED.
!                                     N
!                           I
!     THUS MERGE (U ) AND (U )
!                  S        F
!                                  I
!     FIRST MOVE (U ) DOWN UNDER (U ), THEN DO MERGE.
!                  S               F
 
 1170 ASSIGN 1180 TO iretrn
 GO TO 1010
 
!     OK FORM AND ADD (PRODUCT) IN.
 
 1180 IF (.NOT.nlrad) GO TO 1220
 IF (noqg) GO TO 1250
 CALL OPEN (*1210,rsn,z(buf1),rdrew)
 CALL fwdrec (*1210,rsn)
 DO  i = 1,nsize
   CALL intpk (*1200,rsn,0,precis,0)
   ku = iuniz + i
   ks = isnz + i
   value = (rz(ku) + tabs)**4 - rz(ku)*rz(ks)
   1190 CALL zntpki
   k = iqsz + irow
   rz(k) = rz(k) + ai(1)*value
   IF (eol > 0.0) THEN
     GO TO  1200
   ELSE
     GO TO  1190
   END IF
   1200 CONTINUE
 END DO
 1210 CALL CLOSE (rsn,clsrew)
 
!     (QS) IS COMPLETE AND READY FOR EXPANSION TO GSIZE AND OUTPUT.
 
 1220 CALL makmcb (mcb,qg,gsize,2,precis)
 jrow = 0
 FILE = uset
 IF (ssize == 0) GO TO 1250
 CALL gopen (uset,z(buf1),rdrew)
 iq = iqs
 CALL bldpk (precis,precis,qg,0,0)
 1230 CALL fread (uset,word,1,0)
 jrow = jrow + 1
 IF (MOD(word/2,2) > 0.0) THEN
   GO TO  1240
 ELSE
   GO TO  1230
 END IF
 1240 ao(1) = rz(iq)
 CALL zblpki
 iq = iq + 1
 IF (iq <= nqs) GO TO 1230
 
!     QS HAS NOW BEEN EXPANDED TO GSIZE AND OUTPUT ON QG DATA BLOCK.
 
 CALL bldpkn (qg,0,mcb)
 CALL CLOSE  (qg,clsrew)
 CALL wrttrl (mcb)
 CALL CLOSE  (uset,clsrew)
 
!     PACK OUT (U ) USING THE EQUIVALENCE TABLE TO ORDER
!                G
 
!     THE U  POINTS.
!          N
 
 
!     READ EQUIVALENCE TABLE BACK INTO CORE AT THIS TIME.
 
 1250 FILE = scrt1
 CALL gopen  (scrt1,z(buf1),rdrew)
 CALL skprec (scrt1,1)
 CALL fread  (scrt1,z(ieqiv),gsize,0)
 
 CALL CLOSE  (scrt1,clsrew)
 
!     REPLACE POINTERS WITH THE VALUES.
 
 DO  i = ieqiv,neqiv
   k = iuniz + z(i)
   rz(i) = rz(k)
 END DO
 
!     PACK OUT (U )
!                G
 
 CALL makmcb (mcb,ugv,gsize,2,precis)
 CALL gopen  (ugv,z(buf1),1)
 pkin   = precis
 pkout  = precis
 pkirow = 1
 pknrow = gsize
 pkincr = 1
 CALL pack (z(ieqiv),ugv,mcb)
 CALL CLOSE (ugv,clsrew)
 CALL wrttrl (mcb)
 
!     COMPLETE RULV IF NECESSARY.
 
 IF (.NOT.rulvec) GO TO 1280
 CALL gopen (rulv,z(buf1),3)
 CALL CLOSE (rulv,clsrew)
 CALL wrttrl (rulmcb)
 1280 RETURN
 
!     ERROR CONDITIONS
 
 1290 n = -1
 GO TO 1320
 1300 n = -2
 GO TO 1320
 1320 CALL mesage (n,FILE,subr)
 1330 WRITE  (outpt,1340) ufm,tset
 1340 FORMAT (a23,' 3087, TEMPERATURE SET',i10,' IS NOT PRESENT IN ',  &
     'GPTT DATA BLOCK.')
 CALL mesage (-61,0,subr)
 RETURN
END SUBROUTINE ssght
