SUBROUTINE eqmckm
     
!     THIS SUBROUTINE CALCULATES THE MPC CONSTRAINT FORCES AND CREATES
!     THE OUTPUT FILE FOR OFP.
!     TASKS INCLUDE CREATING THE SCRATCH FILES FOR THE CURRENT SUBCASES
!     (PGG, QG - ALSO USED IN EQUILIBRIUM CHECKS).
!     NOT CODED TO HANDLE CONICAL ELEMENTS OR SORT2.
 
 LOGICAL :: firstc,firsto,lascas,anyout
 INTEGER :: NAME(2),kon(10),idat(3),parm,trl,ug,um,un,mcb(7),  &
     rdnrw,rdrw,wrtnrw,wrtrw,zz,ocb(8)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /names / rdnrw,rdrw,wrtnrw,wrtrw,krw,knrw,knerw
 COMMON /system/ ksystm(80)
 COMMON /bitpos/ um,skps(8),un,ug
 COMMON /BLANK / skpb(2),nskip
 COMMON /eqmk1 / kscc,kexin,kgpl,kbgdt,ksil,kuset,kgg,kgm,kugv,  &
     kpgg,kqg,kcstm,klam,koqm,kscr(7),kmpc,kload,kspc, parm(4),trl(7)
!ZZ   COMMON /ZZSSA2/ ZZ(1)
 COMMON /zzzzzz/ zz(20000)
 COMMON /unpakx/ itypu,inru,ilru,incu
 COMMON /mpyadx/ ma(7),mb(7),mc(7),md(7),mz,mt,msab,msc,mpr,mscr
 COMMON /patx  / lcor,ins(3),luset
 EQUIVALENCE     (mcb(1),ocb(1)),(ksystm(1),isbz),(ksystm(2),nout),  &
     (ksystm(15),itim),(ksystm(16),idat(1))
!WKBI 3/94 SPR93007
 EQUIVALENCE     (ksystm(55),iprec)
 DATA    NAME  / 4HEQMC,2HKM /
 DATA    kon   / 1,20,0,-1,0,0,0,0,1,8 /
 DATA    kg    / 1HG   /
 
!                      KNG
!     PARTITION  KGG = ---- , ONLY KMG SAVED
!                      KMG
 anyout =.false.
 nzz    = korsz (zz(1))
 lcor   = nzz
 luset  = kuset
 IF (kmpc == 0) GO TO 10
 kmpc   = -1
 CALL calcv (kscr(1),ug,un,um,zz)
 CALL ssg2a (kgg,0,kscr(2),kscr(1))
 
!     UNAPPEND FILES
 
 10 CONTINUE
 nzz3 = nzz - 3*isbz + 1
 nzz2 = nzz3 + isbz
 nzz1 = nzz2 + isbz
 nzz4 = nzz3
 IF (nskip <= 0) nzz4 = nzz3 - isbz
 
 IF (nskip <= 1) GO TO 30
 
 IF (kload <= 0) GO TO 20
 trl(1) = kpgg
 mcb(1) = kscr(4)
 CALL curcas (*15,nskip,trl,mcb,zz,nzz2)
 kpgg   = mcb(1)
 GO TO 20
 15 kload  = 0
 
 20 IF (kugv < 0) GO TO 30
 trl(1) = kugv
 mcb(1) = kscr(3)
 CALL curcas (*345,nskip,trl,mcb,zz,nzz2)
 kugv   = mcb(1)
 30 CONTINUE
 IF (kmpc == 0) GO TO 360
 
!                      PN
!     PARTITION  PGG = --- , ONLY PM SAVED
!                      PM
 IF (kload > 0) CALL ssg2a (kpgg,0,kscr(7),kscr(1))
 
!                 M
!     MULTIPLY  QM = -PM + KMG*UGV
 
 md(1) = kscr(5)
 mc(1) = kscr(7)
 CALL rdtrl (mc)
 IF (kload <= 0) mc(1) = 0
 ma(1) = kscr(2)
 mb(1) = kugv
 CALL rdtrl (ma)
 CALL rdtrl (mb)
 md(3) = ma(3)
 md(4) = mb(4)
!WKBR 11/93 SPR93007      MD(5) = 1
 md(5) = iprec
 mz    = nzz
 mt    = 0
 msab  = 1
 msc   =-1
!WKBR 11/93 SPR93007      MPR   = 1
 mpr   = iprec
 mscr  = kscr(1)
 CALL mpyad (zz,zz,zz)
 IF (md(3) == md(2)) md(4) = 1
 CALL wrttrl (md)
 
!                 N       T   M
!     MULTIPLY  QM = - GM * QM
 
 md(1) = kscr(6)
 mc(1) = 0
 ma(1) = kgm
 mb(1) = kscr(5)
 CALL rdtrl (ma)
 CALL rdtrl (mb)
 md(3) = ma(2)
 md(4) = mb(4)
!WKBR SPR93007      MD(5) = 1
 md(5) = iprec
 mt    = 1
 msab  =-1
 CALL mpyad (zz,zz,zz)
 IF (md(3) == md(2)) md(4) = 1
 CALL wrttrl (md)
 
!             N      M
!     MERGE QM AND QM ON SCRATCH 3
 
 trl(1) = kscr(5)
 CALL rdtrl (trl)
 IF (trl(1) < 0) GO TO 345
 kmpc   = 1
 
 CALL sdr1b (kscr(1),kscr(6),kscr(5),kscr(3),ug,un,um,kuset,0,0)
 trl(1) = kscr(3)
 CALL rdtrl (trl)
 msze   = 2*trl(3)
 
!     CREATE MPC-CONSTRAINT OUTPUT FILE
 
 IF (koqm <= 0) GO TO 360
 iapp = 10
 IF (nskip < 0) iapp = 20
 CALL mxcid (*345,zz,kg,msze/2,2,kuset,kgpl,ksil,nzz2)
 DO  i = 2,msze,2
   zz(i) = i/2
 END DO
 
!     SORT ON EXTERNAL ID
 
 nent = 2
 IF (msze == 2) GO TO 70
 
 ifil = ksil
 DO  i = 3,msze,2
   IF (zz(i) > zz(i-2)) GO TO 50
   trl(1) = zz(i)
   trl(2) = zz(i+1)
   CALL bishel (*410,trl,nent,2,zz(1))
   CYCLE
   50 nent   = nent + 2
 END DO
 
 70 CONTINUE
 lkscc  = msze + 1
 IF (lkscc+146 >= nzz4) GO TO 420
 ivec   = 0
 trl(1) = kscr(3)
 CALL rdtrl (trl)
 nvec   = trl(2)
 itypu  = 1
 incu   = 1
 isdone = 0
 lascas = .false.
 CALL gopen (kscc,zz(nzz1),rdrw)
 CALL gopen (kscr(3),zz(nzz3),rdrw)
 
 IF (nskip > 0) GO TO 90
 
!     POSITION LAMA
 
 ifil   = klam
 CALL OPEN (*350,klam,zz(nzz4),rdrw)
 CALL READ (*440,*450,klam,0,0,1,i)
 CALL READ (*440,*450,klam,0,0,1,i)
 90 ifil   = koqm
 CALL OPEN (*430,koqm,zz(nzz2),wrtrw)
 CALL fname (koqm,trl(1))
 trl(3) = itim
 trl(4) = idat(1)
 trl(5) = idat(2)
 trl(6) = idat(3)
 trl(7) = 1
 CALL WRITE (koqm,trl(1),7,1)
 
!     POSITION CASECC.  ASSUME USER WILL MISSET NSKIP
 
 IF (nskip <= 1) GO TO 100
 j    = nskip - 1
 ifil = kscc
 DO  i = 1,j
   CALL fwdrec (*440,kscc)
 END DO
 
!     LOOP ON EACH VECTOR
 
 100 ivec = ivec + 1
 
!     SUBCASE ID
 
 CALL READ (*160,*160,kscc,zz(lkscc),38,0,i)
 isb = zz(lkscc  )
 ild = zz(lkscc+3)
 ieg = 0
 
!     CLEAN UP UNUSED WORDS
 
 i = lkscc + 10
 j = lkscc + 49
 DO  k = i,j
   zz(k) = 0
 END DO
 
!     TITLES
 
 CALL fread (kscc,zz(lkscc+50),96,0)
 CALL fread (kscc,0,-31,0)
 CALL fread (kscc,lcc,1,0)
 CALL fread (kscc,0,-6,0)
 
!     MPCFORCE REQUEST
 
 ngset = 0
 lsetd = lkscc + 146
 CALL fread (kscc,ins(1),3,0)
 IF (ins(1) < 0) THEN
   GO TO   110
 ELSE IF (ins(1) == 0) THEN
   GO TO   120
 ELSE
   GO TO   130
 END IF
 
!     ALL REQUESTED
 
 110 CALL fread (kscc,0,0,1)
 GO TO 180
 
!     NONE REQUESTED
 
 120 ifil = klam
 CALL fread (kscc,0,0,1)
 IF (nskip > 0) GO TO 240
 CALL READ (*350,*350,klam,trl(1),7,0,i)
 GO TO 240
 
!     SET REQUESTED
 
 130 CONTINUE
 CALL fread (kscc,0,-lcc+176,0)
 
!     SKIP SYMMETRY SEQUENCE
 
 CALL fread (kscc,i,1,0)
 IF (i <= 0) GO TO 140
 CALL fread (kscc,0,-i,0)
 140 ifil = kscc
 CALL READ (*440,*450,kscc,trl(1),2,0,i)
 IF (trl(1) == ins(1)) GO TO 150
 CALL fread (kscc,0,-trl(2),0)
 GO TO 140
 150 IF (lsetd+trl(2) > nzz4) GO TO 420
 ngset = trl(2)
 CALL fread (kscc,zz(lsetd),ngset,1)
 GO TO 180
 
!     EOF ON CASE CONTROL.  CHECK IF REALLY DONE
 
 160 CONTINUE
 IF (nskip <   0) GO TO 170
 IF (ivec > nvec) GO TO 350
 ifil = kscc
 GO TO 440
 
 170 IF (ivec > nvec) GO TO 350
 lascas = .true.
 ivec   = ivec - 1
 
!     INITIALIZE
 
 180 CONTINUE
 IF (lascas) ivec = ivec + 1
 firstc = .true.
 firsto = .true.
 IF (nskip > 0) GO TO 190
 CALL READ (*235,*235,klam,trl(1),7,0,i)
 ild = trl(1)
 ieg = trl(3)
 190 mt  = lkscc - 1
 DO  j = 1,10
   i   = j + mt
   zz(i) = kon(j)
 END DO
 
 IF (ins(3) == 1) GO TO 210
 CALL page2 (2)
 WRITE  (6,205) uwm,NAME
 205 FORMAT (a25,' 2373, ONLY SORT1-REAL SUPPORTED IN ',2A4)
 210 zz(mt+1) = ins(2) + iapp
 zz(mt+4) = isb
 zz(mt+5) = ild
 zz(mt+6) = ieg
 lvec = lsetd + ngset - 1
 
!     LOOP ON POINT DATA
!   MT = POINTER TO MATCID GRID ID, MS = POINTER TO CASECC GRID REQUEST.
 
 idg = -1
 mt  = -1
 ms  = lsetd
 220 mt  = mt + 2
 IF (mt > msze) GO TO 240
 IF (zz(mt)/10 == idg) GO TO 220
 idg = zz(mt)/10
 IF (ins(1) < 0) GO TO 300
 
!     LOCATE POINT IN SET
 
 221 i = zz(ms)
 222 IF (ms-lvec < 0) THEN
   GO TO   223
 ELSE IF (ms-lvec == 0) THEN
   GO TO   228
 ELSE
   GO TO   240
 END IF
 223 IF (idg-i < 0) THEN
   GO TO   230
 ELSE IF (idg-i == 0) THEN
   GO TO   300
 END IF
 224 i = zz(ms+1)
 IF (i < 0) THEN
   GO TO   225
 ELSE
   GO TO   227
 END IF
 225 IF (idg+i > 0) THEN
   GO TO   226
 ELSE
   GO TO   300
 END IF
 226 ms = ms+2
 GO TO 221
 227 ms = ms+1
 GO TO 222
 
!     LAST POINT IN SET
 
 228 IF (i < 0 .AND. idg+i <= 0) GO TO 300
 IF (idg-i == 0) THEN
   GO TO   300
 END IF
 
!     NOT IN SET
 
 230 IF (mt+2 < lkscc) GO TO 220
 GO TO 240
 
!     END-OF-FILE
 
 235 isdone = 1
 
!     NO MORE GRIDS IN THIS SET
 
 240 CONTINUE
 IF (.NOT.firsto) CALL WRITE (koqm,0,0,1)
 IF (ivec+1 > nvec) GO TO 350
 IF (isdone /=    0) GO TO 350
 
!     CHECK IF COLUMN NEEDS TO BE SKIPPED
 
 ifil = kscr(3)
 IF (firstc) CALL fwdrec (*440,kscr(3))
 
 IF (lascas) GO TO 180
 GO TO 100
 
!     PROCESS THE GRID FOR OUTPUT
 
 300 mcb(1) = 10*idg + ins(2)
 IF (.NOT.firstc) GO TO 310
 IF (lvec+msze/2 > nzz4) GO TO 420
 inru = 1
 ilru = msze/2
 firstc = .false.
 CALL unpack (*240,kscr(3),zz(lvec+1))
 
 310 CONTINUE
 l      = 1
 nent   = 0
 ocb(3) = 0
 ocb(4) = 0
 ocb(5) = 0
 ocb(6) = 0
 ocb(7) = 0
 ocb(8) = 0
 m = MIN0(mt+10,msze)
 DO  i = mt,m,2
   j = zz(i)/10
   IF (j /= idg) GO TO 335
   k = zz(i+1) + lvec
   j = zz(i) - j *10 + 2
   IF (j > 2) GO TO 320
   
!     SCALAR
   
   l = 2
   j = 3
   
   320 ocb(j) = zz(k)
   IF (zz(k) /= 0) nent = nent + 1
 END DO
 
 335 ocb(2) = l
 IF (nent == 0) GO TO 220
 IF (.NOT.firsto) GO TO 340
 
!     WRITE OUT CONTROL RECORD (ODD NUMBER)
 
 anyout = .true.
 CALL WRITE (koqm,zz(lkscc),146,1)
 firsto = .false.
 
!     WRITE AN ENTRY OUT
 
 340 CALL WRITE (koqm,ocb(1),8,0)
 GO TO 220
 
!     CLOSE FILES
 
 345 CALL CLOSE (koqm,krw)
 koqm = -1
 anyout = .false.
 350 CALL CLOSE (kscr(3),krw)
 CALL CLOSE (kscc,krw)
 CALL CLOSE (klam,krw)
 IF (anyout) CALL eof (koqm)
 CALL CLOSE (koqm,krw)
 IF (.NOT.anyout) GO TO 360
 trl(1) = koqm
 trl(2) = nvec
 trl(3) = msze/2
 trl(4) = 0
 trl(5) = 0
 trl(6) = 1
 trl(7) = 0
 CALL wrttrl (trl)
 360 CONTINUE
 
!     CALCULATE UPDATED QG FILE - SCRATCH 5
 
 IF (kspc  == 0) GO TO 405
 IF (nskip <= 1) GO TO 405
 trl(1) = kqg
 mcb(1) = kscr(5)
 CALL curcas (*407,nskip,trl,mcb,zz,nzz2)
 kqg = mcb(1)
 
 405 RETURN
 
!     ERROR CONDITIONS
 
!     KQG BAD
 
 407 kspc = -1
 GO TO 405
 
 410 i = 7
 GO TO 490
 420 i = 8
 ifil = nzz4
 GO TO 490
 430 i = 1
 GO TO 490
 440 i = 2
 GO TO 490
 450 i = 3
 490 CALL mesage (i,ifil,NAME)
 
!     MPC OUTPUT FILE NOT CREATED, BUT DATA IS ON SCR3 FOR EQMCKS
 
 CALL page2 (2)
 WRITE  (nout,510) uwm,NAME
 510 FORMAT (a25,' 2380, MULTI-POINT CONSTRAINT FORCES NOT OUTPUT IN ',  &
     a4,a2,', SEE QUEUED MESSAGES.')
 GO TO 345
END SUBROUTINE eqmckm
