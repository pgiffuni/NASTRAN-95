SUBROUTINE eqmcks
     
!     THIS SUBROUTINE CALCULATES AND OUTPUTS OVERALL EQUILIBRIUM FORCES
 
!     THE INPUT FILES ARE
!         KSCC - CASE CONTROL     - NOT PREPOSITIONED.
!         KPGG - LOAD VECTORS     - FILE 110 OR SCRATCH4
!         KQG  - SPC CONSTRAINTS  - FILE 111 OR SCRATCH5
!         QMG  - MPC CONSTRAINTS  - SCRATCH3
!         DT   - RIGID BODY TRANS - SCRATCH2
 
 LOGICAL :: lsteig
 INTEGER :: eject    ,NAME(2)  ,parm     , rdnrw    ,rdrw     ,wrtnrw   ,wrtrw
 REAL :: head(2,4),cor1(8,1),cor3(8,3)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /names / rdnrw,rdrw,wrtnrw,wrtrw,krw,knrw,knerw
!WKBR 3/94 SPR93007      COMMON /SYSTEM/ ISBZ,NOUT
 COMMON /system/ isbz,nout,dum(52),iprec
 COMMON /BLANK / iopt,igpt,nskip,skpb(15),core(8,4)
 COMMON /unpakx/ iunpr,iunrw,nunrw,iuninc
 COMMON /mpyadx/ ma(7),mb(7),mc(7),md(7),mz,mt,msab,msc,mpr,mscr
 COMMON /eqmk1 / kscc,keqin(8),kpgg,kqg,kcstm,klama,koqm,kscr(7)  &
     ,               kmpc,kload,kspc,parm(4)
!ZZ   COMMON /ZZEQMS/ ZZ(1)
 COMMON /zzzzzz/ zz(20000)
 EQUIVALENCE     (mb(6),freq), (core(1,1),cor1(1,1),cor3(1,1))
 DATA    NAME  / 4HEQMC,4HKS    /
 DATA    head  / 4HAPPL,4HIED , 4HSPCF,4HORCE,  4HMPCF,4HORCE  &
     ,               4H---t,4HOTAL  /
 
 parm(3) = NAME(1)
 parm(4) = NAME(2)
 nzz   = korsz (zz)
 nzz3  = nzz  - 3*isbz + 1
 nzz2  = nzz3 + isbz
 nzz1  = nzz2 + isbz
 
 nvec  = 0
 ma(1) = kscr(2)
 mc(1) = 0
 CALL rdtrl (ma)
 mz    = nzz
 mt    = 0
 msab  = 1
 msc   = 1
!WKBR 11/93 SPR93007      MPR   = 1
 mpr = iprec
 mscr  = kscr(1)
 
!     CALCULATE  DT*PG  ON SCRATCH7
 
 IF (kload <= 0) GO TO 40
 mb(1) = kpgg
 md(1) = kscr(7)
 CALL rdtrl (mb)
 md(3) = ma(3)
 md(4) = mb(4)
!WKBR 11/93 SPR93007      MD(5) = 1
 md(5) = iprec
 CALL mpyad (zz,zz,zz)
 IF (md(3) == md(2)) md(4) = 1
 CALL wrttrl (md)
 nvec  = md(2)
 
!     CALCULATE DT*QG  ON SCRATCH6
 
 40 IF (kspc <= 0) GO TO 50
 mb(1) = kqg
 md(1) = kscr(6)
 CALL rdtrl (mb)
 md(3) = ma(3)
 md(4) = mb(4)
!WKBR 11/93 SPR93007      MD(5) = 1
 md(5) = iprec
 CALL mpyad (zz,zz,zz)
 IF (md(3) == md(2)) md(4) = 1
 CALL wrttrl (md)
 nvec = MAX0(nvec,md(2))
 
!     CALCULATE  DT*MPC  ON SCRATCH5
 
 50 IF (kmpc <= 0) GO TO 60
 md(1) = kscr(5)
 mb(1) = kscr(3)
 CALL rdtrl (mb)
 parm(2) = mb(1)
 IF (mb(1) <= 0) GO TO 520
 md(3) = ma(3)
 md(4) = mb(4)
!WKBR 11/93 SPR93007      MD(5) = 1
 md(5) = iprec
 CALL mpyad (zz,zz,zz)
 IF (md(3) == md(2)) md(4) = 1
 CALL wrttrl (md)
 nvec  = MAX0(md(2),nvec)
 60 IF (nvec <= 0) GO TO 400
 
!     POSITION CASE CONTROL
 
 CALL gopen (kscc,zz(nzz1),rdrw)
 IF (nskip > 0) GO TO 70
 
!     RESERVE THIRD BUFFER FOR LAMA
 
 ibfl = nzz3
 parm(2) = klama
 CALL gopen (klama,zz(nzz3),rdrw)
 CALL fwdrec (*510,klama)
 GO TO 90
 70 ibfl = nzz2
 IF (nskip <= 1) GO TO 90
 
!     ASSUME USER MAY MALADJUST NSKIP
 
 j = nskip - 1
 parm(2) = kscc
 DO  i = 1,j
   CALL fwdrec (*510,kscc)
 END DO
 
!     READ INTO CORE AS MANY (MAXVEC) VECTORS THAT FIT
 
 90 nentry = 0
 IF (kload > 0) nentry = 6
 IF (kmpc  > 0) nentry = nentry + 6
 IF (kspc  > 0) nentry = nentry + 6
 
 maxvec = (ibfl-1)/nentry
 IF (maxvec >= nvec) GO TO 110
 
!     INSUFFICIENT CORE TO DO ALL VECTORS
 
 CALL page2 (2)
 WRITE  (nout,100) uwm,maxvec,NAME
 100 FORMAT (a25,' 2374, INSUFFICIENT CORE TO PROCESS MORE THAN',i7,  &
     ' VECTORS IN ',2A4)
 
 IF (maxvec <= 0) GO TO 400
 
 110 maxvec = MIN0 (nvec,maxvec)
 l = 1
 ma(1) = 0
 IF (kload <= 0) GO TO 160
 parm(2) = kscr(7)
 ma(1) = 1
 ASSIGN 160 TO iret
 
!     INTERNAL FUNCTION TO LOAD MAXVEC COLUMNS INTO CORE
 
 120 CONTINUE
 CALL gopen (parm(2),zz(nzz2),rdrw)
 iunpr  = 1
 iuninc = 1
 iunrw  = 1
 nunrw  = 6
 
 DO  mt = 1,maxvec
   CALL unpack (*130,parm(2),zz(l))
   GO TO 150
   130 mpr = l - 1
   DO  i = 1,6
     mpr = mpr + 1
     zz(mpr) = 0.0
   END DO
   150 l = l + 6
 END DO
 
 CALL CLOSE (parm(2),krw)
 GO TO iret, (160,170,180)
 
 160 ma(2) = 0
 IF (kspc <= 0) GO TO 170
 parm(2) = kscr(6)
 ma(2) = l
 ASSIGN 170 TO iret
 GO TO 120
 
 170 ma(3) = 0
 IF (kmpc <= 0) GO TO 180
 parm(2) = kscr(5)
 ma(3) = l
 ASSIGN 180 TO iret
 GO TO 120
 
 180 ivec = 0
 lsteig = .false.
 CALL page1
 
!     LOOP ON OUTPUT
 
 200 CONTINUE
 ivec = ivec + 1
 IF (lsteig) GO TO 260
 parm(2) = kscc
 CALL READ (*250,*500,kscc,mb(1),7,1,i)
 i = mb(1)
 IF (ivec == 1 .OR. eject(11) /= 0) WRITE (nout,210) igpt
 210 FORMAT (1H0,20X,'E Q U I L I B R I U M   C H E C K   L O A D S',  &
     /,1H0,16X,'RESULTANT LOADS AT POINT',i7, ' IN BASIC COORDINATE SYSTEM')
 IF (nskip <= 0) GO TO 260
 
!     STATICS SUBCASES
 
 IF (mb(4) == 0) mb(4) = mb(7)
 IF (mb(4) == 0) mb(4) = mb(6)
 WRITE  (nout,220) mb(1),mb(4)
 220 FORMAT (1H0,24X,7HSUBCASE,i8,8H,   load,i8)
 WRITE  (nout,230)
 230 FORMAT (1H0,5X,46H-TYPE-        t1             t2             t3,  &
     13X,32HR1             r2             r3)
 
 240 FORMAT (5X,2A4,1P,6E15.6)
 
 GO TO 300
 
!     EOF FOUND
 
 250 CONTINUE
 IF (ivec > maxvec) GO TO 400
 IF (nskip > 0) GO TO 510
 lsteig = .true.
 
!     EIGENVALUE PROBLEM
 
 260 parm(2) = klama
 CALL READ (*510,*500,klama,mb(2),7,0,i)
 WRITE  (nout,270) mb(1),mb(2),freq
 270 FORMAT (1H0,24X,7HSUBCASE,i8,8H,   mode,i5,13H,   frequency, 1P,e15.6)
 WRITE  (nout,230)
 
!     LOOP ON OUTPUT CATAGORY
 
 300 k = nentry/6 + 1
 ihdcnt = 1
 DO  i = 3,8
   core(i,k) = 0.0E0
 END DO
 
 DO  i = 1,3
   IF (ma(i) == 0) CYCLE
   core(1,ihdcnt) = head(1,i)
   core(2,ihdcnt) = head(2,i)
   j = ma(i) + ivec*6 - 6
   
   DO  l = 3,8
     core(l,ihdcnt) = zz(j)
     core(l,k) = core(l,k) + zz(j)
     j = j + 1
   END DO
   ihdcnt = ihdcnt + 1
 END DO
 
 core(1,k) = head(1,4)
 core(2,k) = head(2,4)
 IF (k == 2) WRITE (nout,240) cor1
 IF (k == 3) WRITE (nout,240) cor3
 IF (k == 4) WRITE (nout,240) core
 IF (ivec < maxvec) GO TO 200
 400 CALL CLOSE (kscc,krw)
 IF (nskip <= 0) CALL CLOSE (klama,krw)
 RETURN
 
!     ERROR MESSAGES
 
!     EOR
 
 500 parm(1) = 3
 GO TO 600
 
!     EOF
 
 510 parm(1) = 2
 GO TO 600
 
!     ILLEGAL INPUT
 
 520 parm(1) = 1
 GO TO 600
 
 600 CALL mesage (parm(1),parm(2),parm(3))
 GO TO 400
END SUBROUTINE eqmcks
