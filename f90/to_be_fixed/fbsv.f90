SUBROUTINE fbsv (*,zs,zd)
     
!     GIVEN A LOWER UNIT TRIANGULAR FACTOR WITH DIAGONAL SUPERIMPOSED
!     AND WRITTEN WITH TRAILING STRING DEFINITION WORDS, FBSV WILL
!     PERFORM THE FORWARD-BACKWARD SUBSTITUTION NECESSARY TO SOLVE A
!     LINEAR SYSTEM OF EQUATIONS.
 
!     THIS FBSV.MDS ROUTINE IS ALMOST SAME AS FBS.MIS
!     IT IS INTENDED TO BE USED FOR VAX, A VIRTUAL MEMORY MACHINE.
!     FBSV.MDS DIFFERS FROM FBS.MIS IN
!       1. THREE BUFFERS ARE USED
!       2. OPEN CORE IS REDUCED. THE INTENTION HERE IS TO AVOID
!          ACCESSIVE SYSTEM PAGING.
 
!     FBSV IS CALLED ONLY BY FBS ROUTINE.
!     IT IS INTRODUCED INTO NASTRAN REPERTORY BY G.CHAN/UNISYS, 10/88
 
 
 , INTENT(OUT)                            :: *
 REAL, INTENT(OUT)                        :: zs(1)
 DOUBLE  PRECISION, INTENT(OUT)           :: zd(1)
 LOGICAL :: ident
 INTEGER :: dbl     ,dbu      ,dbb      ,dbx      ,prec     ,SIGN    ,  &
     sysbuf  ,prc      ,words    ,rlcmpx   ,buf1     ,buf2    ,  &
     buf3    ,typel    ,typeb    ,typex    ,subnam(2),begn    ,  &
END     ,buf(2)   ,rc       ,eol      ,rd       ,rdrew   ,  &
    wrt     ,wrtrew   ,rew      ,eofnrw   ,switch   ,rsp     ,  &
    rdp     ,csp      ,cdp      ,BLOCK(15),hicore   ,sys34
REAL :: xs(4)    ,ys(4)
DOUBLE  PRECISION         xd       ,yd
COMMON /fbsx  / dbl(7)    ,dbu(7)   ,dbb(7)   ,dbx(7)   ,lcore   ,  &
    prec      ,SIGN
COMMON /system/ sysbuf    ,nout     ,sy1(28)  ,hicore   ,sy2(2)  , sys34
COMMON /names / rd        ,rdrew    ,wrt      ,wrtrew   ,rew     ,  &
    norew     ,eofnrw   ,rsp      ,rdp      ,csp     , cdp
COMMON /TYPE  / prc(2)    ,words(4) ,rlcmpx(4)
COMMON /packx / itype1    ,itype2   ,i1       ,j1       ,incr1
COMMON /unpakx/ itype3    ,i2       ,j2       ,incr2
COMMON /zntpkx/ xd(2)     ,ix       ,eol
COMMON /zblpkx/ yd(2)     ,iy
EQUIVALENCE    (dbl(5),typel) , (dbb(5),typeb) , (dbx(5),typex)  ,  &
    (xd(1) ,xs(1)) , (dbl(2),nl   ) , (yd(1) ,ys(1))
DATA    subnam/ 4HFBSV, 1H  / , begn/ 4HBEGN / , END/ 3HEND /


!     CHECK OPEN CORE SITUATION. IF SYS34 (PAGECNTL) IS NOT ZERO, SET
!     OPEN CORE SIZE TO SIZE SPECIFIED BY SYS34 (MINUS OVERHEAD OF 6000
!     WORDS). TOO BIG A CORE SIZE MAY CAUSE EXCESSIVE PAGING AND SLOW
!     DOWN FBS OPERATION


korchg = 0
IF (sys34 == 0) GO TO 50
korchg = hicore - sys34 + 6000
IF (korchg <= 0) GO TO 80
lcore  = lcore- korchg
WRITE (nout,40) sys34
40   FORMAT ('0*** SYSTEM INFORMATION MESSAGE - OPEN CORE FOR FBS IS',  &
    ' SET TO',i7,' WORDS BY PAGECNTL OF /SYSTEM/',/)
GO TO 80
50   IF (hicore > 130001) WRITE (nout,60) hicore
60   FORMAT ('0*** SYSTEM INFORMATION MESSAGE - PRESENT OPEN CORE =',  &
    i7, /5X,'FURTHER INCREASE OF OPEN CORE MAY ACTUALLY ',  &
    'SLOW DOWN FBS''S OPERATION.', //5X,  &
    'SUGGESTION: TO OPTIMIZE CORE USAGE AND USER''S PROBLEM,',  &
    ' AND TO MINIMIZE VAX''S PAGE FAULTS,', /5X,  &
    'CHECK WORKING_SET PAGE LIMIT ASSIGNED TO USER, AND SET ',  &
    'NASTRAN PAGECNTL WORD OF /SYSTEM/ TO MATCH, BUT NOT TO ',  &
    'EXCEED, THE CURRENT SETTING',/)

!     GENERAL INITIALIZATION

80   buf3   = lcore- sysbuf
buf2   = buf3 - sysbuf
buf1   = buf2 - sysbuf
nnn    = buf1 - 1
buf(1) = subnam(1)
buf(2) = begn
CALL conmsg (buf,2,0)
nbrlod = dbb(2)
rc     = rlcmpx(typeb)
i2     = 1
j2     = nl
incr2  = 1
i1     = 1
j1     = nl
incr1  = 1
itype1 = typel
itype2 = typex
itype3 = SIGN*typel
nwds   = words(typel)*nl
nvecs  = nnn/nwds
IF (nvecs == 0) CALL mesage(-8,nwds-nnn,subnam)
switch = 1
IF (typel == rsp .AND. rc == 2) switch = 2
IF (typel == rdp .AND. rc == 2) switch = 3
IF (switch /= 1) nvecs = nvecs/2
k1     = 1
BLOCK(1) = dbl(1)
dbx(2) = 0
dbx(6) = 0
dbx(7) = 0
ident  = .false.
IF (dbb(4) == 8) ident = .true.
nnndbl = nnn/2
nterms = rlcmpx(typel)*nl
IF (ident) nbrlod = nl

!     OPEN OUTPUT FILE (DBX), LOAD VECTORS FILE (DBB), AND LOWER
!     TRIANGULAR FACTOR FILE (DBL)

CALL gopen (dbx,zs(buf3),wrtrew)
CALL gopen (dbl,zs(buf1),rdrew )
IF (.NOT.ident) CALL gopen (dbb,zs(buf2),rdrew)

!     CHECK TIMING AND ISSUE MESSAGE

npass = (nbrlod+nvecs-1)/nvecs
CALL sswtch (11,l11)
IF (npass >= 10) l11=1
IF (l11 /= 1) GO TO 140
CALL page2 (-4)
WRITE (nout,100) typel,npass
100 FORMAT ('0*** USER INFORMATION MESSAGE FROM FBS',i1,' - NO. OF ',  &
    'PASSES NEEDED TO COMPLETE FBS OPERATION =',i5)
IF (npass > 15) WRITE (nout,110)
110 FORMAT (5X,'INCREASE OF OPEN CORE MAY ACTUALLY SLOW DOWN FBS ',  &
    'OPERATION.')
GO TO 140
120 IF (l11 < 0) GO TO 150
CALL cputim (j,t2,1)
t2 = t2-t1
IF (l11 > 0) WRITE (nout,130) t2
130 FORMAT (5X,'TIME TO COMPLETE ONE PASS =',f10.4,' CPU SECONDS',//)
l11 = -1
CALL tmtogo (j)
i = npass*t2
IF (j < i) CALL mesage (-50,i,subnam)
GO TO 150
140 CALL cputim (j,t1,1)

!     COMPUTE EXTENT OF THIS PASS

150 kn   = MIN0(k1+nvecs-1,nbrlod)
last = 1 + (kn-k1)*nwds
IF (ident) GO TO 190
SELECT CASE ( switch )
  CASE (    1)
    GO TO 160
  CASE (    2)
    GO TO 170
  CASE (    3)
    GO TO 180
END SELECT

!     NORMAL CASE - FILL CORE WITH LOAD VECTORS

160 DO  l=1,last,nwds
  CALL unpack (*162,dbb,zs(l))
  CYCLE
  162 ln = l+nwds-1
  DO  ll=l,ln
    zs(ll) = 0.
  END DO
END DO
GO TO 200

!     SPECIAL CASE - FACTOR IS RSP AND VECTORS ARE CSP

170 last = 1 + 2*(kn-k1)*nwds + nwds
l = 0
DO  k=1,nnndbl
  zd(k) = 0.
END DO
DO  k=k1,kn
  icspsg = csp*SIGN
  CALL intpk (*176,dbb,0,icspsg,0)
  172 CALL zntpki
  zs(l+ix   ) = xs(1)
  zs(l+ix+nl) = xs(2)
  IF (eol == 0) GO TO 172
  176 l = l + 2*nl
END DO
GO TO 200

!     SPECIAL CASE - FACTOR IS RDP AND VECTORS ARE CDP

180 last = 1 + 2*(kn-k1)*nwds + nwds
l = 0
DO  k=1,nnndbl
  zd(k) = 0.
END DO
DO  k=k1,kn
  icdpsg = cdp*SIGN
  CALL intpk (*186,dbb,0,icdpsg,0)
  182 CALL zntpki
  zd(l+ix   ) = xd(1)
  zd(l+ix+nl) = xd(2)
  IF (eol == 0) GO TO 182
  186 l = l + 2*nl
END DO
GO TO 200

!     SPECIAL CASE - GENERATE IDENTITY MATRIX

190 l = 0
DO  k=1,nnndbl
  zd(k) = 0.
END DO
DO  k=k1,kn
  SELECT CASE ( typel )
    CASE (    1)
      GO TO 191
    CASE (    2)
      GO TO 192
    CASE (    3)
      GO TO 193
    CASE (    4)
      GO TO 194
  END SELECT
  191 zs(l+k) = 1.0
  GO TO 196
  192 zd(l+k) = 1.0D0
  GO TO 196
  193 zs(l+2*k-1) = 1.0
  GO TO 196
  194 zd(l+2*k-1) = 1.0D0
  196 l = l + nterms
END DO

!     COMPUTE FORWARD-BACKWARD SUBSTITUTION ON LOAD VECTORS NOW IN CORE

200 CALL REWIND (dbl)
CALL fwdrec (*270,dbl)

SELECT CASE ( typel )
  CASE (    1)
    GO TO 201
  CASE (    2)
    GO TO 202
  CASE (    3)
    GO TO 203
  CASE (    4)
    GO TO 204
END SELECT

201 CALL fbs1 (BLOCK,zs,zs(last),nwds)
GO TO 210
202 CALL fbs2 (BLOCK,zs,zs(last),nwds)
GO TO 210
203 CALL fbs3 (BLOCK,zs,zs(last),nwds)
GO TO 210
204 CALL fbs4 (BLOCK,zs,zs(last),nwds)

!     PACK SOLUTION VECTORS ONTO OUTPUT FILE

210 SELECT CASE ( switch )
  CASE (    1)
    GO TO 220
  CASE (    2)
    GO TO 230
  CASE (    3)
    GO TO 240
END SELECT

!     NORMAL CASE - CALL PACK

220 DO  l=1,last,nwds
  CALL pack (zs(l),dbx,dbx)
END DO
GO TO 250

!     SPECIAL CASE - FACTOR IS RSP AND VECTORS ARE CSP, CALL BLDPK

230 l = 0
DO  k=k1,kn
  CALL bldpk (csp,typex,dbx,0,0)
  DO  i=1,nl
    ys(1) = zs(l+i   )
    ys(2) = zs(l+i+nl)
    iy = i
    CALL zblpki
  END DO
  CALL bldpkn (dbx,0,dbx)
  l = l + 2*nl
END DO
GO TO 250

!     SPECIAL CASE - FACTOR IS RDP AND VECTORS ARE CDP, CALL BLDPK

240 l = 0
DO  k=k1,kn
  CALL bldpk (cdp,typex,dbx,0,0)
  DO  i=1,nl
    yd(1) = zd(l+i   )
    yd(2) = zd(l+i+nl)
    iy = i
    CALL zblpki
  END DO
  CALL bldpkn (dbx,0,dbx)
  l = l + 2*nl
END DO

!     TEST FOR MORE PASSES

250 IF (kn == nbrlod) GO TO 300
k1   = kn + 1
GO TO 120

!     ERROR

270 CALL mesage (-2,dbl,subnam)

300 IF (.NOT.ident) CALL CLOSE (dbb,rew)
CALL CLOSE (dbl,rew)
CALL CLOSE (dbx,rew)
buf(2) = END
CALL conmsg (buf,2,0)
RETURN 1
END SUBROUTINE fbsv
