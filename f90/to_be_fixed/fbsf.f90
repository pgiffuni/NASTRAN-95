SUBROUTINE fbsf (zs,zd)
     
!     GIVEN A LOWER TRIANGULAR FACTOR WITH DIAGONAL SUPERIMPOSED, AND
!     WRITTEN WITH TRAILING STRING DEFINITION WORDS, FBS WILL PERFORM
!     THE FORWARD-BACKWARD SUBSTITUTION NECESSARY TO SOLVE A LINEAR
!     SYSTEM OF EQUATIONS.
 
!     OPEN CORE IS DEFINED AS FOLLOWS
 
!     ZS(   1         ) - FIRST RIGHT HAND VECTOR ON FILE DBB
!                         (SIZE = NCOL*NWDS)
!                         NCOL = NUMBER OF COLUMNS (ROWS) IN LOWER
!                                TRIANGULAR MATRIX
!                         NWDS = 1, IF MATRICES ARE REAL SINGLE
!                              = 2, IF MATRICES ARE REAL DOUBLE OR
!                                COMPLEX SINGLE
!                              = 4, IF MATRICES ARE COMPLEX DOUBLE
!     ZS( NCOL*NWDS+1 ) - NEXT RIGHT HAND VECTOR
!         .
!         .               ( "KN" RIGHT HAND VECTORS WILL BE LOADED INTO
!         .               MEMORY)
!         .
!     ZS( BUF1        ) - BUFFER FOR FILE WITH RIGHT HAND VECTORS
!                         AND FOR SOLUTION VECTORS
!     ZS( BUF2        ) - BUFFER FOR FILE WITH TRIANGULAR MATRIX
 
 
 REAL, INTENT(OUT)                        :: zs(1)
 DOUBLE  PRECISION, INTENT(OUT)           :: zd(1)
 IMPLICIT INTEGER (a-z)
 LOGICAL :: ident
INTEGER :: subnam(2) ,BLOCK(15),begn     ,END
REAL :: xs(4)    ,ys(4)
DOUBLE  PRECISION         xd       ,yd
CHARACTER (LEN=29) :: uim
CHARACTER (LEN=25) :: uwm
CHARACTER (LEN=23) :: ufm
COMMON /logout/ lout
COMMON /xmssg / ufm       ,uwm      ,uim
COMMON /fbsx  / dbl(7)    ,dbu(7)   ,dbb(7)   ,dbx(7)   ,lcore   ,  &
    prec      ,SIGN     ,scrx
COMMON /system/ sysbuf    ,nout     ,skip(91) ,ksys94
COMMON /names / rd        ,rdrew    ,wrt      ,wrtrew   ,rew     ,  &
    norew     ,eofnrw   ,rsp      ,rdp      ,csp     , cdp
COMMON /TYPE  / prc(2)    ,words(4) ,rlcmpx(4)
COMMON /packx / itype1    ,itype2   ,i1       ,j1       ,incr1
COMMON /unpakx/ itype3    ,i2       ,j2       ,incr2
COMMON /zntpkx/ xd(2)     ,ix       ,eol
COMMON /zblpkx/ yd(2)     ,iy
EQUIVALENCE     (dbl(2),nl),   (dbb(5),typeb), (dbx(5),typex),  &
    (xd(1),xs(1)), (yd(1),ys(1))
DATA    subnam/ 4HFBSF,4H    /
DATA    begn  / 4HBEGN/
DATA    END   / 4HEND /

!     GENERAL INITIALIZATION

buf2   = lcore - sysbuf
buf1   = buf2  - sysbuf
rc     = rlcmpx(typeb)
typel  = dbl(5)
wds    = words(typel)
nwds   = wds*nl
nbrlod = dbb(2)
ident  = .false.
IF (dbb(4) == 8) ident = .true.
IF (ident) nbrlod = nl
switch = 1
IF (typel == rsp .AND. rc == 2) switch = 2
IF (typel == rdp .AND. rc == 2) switch = 3
dbl1   = dbl(1)
nnn    = buf1 - 1
nvecs  = nnn/nwds
IF (nvecs == 0) CALL mesage (-8,nwds-nnn,subnam)
IF (switch /= 1) nvecs = nvecs/2
npass  = (nbrlod+nvecs-1)/nvecs
subnam(2) = begn
CALL conmsg (subnam,2,0)
40 npass  = (nbrlod+nvecs-1)/nvecs
IF ( npass == 1 ) GO TO 50
need = nwds*nbrlod + 2*sysbuf
WRITE ( lout, 9001 ) npass, need
9001  FORMAT(i4,' PASSES REQUIRED, OPEN CORE NEEDS TO BE ',i7  &
    ,' TO ELIMINATE THIS')
50 CONTINUE
i2     = 1
j2     = nl
incr2  = 1
i1     = 1
j1     = nl
incr1  = 1
itype1 = typel
itype2 = typex
itype3 = SIGN*typel
dbx(2) = 0
dbx(6) = 0
dbx(7) = 0
nnndbl = nnn/2
nterms = rlcmpx(typel)*nl
k1     = 1
oprd   = rdrew
opwrt  = wrtrew
BLOCK(1) = dbl(1)

!     OPEN LOWER TRIANGULAR FACTOR FILE (DBL1)

CALL gopen (dbl1,zs(buf2),rdrew)

!     OPEN RIGHT HAND VECTORS FILE (DBB) AND COMPUTE EXTENT OF THIS PASS

100 kn    = MIN0(k1+nvecs-1,nbrlod)
last  = (kn-k1+1)*nwds
opcls = norew
IF (kn == nbrlod) opcls = rew
IF (ident) GO TO 280
CALL gopen (dbb,zs(buf1),oprd)
SELECT CASE ( switch )
  CASE (    1)
    GO TO 140
  CASE (    2)
    GO TO 180
  CASE (    3)
    GO TO 230
END SELECT

!     NORMAL CASE - FILL CORE WITH RIGHT HAND VECTORS

140 DO  l = 1,last,nwds
  CALL unpack (*150,dbb,zs(l))
  CYCLE
  150 ln = l + nwds - 1
  DO  ll = l,ln
    zs(ll) = 0.
  END DO
END DO
GO TO 390

!     SPECIAL CASE - FACTOR IS RSP AND VECTORS ARE CSP

180 last = 2*(kn-k1+1)*nwds
l = 0
DO  k = 1,nnndbl
  zd(k) = 0.0D+0
END DO
DO  k = k1,kn
  icspsg = csp*SIGN
  CALL intpk (*210,dbb,0,icspsg,0)
  200 CALL zntpki
  zs(l+ix   ) = xs(1)
  zs(l+ix+nl) = xs(2)
  IF (eol == 0) GO TO 200
  210 l = l + 2*nl
END DO
GO TO 390

!     SPECIAL CASE - FACTOR IS RDP AND VECTORS ARE CDP

230 last = 2*(kn-k1+1)*nwds
l = 0
DO  k = 1,nnndbl
  zd(k) = 0.0D+0
END DO
DO  k = k1,kn
  icdpsg = cdp*SIGN
  CALL intpk (*260,dbb,0,icdpsg,0)
  250 CALL zntpki
  zd(l+ix   ) = xd(1)
  zd(l+ix+nl) = xd(2)
  IF (eol == 0) GO TO 250
  260 l = l + 2*nl
END DO
GO TO 390

!     SPECIAL CASE - GENERATE IDENTITY MATRIX

280 DO  k = 1,nnndbl
  zd(k) = 0.0D+0
END DO
l = 0
SELECT CASE ( typel )
  CASE (    1)
    GO TO 300
  CASE (    2)
    GO TO 320
  CASE (    3)
    GO TO 340
  CASE (    4)
    GO TO 360
END SELECT
300 DO  k = k1,kn
  zs(l+k) = 1.0
  l = l + nterms
END DO
GO TO 400
320 DO  k = k1,kn
  zd(l+k) = 1.0D+0
  l = l + nterms
END DO
GO TO 400
340 DO  k = k1,kn
  zs(l+2*k-1) = 1.0
  l = l + nterms
END DO
GO TO 400
360 DO  k = k1,kn
  zd(l+2*k-1) = 1.0D+0
  l = l + nterms
END DO
GO TO 400

!    CLOSE RIGHT HAND VECTORS FILE (DBB).
!    START FORWARD-BACKWARD SUBSTITUTION ON RIGHT HAND VECTORS NOW IN CORE

390 CALL CLOSE  (dbb,opcls)
400 CALL REWIND (dbl1)
CALL fwdrec (*610,dbl1)

j = typel
SELECT CASE ( j )
  CASE (    1)
    GO TO 410
  CASE (    2)
    GO TO 420
  CASE (    3)
    GO TO 430
  CASE (    4)
    GO TO 440
END SELECT
410 CALL fbs1 (BLOCK,zs,zs(last),nwds)
GO TO 500
420 CALL fbs2 (BLOCK,zs,zs(last),nwds)
GO TO 500
430 CALL fbs3 (BLOCK,zs,zs(last),nwds)
GO TO 500
440 CALL fbs4 (BLOCK,zs,zs(last),nwds)
GO TO 500

!     OPEN AND PACK SOLUTION VECTORS ONTO OUTPUT FILE (DBX)

500 CALL gopen (dbx,zs(buf1),opwrt)
SELECT CASE ( switch )
  CASE (    1)
    GO TO 510
  CASE (    2)
    GO TO 530
  CASE (    3)
    GO TO 560
END SELECT

!     NORMAL CASE - CALL PACK

510 DO  l = 1,last,nwds
  CALL pack (zs(l),dbx,dbx)
END DO
GO TO 600

!     SPECIAL CASE - FACTOR IS RSP AND VECTORS ARE CSP, CALL BLDPK

530 l = 0
DO  k = k1,kn
  CALL bldpk (csp,typex,dbx,0,0)
  DO  i = 1,nl
    ys(1) = zs(l+i   )
    ys(2) = zs(l+i+nl)
    iy = i
    CALL zblpki
  END DO
  CALL bldpkn (dbx,0,dbx)
  l = l + 2*nl
END DO
GO TO 600

!     SPECIAL CASE - FACTOR IS RDP AND VECTORS ARE CDP, CALL BLDPK

560 l = 0
DO  k = k1,kn
  CALL bldpk (cdp,typex,dbx,0,0)
  DO  i = 1,nl
    yd(1) = zd(l+i   )
    yd(2) = zd(l+i+nl)
    iy = i
    CALL zblpki
  END DO
  CALL bldpkn (dbx,0,dbx)
  l = l + 2*nl
END DO

!     CLOSE OUTPUT FILE, AND TEST FOR MORE PASSES

600 CALL CLOSE (dbx,opcls)
IF (kn == nbrlod) GO TO 620
k1   = kn + 1
oprd = rd
opwrt= wrt
GO TO 100

!     ERROR

610 CALL mesage (-2,dbl1,subnam)

!     JOB DONE. CLOSE TRIANGULAR FACTOR FILE.

620 CALL CLOSE (dbl1,rew)
subnam(2) = END
CALL conmsg (subnam,2,0)
RETURN
END SUBROUTINE fbsf
