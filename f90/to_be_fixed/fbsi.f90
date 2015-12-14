SUBROUTINE fbsi (zs,zd)
     
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
!         .               ( "NRHV" RIGHT HAND VECTORS WILL BE LOADED INTO
!         .               MEMORY)
!         .
!     ZS( MTRIA       ) - MEMORY FOR STORAGE OF ALL OR PART OF THE LOWER
!                         TRIANGULAR MATRIX.  (SEE SUBROUTINE FBSRDM FOR
!                         FORMAT OF STORAGE OF MATRIX.)
!     ZS( BUF1        ) - BUFFER FOR FILE WITH RIGHT HAND VECTORS
!                         AND FOR SOLUTION VECTORS
!     ZS( BUF2        ) - BUFFER FOR FILE WITH TRIANGULAR MATRIX
 
 
 REAL, INTENT(OUT)                        :: zs(1)
 DOUBLE  PRECISION, INTENT(OUT)           :: zd(1)
 IMPLICIT INTEGER (a-z)
 LOGICAL :: ident
INTEGER :: subnam(2) ,BLOCK(15),begn     ,END      ,iname(2)
REAL :: xs(4)    ,ys(4)
DOUBLE  PRECISION         xd       ,yd
CHARACTER (LEN=29) :: uim
CHARACTER (LEN=25) :: uwm
CHARACTER (LEN=23) :: ufm
COMMON /logout/ lout
COMMON /xmssg / ufm       ,uwm      ,uim
COMMON /fbsx  / dbl(7)    ,dbu(7)   ,dbb(7)   ,dbx(7)   ,lcore   ,  &
    prec      ,SIGN     ,scrx
COMMON /fbsm  / nvec      ,nvecsz   ,nwds     ,lasind   ,ipos(7)
COMMON /system/ sysbuf    ,nout     ,skip(91) ,ksys94
COMMON /names / rd        ,rdrew    ,wrt      ,wrtrew   ,rew     ,  &
    norew     ,eofnrw   ,rsp      ,rdp      ,csp     , cdp
COMMON /TYPE  / prc(2)    ,words(4) ,rlcmpx(4)
COMMON /packx / itype1    ,itype2   ,i1       ,j1       ,incr1
COMMON /unpakx/ itype3    ,i2       ,j2       ,incr2
COMMON /zntpkx/ xd(2)     ,ix       ,eol
COMMON /zblpkx/ yd(2)     ,iy
EQUIVALENCE     (dbl(2),ncol),   (dbb(5),typeb), (dbx(5),typex),  &
    (xd(1),xs(1)), (yd(1),ys(1))
DATA    subnam/ 4HFBSI,4H    /
DATA    begn  / 4HBEGN/
DATA    END   / 4HEND /

!     GENERAL INITIALIZATION

buf2     = lcore - sysbuf
buf1     = buf2  - sysbuf
typel    = dbl(5)
rcb      = rlcmpx( typeb )
rcl      = rlcmpx( typel )
nwds     = words ( typeb )
IF ( rcb == rcl .AND. typel > typeb ) nwds = words( typel )
nrhvwd   = nwds * ncol
nwds     = words ( typel )
nrhv     = dbb(2)
ident    = .false.
IF (dbb(4) == 8) ident = .true.
IF (ident) nrhv = ncol
switch   = 1

! SET SWITCH AS FOLLOWS:
!  =1, IF LOWER TRIANGULAR MATRIX AND RIGHT HAND VECTORS ARE SAME TYPE
!  =2, LOWER TRIANGULAR MATRIX IS REAL SINGLE AND RIGHT HAND VECTOR IS
!      COMPLEX
!  =3, LOWER TRIANGULAR MATRIX IS REAL DOUBLE AND RIGHT HAND VECTOR IS
!      COMPLEX
!  (NOTE, IF SWITCH IS .NE. 1, THEN THE REAL AND IMAGINARY PARTS OF THE
!   THE RIGHT HAND VECTOR ARE TREATED AS TWO SEPARATE VECTORS.  I.E.,
!   THE REAL PART BECOMES ONE VECTOR AND THE IMAGINARY PART BECOMES A
!   SECOND VECTOR.)

IF (typel == rsp .AND. rcb == 2) switch = 2
IF (typel == rdp .AND. rcb == 2) switch = 3
IF (switch == 1 ) GO TO 90
IF (switch == 3 ) GO TO 70
nrhvwd   = 2 * ncol
GO TO 90
70    CONTINUE
nrhvwd   = 4 * ncol
90    CONTINUE
mtria    = nrhv * nrhvwd + 1

! ENSURE DOUBLE WORD BOUNDARY

mtria    = ( mtria/2 ) * 2 + 1
memavl   = buf1 - mtria - 2
subnam(2) = begn
CALL conmsg (subnam,2,0)
CALL fbsrdm ( dbl   , zs(mtria), zs(mtria), zs(mtria)  &
    ,             memavl, zs(buf2 ), lasind   , ipos )
CALL sswtch ( 47, l47 )
CALL fname ( dbl, iname )
IF ( l47 == 0 ) GO TO 100
WRITE ( lout, 9001 ) dbl(1), iname, ipos( 1 ), ncol, lcore, memavl
CALL fname ( dbb, iname )
WRITE ( lout, 9002 ) dbb(1), iname, dbl, dbb
9001  FORMAT(4X  &
    ,      ' FORWARD BACKWARD SUBSTITUTION OF FILE ',i3,'   NAME=',2A4  &
    ,/,4X, ' LAST COLUMN OF TRIANGULAR MATRIX IN MEMORY        =',i8  &
    ,/,4X, ' TOTAL COLUMNS IN TRIANGULAR MATRIX                =',i8  &
    ,/,4X, ' TOTAL OPEN CORE AVAILABLE FOR USE                 =',i8  &
    ,/,4X, ' OPEN CORE AVAILABLE FOR TRIANGULAR MATRIX STORAGE =',i8 )
9002  FORMAT(4X ,      ' RIGHT HAND VECTOR FILE ',i3,'   NAME=',2A4  &
    ,/,4X, ' TRIANGULAR MATRIX TRAILER    =', 7I6  &
    ,/,4X, ' RIGHT HAND VECTOR(S) TRAILER =', 7I6 )
100   CONTINUE
i2       = 1
j2       = ncol
incr2    = 1
i1       = 1
j1       = ncol
incr1    = 1
itype1   = typel
itype2   = typex
itype3   = SIGN*typel
dbx(2)   = 0
dbx(6)   = 0
dbx(7)   = 0
BLOCK(1) = dbl(1)

!     OPEN RIGHT HAND VECTORS FILE (DBB)

last  = nrhv*nrhvwd
IF ( ident ) GO TO 280
CALL gopen ( dbb, zs(buf1), rdrew )
SELECT CASE ( switch )
  CASE (    1)
    GO TO  140
  CASE (    2)
    GO TO  180
  CASE (    3)
    GO TO  230
END SELECT

!     READ RIGHT HAND VECTORS INTO MEMORY

140 DO  l = 1, last, nrhvwd
  CALL unpack ( *150, dbb, zs(l) )
  CYCLE
  150 ln = l + nrhvwd - 1
  DO  ll = l,ln
    zs( ll ) = 0.
  END DO
END DO
GO TO 390

!     SPECIAL CASE - LOWER TRIANGULAR MATRIX IS RSP AND VECTORS ARE CSP

180 last2 = last / 2
l = 0
DO  k = 1,last2
  zd(k)    = 0.0D+0
END DO
DO  k = 1, nrhv
  icspsg   = csp*SIGN
  CALL intpk ( *210, dbb, 0, icspsg, 0 )
  200 CALL zntpki
  zs(l+ix     ) = xs(1)
  zs(l+ix+ncol) = xs(2)
  IF ( eol == 0 ) GO TO 200
  210 l = l + 2*ncol
END DO
GO TO 390

!     SPECIAL CASE - LOWER TRIANGULAR MATRIX IS RDP AND VECTORS ARE CDP

230 last2 = last / 2
l = 0
DO  k = 1,last2
  zd(k)    = 0.0D+0
END DO
DO  k = 1, nrhv
  icdpsg   = cdp*SIGN
  CALL intpk ( *260, dbb, 0, icdpsg, 0 )
  250 CALL zntpki
  zd(l+ix     ) = xd(1)
  zd(l+ix+ncol) = xd(2)
  IF ( eol == 0 ) GO TO 250
  260 l        = l + 2*ncol
END DO
GO TO 390

!     SPECIAL CASE - GENERATE IDENTITY MATRIX

280 last = nrhv * nrhvwd
DO  k = 1,last
  zd(k) = 0.0D+0
END DO
l = 0
SELECT CASE ( typel )
  CASE (    1)
    GO TO  300
  CASE (    2)
    GO TO  320
  CASE (    3)
    GO TO  340
  CASE (    4)
    GO TO  360
END SELECT
300 DO  k = 1, nrhv
  zs(l+k)  = 1.0
  l        = l + nrhvwd
END DO
GO TO 400
320 DO  k = 1, nrhv
  zd(l+k)  = 1.0D+0
  l        = l + nrhvwd
END DO
GO TO 400
340 DO  k    = 1, nrhv
  zs(l+2*k-1) = 1.0
  l           = l + nrhvwd
END DO
GO TO 400
360 DO  k    = 1, nrhv
  zd(l+2*k-1) = 1.0D+0
  l           = l + nrhvwd
END DO
GO TO 400

!    CLOSE RIGHT HAND VECTORS FILE (DBB).
!    START FORWARD-BACKWARD SUBSTITUTION ON RIGHT HAND VECTORS

390 CALL CLOSE  (dbb,rew)
400 CONTINUE
j      = typel
nvec   = nrhv
nvecsz = ncol
IF ( switch > 1 ) nvec = nvec*2
SELECT CASE ( j )
  CASE (    1)
    GO TO  410
  CASE (    2)
    GO TO  420
  CASE (    3)
    GO TO  430
  CASE (    4)
    GO TO  440
END SELECT
410 CONTINUE
CALL fbsi1 ( BLOCK, zs, zs(mtria), zs(mtria), zs(buf2) )
GO TO 500
420 CONTINUE
CALL fbsi2 ( BLOCK, zs, zs(mtria), zs(mtria), zs(buf2) )
GO TO 500
430 CONTINUE
CALL fbsi3 ( BLOCK, zs, zs(mtria), zs(mtria), zs(buf2) )
GO TO 500
440 CONTINUE
CALL fbsi4 ( BLOCK, zs, zs(mtria), zs(mtria), zs(buf2) )
GO TO 500

!     OPEN AND PACK SOLUTION VECTORS ONTO OUTPUT FILE (DBX)

500 CALL gopen ( dbx, zs(buf1), wrtrew)
SELECT CASE ( switch )
  CASE (    1)
    GO TO  510
  CASE (    2)
    GO TO  530
  CASE (    3)
    GO TO  560
END SELECT

!     NORMAL CASE - CALL PACK

510 DO  l = 1, last, nrhvwd
  CALL pack ( zs(l), dbx, dbx )
END DO
GO TO 600

!     SPECIAL CASE - LOWER TRIANGULAR MATRIX IS RSP AND VECTORS ARE CSP

530 l = 0
DO  k = 1, nrhv
  CALL bldpk ( csp, typex, dbx, 0, 0 )
  DO  i = 1, ncol
    ys(1)    = zs(l+i     )
    ys(2)    = zs(l+i+ncol)
    iy       = i
    CALL zblpki
  END DO
  CALL bldpkn ( dbx, 0, dbx )
  l        = l + 2*ncol
END DO
GO TO 600

!     SPECIAL CASE - LOWER TRIANGULAR MATRIX IS RDP AND VECTORS ARE CDP

560 l = 0
DO  k = 1, nrhv
  CALL bldpk ( cdp, typex, dbx, 0, 0 )
  DO  i = 1,ncol
    yd(1)    = zd(l+i     )
    yd(2)    = zd(l+i+ncol)
    iy       = i
    CALL zblpki
  END DO
  CALL bldpkn ( dbx, 0, dbx )
  l        = l + 2*ncol
END DO
GO TO 600

!     JOB DONE. CLOSE TRIANGULAR MATRIX AND SOLUTION FILE.

600 CALL CLOSE ( dbx, rew )
subnam( 2 ) = END
CALL conmsg ( subnam, 2, 0 )
RETURN
END SUBROUTINE fbsi
