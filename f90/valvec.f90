SUBROUTINE valvec
     
!     LARGE ORDER REAL SYMMETRIC EIGENVALUE-EIGENVECTOR PROBLEM
 
 INTEGER :: qr,mcb(8),tri(3),qrx(3),wil(3),val(3)
 DIMENSION       vcom(30)
 COMMON /zzzzzz/ a(1)
 COMMON /givn  / title(150)
!WKBR 2/94 SPR93027 COMMON /SYSTEM/ ISYS
 COMMON /system/ isys, idumm(53), iprec
 EQUIVALENCE     (md   ,title(3)), (vcom(1),title(101)),  &
     (n    ,vcom( 1)), (nv     ,vcom(   7)),  &
     (oeigs,vcom(11)), (nver   ,vcom(  13)),  &
     (never,vcom(14)), (iterm  ,vcom(  16))
 DATA    tri   / 4HTRID, 4HI   , 4H    /
 DATA    qrx   / 4HQRIT, 4HER  , 4H    /
 DATA    wil   / 4HWILV, 4HEC  , 4H    /
 DATA    val   / 4HVALV, 4HEC  , 4H    /
DATA    ibegn , iend  / 4HBEGN, 4HEND /

!     DEFINITION OF VARIABLES AND DATA FORMATS

!     MD       INPUT MATRIX
!     N        SIZE OF MATRIX
!     NV       NUMBER OF EIGENVECTORS DESIRED
!     OEIGS    EIGENVALUE SUMMARY FILE
!     A        OPEN CORE
!     ID       POINTER TO DIAGONALS      -- N OF THEM (D.P.)
!     IO       POINTER TO OFF-DIAGONALS  -- N OF THEM (D.P.)
!     IV       POINTER TO EIGENVALUES    -- N OF THEM (D.P.)
!     IL       POINTER TO ORDER FOUND ARRAY N OF THEM (S.P.)
!     I1 - I6  POINTS TO SCRATCH ARRAYS  -- 2XN LONG
!     NVER     NUMBER OF VECTORS    ERRORS
!     NEVER    NUMBER OF EIGENVALUE ERRORS
!     ITERM    REASON FOR TERMINATION

!     INITIALIZATION FOR VALVEC IN BLOCKDATA ROUTINE READBD

!     DATA
!    1 MO, MD,MR1, M1, M2, M3, M4,LGAMA,OEIGS,PHIA,ORDER,RSTRT,NCOL,MAX/
!    *301,304,202,303,307,308,309,  201,  204, 305,   -2,   0 ,   0,253/


val(3) = ibegn
CALL conmsg (val,3,0)
iterm  = 1
mcb(1) = md
CALL rdtrl (mcb(1))
n  = mcb(2)
n2 = n*iprec
id = 1
io = id + n2
iv = io + n2
il = iv + n2
i1 = il + n
IF ((i1 + 1)/2  == i1/2) i1 = i1 + 1
i2 = i1 + n2
i3 = i2 + n2
i4 = i3 + n2
i5 = i4 + n2
i6 = i5 + n2

!     TRIDIAGONALIZATION.

IF (n > 2) GO TO 101
!WKBD 2/94 SPR93027 CALL SMLEIG (A(ID),A(IO),A(IV))
!WKBNB 2/94 SPR93027
IF ( iprec == 2 ) CALL smleig (a(id),a(io),a(iv))
IF ( iprec == 1 ) CALL smleig1(a(id),a(io),a(iv))
!WKBNE 2/94 SPR93027

IF (n-2 == 0) THEN
  GO TO   200
ELSE
  GO TO   300
END IF
101 tri(3) = ibegn
CALL conmsg (tri,3,0)
!WKBD 2/94 SPR93027 CALL TRIDI (A(ID),A(IO),A(IV),A(IL),A(I1),A(IL))
!WKBNB 2/94 SPR93027
IF ( iprec == 2 ) CALL tridi (a(id),a(io),a(iv),a(il),a(i1),a(il))
!                   D      O    V     A     B
IF ( iprec == 1 ) CALL tridi1(a(id),a(io),a(iv),a(il),a(i1),a(il))
!                   D      O    V     A     B
!WKBNE 2/94 SPR93027
tri(3) = iend
CALL conmsg (tri,3,0)

!     EIGENVALUES

200 qr =  0
IF (n <= 2) qr = 1
qrx(3) = ibegn
CALL conmsg (qrx,3,0)
!WKBD 2/94 SPR93027 CALL QRITER (A(IV),A(I1),A(IL),QR)
!WKBNB 2/94 SPR93027
IF ( iprec == 2 ) CALL qriter (a(iv),a(i1),a(il),qr)
IF ( iprec == 1 ) CALL qriter1(a(iv),a(i1),a(il),qr)
!WKBNE 2/94 SPR93027

qrx(3) = iend
CALL conmsg (qrx,3,0)
rstrt  = 0
wil(3) = ibegn
CALL conmsg (wil,3,0)

!     EIGENVECTORS

!WKBDB 2/94 SPR93027
!     CALL WILVEC (A(ID),A(IO),A(IV),A(IL),A(I1),A(I2),A(I3),A(I4),
!     1             A(I5),A(I6),N,A(I6))
!WKBDE 2/94 SPR93027
!WKBNB 2/94 SPR93027
IF ( iprec == 1 ) THEN
!                    D      0    C    A      B  &
CALL wilvec1(a(id),a(io),a(iv),a(il),a(i1),a(i2),a(i3),a(i4),  &
    a(i5),a(i6),n,a(i6))
ELSEIF ( iprec == 2 ) THEN
!                    D      0    C    A      B  &
CALL wilvec (a(id),a(io),a(iv),a(il),a(i1),a(i2),a(i3),a(i4),  &
    a(i5),a(i6),n,a(i6))
ENDIF    
!WKBNE 2/94 SPR93027
wil(3) = iend
CALL conmsg (wil,3,0)
300 CONTINUE
CALL gopen (oeigs,a(1),1)
mcb(1) = 4
mcb(2) = n
mcb(3) = nv
mcb(4) = never
mcb(5) = nver
mcb(8) = iterm
CALL WRITE (oeigs,mcb,8,1)
CALL CLOSE (oeigs,1)
val(3) = iend
CALL conmsg (val,3,0)
RETURN
END SUBROUTINE valvec
