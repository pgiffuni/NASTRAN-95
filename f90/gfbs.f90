SUBROUTINE gfbs (x,dx)
     
!     GIVEN THE TRIANGULAR FACTORS FOR A GENERAL MATRIX, GFBS WILL
!     PERFORM THE FORWARD-BACKWARD SUBSTITUTION NECESSARY TO SOLVE
!     A SYSTEM OF EQUATIONS
 
!     DEFINITION OF INPUT PARAMETERS
 
!     FILEL    =  MATRIX CONTROL BLOCK FOR THE LOWER TRIANGLE L
!     FILEU    =  MATRIX CONTROL BLOCK FOR THE UPPER TRIANGLE U
!     FILEB    =  MATRIX CONTROL BLOCK FOR THE LOAD   VECTORS B
!     FILEX    =  MATRIX CONTROL BLOCK FOR THE SOLUTION VECTORS X
!     NX       =  NUMBER OF CELLS OF CORE AVAILABLE AT X
!     PREC     =  DESIRED PRECISION OF ARITHMETIC OPERATIONS
!                (1 = SINGLE PRECISION, 2 = DOUBLE PRECISION)
!     ISIGN    =  SIGN TO BE APPLIED TO THE LOAD VECTORS
!     X        =  BLOCK OF CORE AVAILABLE AS WORKING STORAGE
!     DX       =  SAME BLOCK AS X, BUT TYPED DOUBLE PRECISION
 
 
 REAL, INTENT(OUT)                        :: x(1)
 DOUBLE PRECISION, INTENT(OUT)            :: dx(1)
 INTEGER :: filel     ,fileu    ,fileb    ,filex    ,  &
     typea     ,type1    ,type2    ,formb    ,  &
     sysbuf    ,prec     ,eol      ,typear   ,  &
     typex     ,typel    ,rc       ,rew      ,  &
     typeb     ,tra1     ,tra2     ,tra3     ,  &
     tra4      ,tra5     ,parm(4)  ,CMPLX    ,  &
     eofnrw    ,col      ,fstcol   ,clsop
 REAL :: zeros(4)  , subnam(2) ,buf(2)   ,begn     ,END
 DOUBLE PRECISION :: da(2)    ,dtemp
 
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON   /xmssg /  ufm       ,uwm      ,uim      ,sfm
 COMMON   /system/  sysbuf    ,nout
 COMMON   /TYPE  /  prc(2)    ,nwds(4)  ,rc(10)
!     COMMON   /DESCRP/  LENGTH    ,MAJOR
 COMMON   /names /  rd        ,rdrew    ,wrt      ,wrtrew   ,  &
     rew       ,norew    ,eofnrw   ,rsp      ,  &
     rdp       ,csp      ,cdp      ,sqr      ,  &
     rect      ,diag     ,lowtri   ,uprtri   , sym       ,row      ,identy
 COMMON   /unpakx/  typea     ,ixy      ,jxy      ,incry
 COMMON   /packx /  type1     ,type2    ,iy       ,jy       , incrx
 COMMON   /zntpkx/  a(4)      ,ii       ,eol
 COMMON   /gfbsx /  filel(7)  ,fileu(7) ,fileb(7) ,filex(7) ,  &
     nx        ,prec     ,ISIGN
 EQUIVALENCE        (a(1),da(1))        ,(filel(5),typel)   ,  &
     (filel(3),nrow)     ,(filex(5),typex)   ,  &
     (fileb(4),formb)    ,(fileb(5),typeb)
 DATA      parm(3), parm(4)  /4HGFBS,4H     /
 DATA      zeros /  0., 0., 0., 0. /
DATA      subnam/  4HGFBS,4H      /, begn/ 4HBEGN/, END/ 4HEND /

buf(1) = subnam(1)
buf(2) = begn
CALL conmsg (buf,2,0)

!     INITIALIZE

IF (formb == identy)typeb = 1
typear = prec
IF (rc(typel)+rc(typeb)-1 > 1) typear = prec + 2
incr  = nwds(typear)*nrow
typea = typear*ISIGN
type1 = typear
type2 = typex
incrx = 1
incry = 1
CMPLX = rc(typear)
iobuf = nx - sysbuf
icol  = iobuf - 1
col   = 1
clsop = eofnrw

!     SET UP TRANSFER VECTORS FOR THE ARITHMETIC TYPES

SELECT CASE ( typear )
  CASE (    1)
    GO TO 10
  CASE (    2)
    GO TO 20
  CASE (    3)
    GO TO 30
  CASE (    4)
    GO TO 40
END SELECT
10 ASSIGN 120 TO tra1
ASSIGN 240 TO tra2
ASSIGN 330 TO tra3
ASSIGN 430 TO tra4
ASSIGN 540 TO tra5
GO TO 50
20 ASSIGN 130 TO tra1
ASSIGN 250 TO tra2
ASSIGN 340 TO tra3
ASSIGN 440 TO tra4
ASSIGN 550 TO tra5
GO TO 50
30 ASSIGN 140 TO tra1
ASSIGN 260 TO tra2
ASSIGN 350 TO tra3
ASSIGN 450 TO tra4
ASSIGN 560 TO tra5
GO TO 50
40 ASSIGN 150 TO tra1
ASSIGN 270 TO tra2
ASSIGN 360 TO tra3
ASSIGN 460 TO tra4
ASSIGN 570 TO tra5
50 CONTINUE
nm = (iobuf-1)/incr
IF (nm <= 0) GO TO 640
noload = fileb(2)
IF (formb == identy)noload = nrow
ident  = 1
lstlod = noload

!     WRITE OUTPUT HEADER RECORDS AND INITIALIZE MATRIX CONTROL BLOCKS

CALL gopen (filex,x(iobuf),1)
CALL CLOSE (filex(1),norew)
filex(2) = 0
filex(6) = 0
filex(7) = 0
IF (formb == identy) GO TO 100

!     OPEN THE LOAD FILE AND FILL CORE WITH LOAD VECTORS

CALL gopen (fileb,x(iobuf),0)
60 nn  = 0
khr = icol
fstcol = col
l   = 1
ixy = 1
jxy = nrow
70 IF (l+incr >= khr) GO TO 85
CALL unpack (*80,fileb,x(l))
nn  = nn + 1
x(khr) = col
khr = khr - 1
l   = l + incr
80 IF (col == lstlod) GO TO 90
col = col + 1
GO TO 70
85 col  = col - 1
90 ncol = khr
x(ncol) = lstlod + 1
lstcol  = col
IF (lstcol == lstlod) clsop = rew
CALL CLOSE (fileb,clsop)
IF (nn == 0) GO TO 592
GO TO 180

!     GENERATE COLUMNS OF THE IDENTITY MATRIX

100 nn = MIN0(nm,noload)
l  = 1
DO  i = 1,nn
  j1 = l
  j2 = j1 + incr - 1
  DO  k = j1,j2
    x(k) = 0.
  END DO
  k = l + ident - 1
  GO TO tra1, (120,130,140,150)
  120 x(k) = 1.
  GO TO 160
  130 k = (l-1)/2 + ident
  dx(k) = 1.d0
  GO TO 160
  140 kk = k + ident - 1
  x(kk) = 1.
  GO TO 160
  150 kk = (l-1)/2 + 2*ident - 1
  dx(kk) = 1.d0
  160 ident  = ident + 1
  l = l + incr
END DO
fstcol = col
col    = ident - 1
lstcol = col
180 ijk    = 0

!     OPEN FILE FOR THE LOWER TRIANGLE

parm(2) = filel(1)
CALL gopen (filel,x(iobuf),0)

!     BEGIN FORWARD PASS

j = 1
190 CALL intpk (*380,filel(1),0,typear,0)
200 IF (eol == 0.0) THEN
  GO TO   210
ELSE
  GO TO   650
END IF
210 CALL zntpki
IF (j-ii < 0) THEN
  GO TO   310
ELSE IF (j-ii == 0) THEN
  GO TO   220
ELSE
  GO TO   200
END IF

!     PERFORM THE REQUIRED ROW INTERCHANGE

220 intchn = a(1)
k = 0
IF (prec == 2) intchn = da(1)
in1 = j*CMPLX
in2 = in1 + intchn*CMPLX
230 GO TO tra2, (240,250,260,270)
240 temp     = x(in1)
x(in1)   = x(in2)
x(in2)   = temp
GO TO 280
250 dtemp    = dx(in1)
dx(in1)  = dx(in2)
dx(in2)  = dtemp
GO TO 280
260 temp     = x(in1)
x(in1)   = x(in2)
x(in2)   = temp
temp     = x(in1-1)
x(in1-1) = x(in2-1)
x(in2-1) = temp
GO TO 280
270 dtemp    = dx(in1)
dx(in1)  = dx(in2)
dx(in2)  = dtemp
dtemp    = dx(in1-1)
dx(in1-1) = dx(in2-1)
dx(in2-1) = dtemp
280 in1 = in1 + nrow*CMPLX
in2 = in2 + nrow*CMPLX
k   = k + 1
IF (k-nn < 0) THEN
  GO TO   230
END IF
290 IF (eol == 0.0) THEN
  GO TO   300
ELSE
  GO TO   380
END IF
300 CALL zntpki
310 k = 0
in2 = j*CMPLX
in1 = ii*CMPLX
320 k = k + 1
GO TO tra3, (330,340,350,360)
330 x(in1) = x(in1) - x(in2)*a(1)
GO TO 370
340 dx(in1) = dx(in1) - dx(in2)*da(1)
GO TO 370
350 x(in1-1) = x(in1-1) - a(1)*x(in2-1) + a(2)*x(in2  )
x(in1  ) = x(in1  ) - a(1)*x(in2  ) - a(2)*x(in2-1)
GO TO 370
360 dx(in1-1) = dx(in1-1) - da(1)*dx(in2-1) + da(2)*dx(in2  )
dx(in1  ) = dx(in1  ) - da(1)*dx(in2  ) - da(2)*dx(in2-1)
370 in1 = in1 + nrow *CMPLX
in2 = in2 + nrow *CMPLX
IF (k-nn < 0) THEN
  GO TO   320
ELSE
  GO TO   290
END IF
380 j = j + 1
IF (j < nrow) GO TO 190
CALL CLOSE (filel(1),rew)

!     BEGIN BACKWARD PASS

ioff = fileu(7)-1
parm(2) = fileu(1)
CALL gopen (fileu,x(iobuf),0)
j = nrow
390 CALL intpk (*650,fileu(1),0,typear,0)
IF (eol == 0.0) THEN
  GO TO   410
ELSE
  GO TO   650
END IF
410 CALL zntpki
i = nrow - ii + 1
IF (i /= j) GO TO 510

!     DIVIDE BY THE DIAGONAL

in1 = i*CMPLX
k   = 0
420 GO TO tra4, (430,440,450,460)
430 x(in1) = x(in1)/a(1)
GO TO 470
440 dx(in1) = dx(in1)/da(1)
GO TO 470
450 temp   = (a(1)*x(in1-1) + a(2)*x(in1  ))/(a(1)*a(1) + a(2)*a(2))
x(in1) = (a(1)*x(in1  ) - a(2)*x(in1-1))/(a(1)*a(1) + a(2)*a(2))
x(in1-1) = temp
GO TO 470
460 dtemp   = (da(1)*dx(in1-1) + da(2)*dx(in1  ))/(da(1)**2 +da(2)**2)
dx(in1) = (da(1)*dx(in1  ) - da(2)*dx(in1-1))/(da(1)**2 +da(2)**2)
dx(in1-1) = dtemp
470 k = k + 1
in1 = in1 + nrow*CMPLX
IF (k-nn < 0) THEN
  GO TO   420
ELSE
  GO TO   490
END IF

!     SUBTRACT OFF REMAINING TERMS

480 IF (i > j) GO TO 410
490 IF (eol == 0.0) THEN
  GO TO   500
ELSE
  GO TO   590
END IF
500 CALL zntpki
i   = nrow - ii + 1
510 in1 = i*CMPLX
in2 = j*CMPLX
IF (i < j) GO TO 520
k   = in1
in1 = in2 - ioff*CMPLX
in2 = k
520 k   = 0
530 GO TO tra5, (540,550,560,570)
540 x(in1) = x(in1) - a(1)*x(in2)
GO TO 580
550 dx(in1) = dx(in1) - dx(in2)*da(1)
GO TO 580
560 x(in1-1) = x(in1-1) - a(1)*x(in2-1) + a(2)*x(in2  )
x(in1  ) = x(in1  ) - a(1)*x(in2  ) - a(2)*x(in2-1)
GO TO 580
570 dx(in1-1) = dx(in1-1) - da(1)*dx(in2-1) + da(2)*dx(in2  )
dx(in1  ) = dx(in1  ) - da(1)*dx(in2  ) - da(2)*dx(in2-1)
580 in1 = in1 + nrow*CMPLX
in2 = in2 + nrow*CMPLX
k   = k + 1
IF (k-nn < 0) THEN
  GO TO   530
ELSE
  GO TO   480
END IF
590 j = j - 1
IF (j > 0) GO TO 390
CALL CLOSE (fileu(1),rew)

!     OUTPUT LOAD VECTORS

592 CONTINUE
CALL gopen (filex,x(iobuf),wrt)
l  = 1
iy = 1
IF (formb /= identy) nxtnz = x(icol)
khr = icol
DO  col = fstcol,lstcol
  IF (formb == identy) GO TO 595
! 593 CONTINUE
  IF (col-nxtnz < 0.0) THEN
    GO TO   594
  ELSE IF (col-nxtnz == 0.0) THEN
    GO TO   595
  ELSE
    GO TO   901
  END IF
  594 jy = 1
  CALL pack (zeros,filex,filex)
  CYCLE
  595 jy = nrow
  CALL pack (x(l),filex,filex)
  l   = l + incr
  khr = khr - 1
  IF (formb /= identy) nxtnz = x(khr)
END DO
IF (formb /= identy .AND. khr /= ncol) GO TO 902
IF (lstcol == lstlod) clsop = rew
CALL CLOSE (filex,clsop)
noload = noload - (lstcol-fstcol+1)
IF (lstcol == lstlod) GO TO 670
col = lstcol + 1
IF (formb == identy) GO TO 100
CALL gopen (fileb,x(iobuf),rd)
GO TO 60
640 parm(1) = -8
GO TO 660
650 parm(1) = -5
660 CALL mesage (parm(1),parm(2),parm(3))
670 IF (filex(2) /= lstlod) GO TO 903
buf(1) = subnam(1)
buf(2) = END
CALL conmsg (buf,2,0)
RETURN

!     LOGIC ERRORS LAND HERE

901 kerr = 593
GO TO 997
902 kerr = 600
GO TO 997
903 kerr = 670
GO TO 997
997 WRITE  (nout,998) sfm,kerr
998 FORMAT (a25,i4,' - LOGIC ERROR IN GFBS')
CALL mesage (-61,0,0)
RETURN
END SUBROUTINE gfbs
