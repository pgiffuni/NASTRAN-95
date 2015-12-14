SUBROUTINE bufchk
     
 IMPLICIT INTEGER (a-z)
 LOGICAL :: od,os,OR,strdat
 INTEGER :: i(4)
 REAL :: s(4)
 DOUBLE PRECISION :: d(2)
 CHARACTER (LEN=8) :: nam,iibl
 COMMON /zzzzzz/  b(1)
 COMMON /bufcom/  offset,outbgn,outend,outlvl
 EQUIVALENCE      (i(1),s(1),d(1))
 DATA             nout    ,  nam      , iibl      /  &
     6       ,  'BUFCHK@','II,B(L)=' /
! VAX:
 DATA             rectrl  ,  rctrll   , rctrlc  , colhdr, coltrl /  &
     '1'x    ,  '2'x     , '3'x    , '4'x  , '8'x   /
 DATA             rechdr  ,  rchdst   , strdum  , eobstr   /  &
     'F1111'x,  'F2222'x , 'FAAAA'x, 'FBBBB'x /
 DATA             eob     ,  eof      , strhdr  , strtrl   /  &
     'F5555'x,  'F7777'x , 'F8888'x, 'F9999'x /
! UNIX:
!     DATA             RECTRL  ,  RCTRLL   , RCTRLC  , COLHDR, COLTRL /
!    1                 X'1'    ,  X'2'     , X'3'    , X'4'  , X'8'   /
!     DATA             RECHDR  ,  RCHDST   , STRDUM  , EOBSTR   /
!    1                 X'F1111',  X'F2222' , X'FAAAA', X'FBBBB' /
!     DATA             EOB     ,  EOF      , STRHDR  , STRTRL   /
!    1                 X'F5555',  X'F7777' , X'F8888', X'F9999' /
 
!*****
 lshift(k,j) = ishft(k, j)
 rshift(k,j) = ishft(k,-j)
!     WHERE         ISHFT(K,+J) IS  LEFT-SHIFT K BY J BITS, ZERO FILL
!                   ISHFT(K,-J) IS RIGHT-SHIFT K BY J BITS, ZERO FILL
!     AND           ISHFT IS SYSTEM ROUTINE
 
! UNIX:
!     REMOVE ABOVE 2 ON-LINE FUNCTIONS IF THE SYSTEM ISHFT FUNCTION IS
!     NOT AVAILABLE.  LSHIFT AND RSHIFT ARE ALREADY ENTRY POINTS IN
!     SUBROUTINE MAPFNS.
!*****
 
 l      = offset - 1
 datbgn = l + 1
 datend = l + 8
 dattyp = 5
 skp    = 0
 
 DO  ii = 1,10000
   IF (ii <= outend) GO TO 110
   IF (od) WRITE (nout,100) nam,ii
   100 FORMAT (5X,a7,'100  LIMIT POINTER REACHED.  II=',i6)
   GO TO 900
   110 l   = l + 1
   od  = ii >= outbgn .AND. outlvl > 0
   os  = ii >= outbgn .AND. outlvl > 1
   OR  = ii >= outbgn .AND. outlvl > 2
   skp = skp - 1
   IF (skp > 0) CYCLE
   IF (l >= datbgn .AND. l <= datend) GO TO 500
   w   = b(l)
   f1  = rshift(       w    ,28)
   f2  = rshift(lshift(w, 4),16)
   f3  = rshift(lshift(w,20),20)
   f12 = rshift(       w    ,12)
   f31 = rshift(lshift(w,20),24)
   f32 = rshift(lshift(w,28),28)
   IF (f12 /= rechdr) GO TO 160
   datbgn = l + 1
   datend = l + f3
   dattyp = 5
   IF (od) WRITE (nout,150) nam,iibl,ii,b(l)
   150 FORMAT (5X,a7,'150  REC HDR - NO STRING.  ',a8,i6,z8)
   CYCLE
   160 IF (f12 /= rchdst) GO TO 180
   IF (od) WRITE (nout,170) nam,iibl,ii,b(l)
   170 FORMAT (5X,a7,'170  REC HDR - STRING.  ',a8,i6,z8)
   strdat = .true.
   CYCLE
   180 IF (.NOT.(f1 == rectrl .OR. f1 == rctrll)) GO TO 210
   IF (od) WRITE (nout,200) nam,iibl,ii,b(l)
 200 FORMAT (5X,a7,'200  REC TRAILR - END OF RECORD.  ',a8,i6,z8)
 CYCLE
 210 IF (f1 /= rctrlc) GO TO 240
 IF (od) WRITE (nout,230) nam,iibl,ii,b(l)
 230 FORMAT (5X,a7,'230  REC TRAILER - RECORD CONTINUES.  ',a8,i6,z8)
 CYCLE
 240 IF (f12 /= eob) GO TO 260
 IF (od) WRITE (nout,250) nam,iibl,ii,b(l)
250 FORMAT (5X,a7,'250  END OF BLOCK.  ',a8,i6,z8)
CYCLE
260 IF (f12 /= eof) GO TO 280
IF (od) WRITE (nout,270) nam,iibl,ii,b(l)
270 FORMAT (5X,a7,'270  END OF FILE.  ',a8,i6,z8)
CYCLE
280 IF (f1 /= colhdr) GO TO 310
IF (od) WRITE (nout,300) nam,iibl,ii,b(l),f2
300 FORMAT (5X,a7,'300  COLUMN HEADER.  ',a7,',F2=',i6,z8,i8)
dattyp = f31
CYCLE
310 IF (f1 /= coltrl) GO TO 340
IF (od) WRITE (nout,330) nam,iibl,ii,b(l),f2
330 FORMAT (5X,a7,'330  COLUMN TRAILER.  ',a7,',F2=',i6,z8,i8)
CYCLE
340 IF (f12 /= strhdr) GO TO 360
IF (od) WRITE (nout,350) nam,iibl,ii,b(l),b(l+1),f3
350 FORMAT (5X,a7,'350  STRING HEADER.  ',a7,',B(L+1),F3=',i6,2Z8,i8)
skp    = 2
datbgn = l + 1
datend = l + f3
CYCLE
360 IF (f12 /= strtrl) GO TO 380
IF (od) WRITE (nout,370) nam,iibl,ii,b(l),b(l+1),f3
370 FORMAT (5X,a7,'370  STRING TRAILER.  ',a7,',B(L+1),F3=',i6,2Z8,i8)
skp = 2
CYCLE
380 IF (f12 /= strdum) GO TO 410
IF (od) WRITE (nout,400) nam,iibl,ii,b(l)
400 FORMAT (5X,a7,'400  STRING DUMMY WORD.  ',a8,i6,z8)
CYCLE
410 WRITE (nout,430) nam,iibl,ii,b(l)
430 FORMAT (5X,a7,'430  INVALID CONTROL WORD.  ',a8,i6,z8)
GO TO 800

500 CONTINUE
DO  j = 1,4
  i(j) = b(l+j-1)
END DO
SELECT CASE ( dattyp )
  CASE (    1)
    GO TO 610
  CASE (    2)
    GO TO 620
  CASE (    3)
    GO TO 630
  CASE (    4)
    GO TO 640
  CASE (    5)
    GO TO 650
END SELECT
WRITE (nout,550) nam,dattyp
550 FORMAT (5X,a7,'550  BAD DATA TYP,  DATTYP=',i6)
GO TO 800
610 IF (os) WRITE (nout,615) nam,ii,s(1)
615 FORMAT (5X,a7,'615  II,S(1)=',i6,e13.6)
CYCLE
620 IF (os) WRITE (nout,625) nam,ii,d(1)
625 FORMAT (5X,a7,'625  II,D(1)=',i6,d17.9)
skp = 2
CYCLE
630 IF (os) WRITE (nout,635) nam,ii,s(1),s(2)
635 FORMAT (5X,a7,'635  II,S(1),S(2)=',i6,2E13.6)
skp = 2
CYCLE
640 IF (os) WRITE (nout,645) nam,ii,d(1),d(2)
645 FORMAT (5X,a7,'645  II,D(1),D(2)=',i6,2D17.9)
skp = 4
CYCLE
650 IF (OR) WRITE (nout,660) nam,ii,i(1)
660 FORMAT (5X,a7,'615  II,I(1)=',i6,i14)
END DO

WRITE (nout,750) nam
750 FORMAT (5X,a7,'750  BLOCK TOO LONG')
800 CALL vaxend
900 RETURN
END SUBROUTINE bufchk
