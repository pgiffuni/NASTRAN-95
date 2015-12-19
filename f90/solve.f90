SUBROUTINE solve
     
!     SOLVE IS A DMAP DRIVER TO SOLVE THE MATRIX EQUATION AX=B
 
!     SOLVE   A,B/X/SYM/SIGN/PREC/TYPE $
 
!     SYM     =  1 - USE SYMETRIC DECOMPOSITION
!                0 - CHOOSE WHICH DECOMPOSITION BASED ON INPUT MATRIX
!               -1 - USE UNSYMETRIC DECOMPOSITION
!     ISIGN   =  1   SOLVE AX = B
!               -1   SOLVE AX =-B
!     IPREC   =  PRECISION USED IN THE FBS PASS
!     ITYPE   =  DESIRED TYPE OF THE OUTPUT MATRIX X
 
 
 INTEGER :: NAME(2)   ,rect     ,a        ,b        ,  &
     cdp       ,rdp      ,sym      ,sqr      , dosi(3)   ,refus(3) ,outpt    ,x
 REAL :: zz(1)     ,zzz(1)   ,zzzz(1)  ,zzzzz(1)
 DOUBLE PRECISION :: det       ,dett     ,mindia   ,cdet     ,  &
     cmndia    ,detc     ,minds
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm 
 COMMON /xmssg /  ufm       ,uwm      ,uim      ,sfm      , swm
 COMMON /BLANK /  isym      ,ksign    ,iprec    ,itype
 COMMON /system/  ksystm(65)
 COMMON /sfact /  ifila(7)  ,ifill1(7),ifilc(7) ,iscr11   ,  &
     iscr22    ,nz       ,det      ,detc     ,  &
     ipower    ,iscr33   ,minds    ,ichol
 COMMON /fbsx  /  ifill(7)  ,ifillt(7),ifilb(7) ,ifilx(7) ,  &
     nx        ,iprec1   ,isign1   ,iscr
 COMMON /dcompx/  ia(7)     ,il(7)    ,iu(7)    ,isr1     ,  &
     isr2      ,isr3     ,dett     ,ipow     ,  &
     nzz       ,mindia   ,ib       ,ibbar
 COMMON /cdcmpx/  ja(7)     ,kl(7)    ,ku(7)    ,jscr1    ,  &
     jscr2     ,jscr3    ,cdet(2)  ,jpow     ,  &
     nzzzz     ,cmndia   ,jbb      ,jbbar
 COMMON /gfbsx /  jl(7)     ,ju(7)    ,jb(7)    ,jx(7)    ,  &
     nzzz      ,ipr      ,isgn
 COMMON /names /  rd        ,rdrew    ,wrt      ,wrtrew   ,  &
     rew       ,norew    ,eofnrw   ,rsp      ,  &
     rdp       ,csp      ,cdp      ,sqr      ,  &
     rect      ,diag     ,lower    ,upper    , sym       ,row      ,identy
 COMMON /zzzzzz/  z(1)
 EQUIVALENCE      (zz(1),z(1))
 EQUIVALENCE      (zzz(1),z(1))
 EQUIVALENCE      (zzzz(1),z(1))
 EQUIVALENCE      (zzzzz(1),z(1))
 EQUIVALENCE      (ksystm(55),kprec)  ,(ksystm(2),outpt)
 DATA    a,b,x /  101,102,201/,     NAME  / 4HSOLV,4HE    /
 DATA    iscr1 ,  iscr2,iscr3,iscr4,iscr5 / 301   ,  302  ,303  ,304  ,305   /
 DATA    dosi  /  4HSING , 4HDOUB , 4HMLTP/,  &
     refus /  2*3H   , 3HREF          /
 
 
 ja(1) = a
 CALL rdtrl (ja)
 
 iform = ja(4)
 IF (isym < 0) THEN
   GO TO     1
 ELSE IF (isym == 0) THEN
   GO TO     5
 ELSE
   GO TO     3
 END IF
 1 IF (iform == sym) WRITE (outpt,2) uwm,NAME
 2 FORMAT (a25,' 2340, MODULE ',2A4,' HAS BEEN REQUESTED TO DO ',  &
       'UNSYMETRIC DECOMPOSITION OF A SYMETRIC MATRIX.')
   iform = rect
   IF (ja(2) == ja(3)) iform = sqr
   GO TO 5
   3 IF (ja(2) == ja(3) .AND. iform /= sym) WRITE (outpt,4) swm,NAME
   4 FORMAT (a27,' 2341, MODULE ',2A4,' HAS BEEN FURNISHED A SQUARE ',  &
       'MATRIX MARKED UNSYMETRIC FOR SYMETRIC DECOMPOSITION.')
   iform = sym
   5 isym  = -1
   IF (iform == sym) isym = 1
   ja(4) = iform
   IF (isym < 0) GO TO 30
   
!     SET UP CALL TO SDCOMP AND FBS
   
   INDEX = 1
   ichol = 0
   DO    i = 1,7
     ifila(i) = ja(i)
   END DO
   n = ifila(2)
   ifill1(1) = iscr1
   ifilc(1)  = iscr2
   iscr11 = iscr3
   iscr22 = iscr4
   iscr33 = iscr5
   nz = korsz(z)
   ifill1(5) = ifila(5)
   CALL sdcomp (*20,z,z,z)
   ifill1(3) = ifill1(2)
   ifill1(4) = lower
   CALL wrttrl (ifill1)
   ifill(1) = iscr1
   CALL rdtrl (ifill)
   ifilb(1) = b
   CALL rdtrl (ifilb)
   
!     IF THE B MATRIX IS PURGED, ASSUME AN IDENTITY MATRIX IN ITS PLACE
   
   IF (ifilb(1) <= 0) CALL makmcb (ifilb,b,n,identy,ja(5))
   isign1 = ksign
   ia5 = ifila(5)
   ib5 = ifilb(5)
   
!     DETERMINE THE PRECISION FOR THE CALCULATIONS
!     AND THE TYPE OF THE OUTPUT MATRIX
   
   200 iprec1 = kprec
   IF ((ia5 > 0 .AND. ia5 <= 4) .OR. (ib5 > 0 .AND. ib5 <= 4)) iprec1 = 1
   IF (ia5 == 2 .OR. ia5 == 4 .OR. ib5 == 2 .OR. ib5 == 4) iprec1 = 2
   IF (iprec == iprec1 .OR. iprec == 0) GO TO 222
   IF (iprec < 1 .OR. iprec > 2) iprec = 3
   WRITE  (outpt,221) swm,dosi(iprec),refus(iprec),NAME,dosi(iprec1)
   221 FORMAT (a27,' 2163, REQUESTED ',a4,'LE PRECISION ',a3,'USED BY ',  &
       2A4,2H. ,a4,'LE PRECISION IS LOGICAL CHOICE')
   IF (iprec /= 3 ) iprec1 = iprec
   222 iprec = iprec1
   ltype = iprec1
   IF (ia5 == 3 .OR. ia5 == 4 .OR. ib5 == 3 .OR. ib5 == 4) ltype = iprec1 + 2
   IF (itype == 0 .OR. itype == ltype) GO TO 224
   jj = 1
   IF (itype < 1 .OR. itype > 4 ) jj = 3
   WRITE  (outpt,223) sfm,itype,refus(jj),NAME,ltype
   223 FORMAT (a27,' 2164, REQUESTED TYPE ',i4,2H, ,a3,'USED BY ',2A4,  &
       '. TYPE ',i4,' IS LOGICAL CHOICE.')
   IF (jj /= 3 ) ltype = itype
   224 itype = ltype
   IF (INDEX == 2) GO TO 45
   
!     DEFINE THE MATRIX CONTROL BLOCK FOR THE OUTPUT MATRIX
   
   CALL makmcb (ifilx,x,n,rect,itype)
   nx = korsz(zz)
   IF (ifilb(4) == identy) ifilb(5) = iprec
   iscr = iscr1
   CALL fbs (zz,zz)
   IF (ifilx(2) == n) ifilx(4) = sqr
   CALL wrttrl (ifilx)
   RETURN
   
   20 no = ISIGN(5,isym)
   isym = -1
   CALL mesage (no,a,NAME)
   
!     SET UP THE CALL TO DECOMP AND GFBS
   
   30 CONTINUE
   INDEX = 2
   IF (ja(5) > 2) GO TO 80
   ia(1) = a
   il(1) = iscr1
   iu(1) = iscr2
   isr1  = iscr3
   isr3  = iscr5
   isr2  = iscr4
   nzz   = korsz(zzz)
   CALL rdtrl (ia)
   ia(4) = sqr
   n     = ia(2)
   il(5) = ja(5)
   ib    = 0
   ibbar = 0
   CALL decomp (*20,zzz,zzz,zzz)
   DO  i = 1,7
     jl(i) = il(i)
     ju(i) = iu(i)
   END DO
   40 jb(1) = b
   CALL rdtrl (jb)
   
!     IF THE B MATRIX IS PURGED, ASSUME AN IDENTITY MATRIX IN ITS PLACE
   
   IF (jb(1) <= 0) CALL makmcb (jb,b,n,identy,ja(5))
   ia5  = ja(5)
   ib5  = jb(5)
   isgn = ksign
   
!     DETERMINE THE PRECISION FOR THE CALCULATIONS
!     AND THE TYPE OF THE OUTPUT MATRIX
   
   GO TO 200
   45 ipr = iprec
   
!     DEFINE THE MATRIX CONTROL BLOCK FOR THE OUTPUT MATRIX
   
   CALL makmcb (jx,x,n,rect,itype)
   nzzz = korsz(zzzz)
   IF (jb(4) == identy) jb(5) = iprec
   CALL gfbs (zzzz,zzzz)
   IF (jx(2) == n) jx(4) =  sqr
   CALL wrttrl (jx)
   RETURN
   
!     SET UP CALL TO CDCOMP AND GFBS
   
   80 CONTINUE
   kl(1) = iscr1
   ku(1) = iscr2
   jscr1 = iscr3
   jscr2 = iscr4
   jscr3 = iscr5
   nzzzz = korsz(zzzzz)
   ja(4) = sqr
   n     = ja(2)
   kl(5) = ja(5)
   jbb   = 0
   jbbar = 0
   CALL cdcomp (*20,zzzzz,zzzzz,zzzzz)
   DO  i = 1, 7
     jl(i) = kl(i)
     ju(i) = ku(i)
   END DO
   GO TO 40
 END SUBROUTINE solve
