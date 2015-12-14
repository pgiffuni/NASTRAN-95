SUBROUTINE dmpyad
     
!     DMPYAD IS THE DMAP DRIVER FOR MATRIX MULTIPLICATION.
 
!     COMMENTS FROM G.CHAN/UNISYS ABOUT PREC1 IN /MPYADX/     1/91
!     ACCORDING TO THE USER'S MANUAL ON P. 3.5-18
!       PREC1 = 0, PERFORM ARITHMETIC IN D.P. IF A,B OR C IS IN D.P.
!             = 1, PERFORM ARITHMETIC IN S.P.
!             = 2, PERFORM ARITHMETIC IN D.P.
!     HOWEVER, THE CODE BELOW ALWAYS SETS
!       PREC1 TO 2, IF ANY OF THE A,B OR C IS IN D.P. AND 1 OTHERWISE
!       IN SUBROUTINE MPYAD, PREC1 IS ALWAYS SET TO 1 FOR CDC MACHINE
 
!     IF ITYPE IN /BLANK/ IS 1 OR 3, MPYAD PRODUCT WILL BE OUPUT IN S.P.
!     AND IN D.P. OF IT IS 2 OR 4
!     IF ITYPE IS 0, MPYAD PRODUCT WILL BE IN S.P. ONLY IF ALL A, B, AND
!     C MATRICES ARE IN S.P. OTHERWISE, THE PRODUCT WILL BE IN D.P.
 
 IMPLICIT INTEGER (a-z)
 INTEGER :: NAME(2)      ,dosi(3)      ,refus(3)      ,  &
     p(4)  ,q(4)  ,r(4)  ,zz(1) ,zzz(1)
 REAL :: alp(1),bet(1)
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm   ,uwm   ,uim   ,sfm   ,swm
 COMMON /system/ ksystm(65)
 COMMON /BLANK / t     ,signab,signc ,itype
 COMMON /mpyadx/ a(7)  ,b(7)  ,c(7)  ,d(7)  ,nz    ,trnsp  ,  &
     sab   ,sc    ,prec1 ,scr
 COMMON /dmpyx / e(7)  ,f(7)  ,g(7)  ,nzz   ,flag  ,sgn
 COMMON /saddx / nomat ,nzzz  ,mcbs(67)
 COMMON /zzzzzz/ z(1)
 EQUIVALENCE     (zz(1),z(1))
 EQUIVALENCE     (zzz(1),z(1))
 EQUIVALENCE     (ksystm(55),kprec), (ksystm(2),outpt)
 EQUIVALENCE     (mcbs( 1),p(1)   ), (mcbs( 8),typa  ),  &
     (mcbs( 9),alp(1) ), (mcbs(13),q(1)  ),  &
     (mcbs(20),typb   ), (mcbs(21),bet(1)), (mcbs(61),r(1)   )
 DATA    filea , fileb,filec,filed,scrtch  /  &
     101   , 102  ,103  ,201  ,301     /
 DATA    NAME  / 4HMPYA,4HD   /
 DATA    dosi  / 4HSING,4HDOUB,4HMLTP /, refus/ 2*3H   ,3HREF /
 DATA    square, rect,diag,symm,ident /  1,2,3,6,8  /
 
 
!     READ TRAILERS FOR A, B AND C MATRICES.
 
 nz   = korsz(z)
 a(1) = filea
 CALL rdtrl (a)
 IF (a(1) /= filea) GO TO 230
 b(1) = fileb
 CALL rdtrl (b)
 IF (b(1) /= fileb) GO TO 230
 c(1) = filec
 c(5) = 0
 CALL rdtrl (c)
 IF (c(1) < 0) c(1) = 0
 d(1) = filed
 d(3) = a(3)
 IF (t /= 0) d(3) = a(2)
 d(4) = rect
 
!     CHECK FOR CONFORMABLE MATRICIES
 
 IF (((c(2) /= b(2) .OR.  c(3) /= d(3)) .AND. c(1) /= 0) .OR.  &
     (b(3) /= a(2) .AND. t == 0) .OR. (b(3) /= a(3) .AND. t /= 0))  &
     CALL mesage (-55,0,NAME)
 trnsp = t
 sab   = signab
 sc    = signc
 prec  = 1
 IF (itype == 0) prec = 0
 IF (itype == 2 .OR. itype == 4) prec = 2
 prec1 = MAX0(a(5),b(5),c(5))
 IF (prec1 > 2) prec1 = prec1 - 2
 IF (prec1 < 1 .OR. prec1 > 2) prec1 = kprec
 IF (prec == prec1 .OR. prec == 0) GO TO 20
 IF (prec < 1 .OR. prec > 2) prec = 3
 WRITE  (outpt,10) swm,dosi(prec),refus(prec),NAME,dosi(prec1)
 10 FORMAT (a27,' 2430, REQUESTED ',a4,'LE PRECISION ',a3,'USED BY ',  &
     2A4,2H. ,a4,'LE PRECISION IS LOGICAL CHOICE')
 IF (prec /= 3) prec1 = prec
 20 ltype = prec1
 IF (a(5) == 3 .OR. a(5) == 4 .OR. b(5) == 3 .OR. b(5) == 4 .OR.  &
     c(5) == 3 .OR. c(5) == 4) ltype = prec1 + 2
 IF (itype == 0 .OR. itype == ltype) GO TO 40
 jj = 1
 IF (itype < 1 .OR. itype > 4) jj = 3
 WRITE  (outpt,30) swm,itype,refus(jj),NAME,ltype
 30 FORMAT (a27,' 2431, REQUESTED TYPE ',i4,2H, ,a3,'USED BY ',2A4,  &
     7H. TYPE ,i4,'  IS LOGICAL CHOICE.')
 IF (jj /= 3) ltype = itype
 40 itype = ltype
 d(5)  = itype
 scr   = scrtch
 
!     IF NEITHER A NOR B IS DIAGONAL, CALL MPYAD AND RETURN.
 
 IF (a(4) == diag .OR. b(4) == diag) GO TO 100
 CALL mpyad (z,z,z)
 IF (d(2) /= d(3)) GO TO 60
 d(4) = square
 IF (c(4) == 0) c(4) = diag
 k = 0
 DO  i = 4,21,7
   j = a(i)
   IF (j /= symm .AND. j /= diag .AND. j /= ident) GO TO 60
   IF (j == symm) k = 1
 END DO
 IF (k == 1) d(4) = symm
 60 CALL wrttrl (d)
 RETURN
 
!     OTHERWISE, CALL DMPY FOR DIAGONAL MULTIPLICATION.
 
 100 DO  i = 1,7
   e(i) = a(i)
   f(i) = b(i)
   IF (a(4) == diag) GO TO 110
   e(i) = b(i)
   f(i) = a(i)
   110 g(i) = d(i)
 END DO
 nzz  = korsz(zz)
 sgn  = signab
 flag = 0
 IF (b(4) == diag) flag = 1
 IF (c(1) /=    0) g(1) = scrtch
 CALL dmpy (zz,zz)
 IF (g(2) /= g(3)) GO TO 130
 g(4) = square
 k = 0
 DO  i = 4,14,7
   j = e(i)
   IF (j /= symm .AND. j /= diag .AND. j /= ident) GO TO 130
   IF (j == symm) k = 1
 END DO
 IF (k == 1) g(4) = symm
 130 CALL wrttrl (g)
 
!     IF ADDITION REQUIRED, CALL ADD ROUTINE.
 
 IF (c(1) == 0) RETURN
 DO  i = 1,7
   p(i)  = g(i)
   q(i)  = c(i)
   r(i)  = d(i)
 END DO
 DO  i = 2,4
   alp(i)= 0.0
   bet(i)= 0.0
 END DO
 typa  = 1
 alp(1)= 1.0
 typb  = 1
 bet(1)= 1.0
 IF (signc < 0) bet(1) =-1.0
 nzzz  = korsz(zzz)
 nomat = 2
 CALL sadd (zzz,zzz)
 IF (r(2) /= r(3)) GO TO 230
 r(4) = square
 IF (p(4) == symm .AND. (q(4) == symm .OR. q(4) == diag .OR.  &
     q(4) == ident)) r(4) = symm
 CALL wrttrl (r)
 230 RETURN
END SUBROUTINE dmpyad
