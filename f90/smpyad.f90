SUBROUTINE smpyad
     
 INTEGER :: sreslt   ,sadd     ,tmat     ,trlra    ,trlrb  ,  &
     trlrc    ,trlrd    ,trnsp    ,signab   ,signc  ,  &
     scrtch   ,trlr(7,5),mat(5)   ,addmat   ,resmat ,  &
     recmat   ,dosi(3)  ,refus(3) ,outpt    ,NAME(2)
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm      ,uwm      ,uim      ,sfm      ,swm
 COMMON /mpyadx/ trlra(7) ,trlrb(7) ,trlrc(7) ,trlrd(7) ,na     ,  &
     trnsp    ,signab   ,signc    ,iprec1   ,scrtch /zzzzzz/ a(1)  &
     /BLANK / n        ,sreslt   ,sadd     ,iprec    ,tmat(4) /system/ ksystm(65)
 EQUIVALENCE     (ksystm(55),kprec) ,(ksystm(2),outpt)
 DATA    mat   / 101      ,102      ,103      ,104      ,105 /
 DATA    addmat, resmat   ,intres   ,mpyads   ,recmat        /  &
     106   , 201      ,301      ,302      ,2             /
 DATA    dosi  / 4HSING   ,4HDOUB   ,4HMLTP   /
 DATA    refus / 2*3H     ,3HREF    /
 DATA    NAME  / 4HSMPY   ,4HAD     /
 
 IF (n <= 1) GO TO 200
 IF (n > 5) n = 5
 iprec1 = 1
 itype  = 0
 
!     IF ONE OF THE -N- MATRICES IN THE PRODUCT DOES NOT EXIST,
!     SKIP THE ENTIRE CALCULATION.
 
 DO  i = 1,n
   trlr(1,i) = mat(i)
   CALL rdtrl (trlr(1,i))
   IF (trlr(1,i) <= 0 .OR. trlr(2,i) <= 0 .OR. trlr(3,i) <= 0) GO TO 200
   IF (trlr(5,i) == 2 .OR. trlr(5,i) == 4) iprec1 = 2
   IF (trlr(5,i) == 3 .OR. trlr(5,i) == 4) itype  = 2
 END DO
 
!     CHECK TO SEE IF THE INPUT MATRICES ARE CONFORMABLE
 
 nm1  = n - 1
 nogo = 0
 DO  i = 1,nm1
   icol = trlr(2,i)
   IF (tmat(i) /= 0) icol = trlr(3,i)
   irow = trlr(3,i+1)
   IF (i == nm1) GO TO 160
   IF (tmat(i+1) /= 0) irow = trlr(2,i+1)
   160 IF (icol /= irow) nogo = 1
 END DO
 trlrc(1) = addmat
 CALL rdtrl (trlrc)
 IF (trlrc(1) <= 0) GO TO 180
 irow = trlr(3,1)
 IF (tmat(1) /= 0) irow = trlr(2,1)
 icol = trlr(2,n)
 IF (irow /= trlrc(3) .OR. icol /= trlrc(2)) nogo = 1
 180 IF (nogo == 1) CALL mesage (-55,0,NAME)
 
 IF (iprec1 < 1 .OR. iprec1 > 2) iprec1 = kprec
 IF (iprec == iprec1 .OR. iprec == 0) GO TO 222
 IF (iprec < 1 .OR. iprec > 2) iprec = 3
 WRITE  (outpt,221) swm,dosi(iprec),refus(iprec),NAME,dosi(iprec1)
 221 FORMAT (a27,' 2163, REQUESTED ',a4,'LE PRECISION ',a3,'USED BY ',  &
     2A4,2H. ,a4,'LE PRECISION IS LOGICAL CHOICE')
 IF (iprec /= 3) iprec1 = iprec
 222 iprec = iprec1
 itype = itype + iprec1
 
!     SETUP THE MPYADX COMMON BLOCK.
 
 IF ((n+1)/2 == n/2) GO TO 105
 trlrb(1) = intres
 m = resmat
 GO TO 106
 105 trlrb(1) = resmat
 m = intres
 106 trlrc(1) = 0
 DO  i = 1,7
   trlrd(i) = trlr(i,n)
 END DO
 trlrd(4) = recmat
 na = korsz(a)
 signab = 1
 signc  = sadd
 scrtch = mpyads
 
!     DO THE N-1 MULTIPLICATIONS.
 
 DO  k = 2,n
   j = n - k + 1
   trlra(1) = trlr(1,j)
   IF (k /= 3) l = trlrb(1)
   IF (k == 3) l = m
   trlrb(1) = trlrd(1)
   trlrd(1) = l
   DO  i = 2,7
     trlra(i) = trlr(i,j)
     trlrb(i) = trlrd(i)
   END DO
   IF (k /= n) GO TO 111
   trlrc(1) = addmat
   CALL rdtrl (trlrc)
   IF (trlrc(1) < 0) trlrc(1) = 0
   trlrd(5) = itype
   signab   = sreslt
   GO TO 115
   111 trlrd(5) = iprec1
   IF (trlra(5) > 2 .OR. trlrb(5) > 2) trlrd(5) = iprec1 + 2
   115 trnsp = tmat(j)
   trlrd(3) = trlra(3)
   IF (trnsp /= 0) trlrd(3) = trlra(2)
   trlrd(2) = trlrb(2)
   CALL mpyad (a,a,a)
 END DO
 IF (trlrd(2) == trlrd(3)) trlrd(4) = 1
 CALL wrttrl (trlrd)
 200 RETURN
END SUBROUTINE smpyad
