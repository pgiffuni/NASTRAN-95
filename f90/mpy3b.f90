SUBROUTINE mpy3b (z,iz,dz)
!*****
!    PROCESSES A AND PERFORMS FIRST PART OF PRODUCT.
!*****
 
 REAL, INTENT(OUT)                        :: z(1)
 INTEGER, INTENT(OUT)                     :: iz(1)
 DOUBLE PRECISION, INTENT(OUT)            :: dz(1)
 DOUBLE PRECISION :: da
 
 
 
 INTEGER :: filea,scr1,scr2,FILE
 INTEGER :: prec,precn
 INTEGER :: zpntrs
 INTEGER :: utyp,urow1,urown,uincr
 INTEGER :: eol,eor
 INTEGER :: preck
 
 
 
 LOGICAL :: first1,first2
 
 
 
 
 DIMENSION NAME(2)
 
 
! FILES
 COMMON / mpy3tl / filea(7),fileb(7),filee(7),filec(7),scr1,scr2,  &
     scr,lkore,code,prec,lcore,scr3(7),buf1,buf2, buf3,buf4,e
! SUBROUTINE CALL PARAMETERS
 COMMON / mpy3cp / itrl,icore,n,ncb,m,nk,dumcp(2),zpntrs(22),laend,  &
     first1,first2,k,k2,kcount,iflag,ka,kb,j,i,ntbu
! UNPACK
 COMMON / unpakx / utyp,urow1,urown,uincr
! TERMWISE MATRIX READ
 COMMON / zntpkx / a(2),dum(2),irow,eol,eor
 
 
 
 EQUIVALENCE     (a(1),da)
! OPEN CORE POINTERS
 EQUIVALENCE     (ibcols,zpntrs(11)),     (ibcid,zpntrs(13)),  &
     (ibntu,zpntrs(15)),      (iktbp,zpntrs(17)), (iakj,zpntrs(21))
 
 
 
 DATA NAME / 4HMPY3,4HB    /
!*****
!    INITIALIZATION.
!*****
 FILE = scr1
 utyp = prec
 urow1 = 1
 urown = n
 uincr = 1
 precn = prec*n
!*****
!    READ AND STORE COLUMN OF A.
!*****
 k = 0
 kt = iktbp - 1
 IF (prec == 2) GO TO 20
! SINGLE PRECISION CASE
 kj = iakj - 1
 CALL intpk(*120,filea,0,1,0)
 10 CALL zntpki
 k = k + 1
 kt = kt + 1
 iz(kt) = irow
 kj = kj + 1
 z(kj) = a(1)
 IF (eol == 1) GO TO 30
 GO TO 10
! DOUBLE PRECISION CASE
 20 kj = (iakj - 1)/2
 CALL intpk(*120,filea,0,2,0)
 25 CALL zntpki
 k = k + 1
 kt = kt + 1
 iz(kt) = irow
 kj = kj + 1
 dz(kj) = da
 IF (eol == 1) GO TO 30
 GO TO 25
 30 IF (.NOT. first1) GO TO 80
!*****
!    READ COLUMNS OF B INTO CORE.
!*****
 first1 = .false.
 IF (k > nk) GO TO 40
 k2 = k
 GO TO 50
 40 k2 = nk
 50 kt = iktbp - 1
 kb = ibcols - precn
 kbc = ibcid - 1
 DO  kk=1,k2
   kt = kt + 1
   kkk = iz(kt)
   CALL filpos (scr1,iz(kkk))
   kb = kb + precn
   CALL unpack(*55,scr1,z(kb))
   GO TO 59
   55 ib = kb - 1
   DO  l=1,precn
     ib = ib + 1
     z(ib) = 0.
   END DO
   59 CONTINUE
   kbc = kbc + 1
   iz(kbc) = kkk
 END DO
!*****
!    BEGIN CALCULATING MATRIX PRODUCT.
!*****
 80 kt = iktbp - 1
 kcount = 0
 preck = prec*k
 DO  ka=1,k
   kt = kt + 1
   kbc = ibcid - 1
   DO  kb=1,k2
     kbc = kbc + 1
     IF (iz(kt) == iz(kbc)) GO TO 100
   END DO
   CYCLE
   100 kkb = kb
   CALL mpy3p (z,z,z)
   iz(kt) = 0
   kcount = kcount + 1
   IF (first2 .OR. icore == 1) CYCLE
   i = iz(kbc)
   CALL mpy3nu (z)
   kbn = ibntu + kkb - 1
   iz(kbn) = ntbu
 END DO
!*****
!    SET RETURN FLAG.
!*****
 IF (kcount == k) GO TO 120
 iflag = 1
 GO TO 9999
 120 iflag = 0
 IF (icore /= 1 .OR. first2) GO TO 9999
 IF (j == m) GO TO 9999
 FILE = scr2
 CALL fwdrec(*902,scr2)
 GO TO 9999
!*****
!    ERROR MESSAGES.
!*****
 902 nerr = -2
 CALL mesage (nerr, FILE, NAME)
 
 9999 RETURN
END SUBROUTINE mpy3b
