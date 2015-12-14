SUBROUTINE mpy3c (z,iz,dz)
!*****
!    PERFORMS MULTIPLICATION AND SUMMATION FOR REMAINING TERMS OF COLUMN
!    OF A.
!*****
 
 REAL, INTENT(OUT)                        :: z(1)
 INTEGER, INTENT(IN OUT)                  :: iz(1)
 REAL, INTENT(IN OUT)                     :: dz
 INTEGER :: prec
 INTEGER :: scr1,scr2,FILE
 INTEGER :: zpntrs
 INTEGER :: utyp,urow1,urown,uincr
 INTEGER :: precn
 
 
 
 LOGICAL :: first2
 
 
 
 
!     DIMENSION NAME(2)
 
 
! FILES
 COMMON / mpy3tl / filea(7),fileb(7),filee(7),filec(7),scr1,scr2,  &
     scr,lkore,code,prec,lcore,scr3(7),buf1,buf2, buf3,buf4,e
! SUBROUTINE CALL PARAMETERS
 COMMON / mpy3cp / itrl,icore,n,ncb,m,nk,dum1(2),zpntrs(22),  &
     dum2(2),first2,k,k2,kcount,iflag,ka,ltbc,j, ltac,ntbu
! UNPACK
 COMMON / unpakx / utyp,urow1,urown,uincr
 
 
 
 EQUIVALENCE     (ibcols,zpntrs(11)),     (ibcid,zpntrs(13)),  &
     (ibntu,zpntrs(15)),      (iktbp,zpntrs(17)), (iantu,zpntrs(19))
 
 
 
!     DATA NAME / 4HMPY3,4HC    /
!*****
!    INITIALIZATION.
!*****
 utyp  = prec
 urow1 = 1
 urown = n
 uincr = 1
 precn = prec*n
 FILE  = scr1
!*****
!    TEST TO SEE IF LESS THAN NK COLUMNS OF B IN CORE.
!*****
 IF (first2) GO TO 30
!*****
!    DETERMINE WHICH COLUMN OF B TO BE PUT INTO CORE.
!*****
 lta = 0
 ia  = iantu - 1
 DO  i=1,k
   ia = ia + 1
   IF (lta >= iz(ia)) CYCLE
   lta = iz(ia)
   ik  = iktbp + i - 1
   ltac= iz(ik)
   ka  = i
 END DO
!*****
!    DETERMINE WHICH COLUMN OF B TO BE REPLACED.
!*****
 ltb = 0
 ib  = ibntu - 1
 DO  i=1,nk
   ib = ib + 1
   IF (ltb >= iz(ib)) CYCLE
   ltb  = iz(ib)
   ltbc = i
 END DO
 GO TO 50
!*****
!    LESS THAN NK COLUMNS OF B IN CORE.
!*****
 30 k2  = k2 + 1
 ltbc= k2
 kk  = iktbp - 1
 DO  ka=1,k
   kk  = kk + 1
   IF (iz(kk) == 0) CYCLE
   ltac = iz(kk)
   EXIT
 END DO
!*****
!    ADD OR REPLACE COLUMN OF B INTO CORE.
!*****
 50 CALL filpos (scr1,iz(ltac))
 kk = ibcols + precn*(ltbc - 1)
 CALL unpack(*55,scr1,z(kk))
 GO TO 59
 55 ik = kk - 1
 DO  l=1,precn
   ik = ik + 1
   z(ik) = 0.
 END DO
 59 CONTINUE
 IF (first2) GO TO 70
 IF (icore == 1) GO TO 60
 CALL mpy3nu (z)
 kk = ibntu + ltbc - 1
 iz(kk) = ntbu
 60 kk = iantu + ka - 1
 iz(kk) = 0
 70 kk = ibcid + ltbc - 1
 iz(kk) = ltac
 kb = ltbc
!*****
!    PERFORM COMPUTATION.
!*****
 CALL mpy3p (z,z,z)
 ltbc = kb
 kk = iktbp + ka - 1
 iz(kk) = 0
 kcount = kcount + 1
 RETURN
END SUBROUTINE mpy3c
