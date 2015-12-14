SUBROUTINE mpy3p (z,iz,dz)
!*****
!    PERFORMS MULTIPLICATION AND SUMMATION.
!*****
 
 REAL, INTENT(IN OUT)                     :: z(1)
 INTEGER, INTENT(IN)                      :: iz(1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: dz(1)
 DOUBLE PRECISION :: dfact
 
 
 
 INTEGER :: code,prec
 INTEGER :: zpntrs
 
 
 
 
 
 
! SUBROUTINE CALL PARAMETERS
 COMMON / mpy3cp / itrl,icore,n,ncb,m,dum1(3),zpntrs(22),laend, dum2(6),ka,kb
! FILES
 COMMON / mpy3tl / filea(7),fileb(7),filee(7),filec(7),scr1,scr2,  &
     scr,lkore,code,prec,lcore,scr3(7),buf1,buf2, buf3,buf4,e
 
 
 
 EQUIVALENCE     (fact,dfact)
! OPEN CORE POINTERS
 EQUIVALENCE     (ipoint,zpntrs(3)),      (iacols,zpntrs(5)),  &
     (itrans,zpntrs(7)),      (ic,zpntrs(9)),  &
     (ibcols,zpntrs(11)),     (iakj,zpntrs(21))
!*****
!    LOOP FOR ACCUMULATING SUMS.
!*****
 kj = iakj + ka - 1
 kj2 = (iakj - 1)/2 + ka
 kb = ibcols + prec*((kb - 1)*n - 1)
 IF (code == 2 .OR. icore == 1) GO TO 100
!*****
!    A(T)BA CASE.
!*****
 lp = ipoint - 1
 DO  l=1,n
! CALCULATE FACTOR = B(LK)*A(KJ) TO BE MULTIPLIED TO NON-ZERO TERMS IN
! LTH COLUMN OF A(T)
   kb = kb + prec
   lp = lp + 1
   IF (iz(lp) == 0) CYCLE
   IF (prec == 2) GO TO 10
   IF (z(kb) == 0.0) CYCLE
   fact = z(kb)*z(kj)
   GO TO 20
   10 kb2 = (kb + 1)/2
   IF (dz(kb2) == 0.0D0) CYCLE
   dfact = dz(kb2)*dz(kj2)
   20 i1 = iz(lp)
   IF (l == n) GO TO 40
! ACCUMULATE SUMS FOR NON-ZERO TERMS IN COLUMN L OF A(T)
   l1 = l + 1
   llp = lp
   DO  ll=l1,n
     llp = llp + 1
     IF (iz(llp) /= 0) GO TO 50
   END DO
   40 i2 = laend
   GO TO 60
   50 i2 = iz(llp) - 1
   60 iac = iacols + i1 - 2
   IF (prec == 2) GO TO 80
! SINGLE PRECISION CASE
   iat = itrans + i1 - 2
   DO  i=i1,i2
     iac = iac + 1
     iat = iat + 1
     ii = ic + iz(iac) - 1
     z(ii) = z(ii) + z(iat)*fact
   END DO
   CYCLE
! DOUBLE PRECISION CASE
   80 iat = (itrans - 3)/2 + i1
   DO  i=i1,i2
     iac = iac + 1
     iat = iat + 1
     ii = (ic - 1)/2 + iz(iac)
     dz(ii) = dz(ii) + dz(iat)*dfact
   END DO
   iii = (ic - 1)/2 + 1
 END DO
 GO TO 999
!*****
!    BA CASE.
!*****
 100 IF (prec == 2) GO TO 140
! SINGLE PRECISION CASE
 ii = ic - 1
 DO  i=1,n
   ii = ii + 1
   kb = kb + 1
   IF (z(kb) == 0.0) CYCLE
   z(ii) = z(ii) + z(kb)*z(kj)
 END DO
 GO TO 999
! DOUBLE PRECISION CASE
 140 ii = (ic - 1)/2
 kb = (kb + 1)/2
 DO  i=1,n
   ii = ii + 1
   kb = kb + 1
   IF (dz(kb) == 0.0D0) CYCLE
   dz(ii) = dz(ii) + dz(kb)*dz(kj2)
 END DO
 
 999 RETURN
END SUBROUTINE mpy3p
