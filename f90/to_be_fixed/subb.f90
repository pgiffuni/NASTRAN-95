SUBROUTINE subb(kb,ks,i,j,jb,lb,ls,ndy,nyfl,pi,eps,sgr,cgr,  &
        ar,beta,sum,ria,delx,yb,zb,ys,zs,x)
!   ***   COMPUTES ELEMENTS OF THE SUBMATRICES  DZP, DZZ, DZY, DYP,
!         DYZ  AND  DYY  USING  SUBROUTINE  DZY
 
 INTEGER, INTENT(IN OUT)                  :: kb
 INTEGER, INTENT(IN OUT)                  :: ks
 INTEGER, INTENT(IN)                      :: i
 INTEGER, INTENT(IN)                      :: j
 INTEGER, INTENT(IN OUT)                  :: jb
 INTEGER, INTENT(IN OUT)                  :: lb
 INTEGER, INTENT(IN OUT)                  :: ls
 INTEGER, INTENT(IN)                      :: ndy
 INTEGER, INTENT(IN OUT)                  :: nyfl
 REAL, INTENT(IN)                         :: pi
 REAL, INTENT(IN OUT)                     :: eps
 REAL, INTENT(IN OUT)                     :: sgr
 REAL, INTENT(IN OUT)                     :: cgr
 REAL, INTENT(IN)                         :: ar
 REAL, INTENT(IN OUT)                     :: beta
 COMPLEX, INTENT(OUT)                     :: sum
 REAL, INTENT(IN)                         :: ria(1)
 REAL, INTENT(IN)                         :: delx(1)
 REAL, INTENT(IN)                         :: yb(1)
 REAL, INTENT(IN)                         :: zb(1)
 REAL, INTENT(IN)                         :: ys(1)
 REAL, INTENT(IN)                         :: zs(1)
 REAL, INTENT(IN)                         :: x(1)
 REAL :: kd1r,kd1i, kd2r,kd2i
 COMPLEX :: dpur,dpul,dplr,dpll,dp
 
 COMMON /amgmn/ mcb(7),nrow,nd,NE,refc,fmach,kr
 COMMON     /kds/ ind,kd1r,kd1i,kd2r,kd2i
 
 flnd = FLOAT(nd)
 flne = FLOAT(NE)
 ind  = 0
 dpur = (0.0,0.0)
 dpul = (0.0,0.0)
 dplr = (0.0,0.0)
 dpll = (0.0,0.0)
 anot = ria(jb)
 dxs  = delx(j)
 absyb= ABS(yb(lb))
 abszb= ABS(zb(lb))
 iflag = 0
 idflag = 0
 IF (kb == 0)  GO TO  20
 test1= ABS(yb(lb) -yb(kb))
 test2= ABS(zb(lb) -zb(kb))
 IF  (test1 > eps. OR .test2 > eps)  GO TO  20
 iflag = 1
 IF(ndy /= nyfl) GO TO 20
 IF( i /= j ) GO TO 20
 idflag = 1
 d2d  =       1.0 /(2.0*pi*anot*anot*(1.0+ar))
 IF    (ndy /= 0)  d2d=d2d/ar
 sum  = CMPLX(d2d,0.0)
 sign1 = 1.0
 IF(ndy /= 0) sign1 = -1.0
 IF(absyb < eps) sum=(1.0+sign1*flnd)*sum
 IF(abszb < eps) sum=(1.0+sign1*flne)*sum
 dpur = sum
 20 CONTINUE
 xx   = x(i)
 y    = ys(ks)
 z    = zs(ks)
 xi1  = x(j) - 0.5*dxs
 xi2  = x(j) + 0.5*dxs
 eta  = ys(ls)
 zeta = zs(ls)
 ao   = anot
 idzdy= ndy
 igo  = 1
 lhs = 0
 IF(iflag == 1) GO TO 45
 30 CONTINUE
 CALL        dzy  (xx, y, z, sgr, cgr, xi1, xi2, eta, zeta, ar, ao,  &
     kr, refc, beta, fmach, lhs, idzdy ,   dzdyr ,   dzdyi )
 dp   = CMPLX(dzdyr,dzdyi)
 SELECT CASE ( igo )
   CASE (    1)
     GO TO 40
   CASE (    2)
     GO TO 50
   CASE (    3)
     GO TO 70
   CASE (    4)
     GO TO 80
 END SELECT
 40 CONTINUE
!  UPPER RIGHT-HAND SIDE CONTRIBUTION
 dpur = dp
 IF (kb == lb)  GO TO 100
 45 CONTINUE
 IF (nd == 0) GO TO 60
 IF (idflag == 1.AND.absyb < eps) GO TO 60
 igo  = 2
 eta  = -ys(ls)
 lhs = 1
 GO TO  30
 50 CONTINUE
!  UPPER LEFT-HAND SIDE CONTRIBUTION
 dpul = dp
 60 CONTINUE
 IF (NE == 0) GO TO 90
 IF(idflag == 1.AND.abszb < eps) GO TO 90
 igo  = 3
 eta  =  ys(ls)
 zeta = -zs(ls)
 lhs = 1
 GO TO  30
 70 CONTINUE
!  LOWER RIGHT-HAND SIDE CONTRIBUTION
 dplr = dp
 IF (nd == 0) GO TO 90
 IF(idflag == 1.AND.absyb < eps) GO TO 90
 igo  = 4
 eta  = -ys(ls)
 zeta = -zs(ls)
 lhs = 0
 GO TO  30
 80 CONTINUE
!  LOWER  LEFT-HAND SIDE CONTRIBUTION
 dpll = dp
 90 CONTINUE
 sum  = dpur + flnd*dpul + flne*dplr + flnd*flne*dpll
 100 CONTINUE
 RETURN
END SUBROUTINE subb
