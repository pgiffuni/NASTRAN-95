SUBROUTINE subpb (i,l,ls,j,sgr,cgr,yrec,zrec,sum,xic,delx,ee,xlam,  &
        sg,cg,ys,zs,nas,nasb,avr,zb,yb,arb,xle,xte,x,nb)
     
!     COMPUTES ELEMENTS OF THE SUBMATRICES  DPP, DPZ  AND  DPY
!     USING  SUBROUTINES  SNPDF, INCRO AND SUBI
 
 
 INTEGER, INTENT(IN OUT)                  :: i
 INTEGER, INTENT(IN OUT)                  :: l
 INTEGER, INTENT(IN OUT)                  :: ls
 INTEGER, INTENT(IN OUT)                  :: j
 REAL, INTENT(IN OUT)                     :: sgr
 REAL, INTENT(IN OUT)                     :: cgr
 REAL, INTENT(IN)                         :: yrec
 REAL, INTENT(IN)                         :: zrec
 COMPLEX, INTENT(OUT)                     :: sum
 REAL, INTENT(IN)                         :: xic(1)
 REAL, INTENT(IN)                         :: delx(1)
 REAL, INTENT(IN)                         :: ee(1)
 REAL, INTENT(IN)                         :: xlam(1)
 REAL, INTENT(IN)                         :: sg(1)
 REAL, INTENT(IN)                         :: cg(1)
 REAL, INTENT(IN)                         :: ys(1)
 REAL, INTENT(IN)                         :: zs(1)
 INTEGER, INTENT(IN)                      :: nas(1)
 INTEGER, INTENT(IN)                      :: nasb(1)
 REAL, INTENT(IN)                         :: avr(1)
 REAL, INTENT(IN)                         :: zb(1)
 REAL, INTENT(IN)                         :: yb(1)
 REAL, INTENT(IN)                         :: arb(1)
 REAL, INTENT(IN)                         :: xle(1)
 REAL, INTENT(IN)                         :: xte(1)
 REAL, INTENT(IN)                         :: x(1)
 INTEGER, INTENT(IN)                      :: nb
 REAL :: kr,m
 COMPLEX :: dpur,dpul,dplr,dpll,dp
 
 COMMON /amgmn/ mcb(7),nrow,nd,NE,refc,fmach,kr
 
 eps  = 0.00001
 m    = fmach
 beta = SQRT(1.0-m*m)
 fl   = refc
 flnd = FLOAT(nd)
 flne = FLOAT(NE)
 sgs  = sg(ls)
 cgs  = cg(ls)
 dpur = (0.0,0.0)
 dpul = (0.0,0.0)
 dplr = (0.0,0.0)
 dpll = (0.0,0.0)
 dij  = 0.0
 delr = 0.0
 deli = 0.0
 diji = 0.0
 delri= 0.0
 delii= 0.0
 infl = 0
 ioutfl = 0
 
!     UPPER RIGHT SENDING POINT
 
 igo  = 1
 tl   = xlam(j)
 sqtl = SQRT(1.0+tl**2)
 sl   = tl/sqtl
 cl   = 1.0/sqtl
 x0   = x(i) - xic(j)
 y0   = yrec - ys(ls)
 z0   = zrec - zs(ls)
 es   = ee(ls)
 dxs  = delx(j)
 ax   = x0
 ay   = y0
 az   = z0
 cv   = dxs
 
 30 nobi = 1
 na2  = 0
 CALL snpdf (sl,cl,tl,sgs,cgs,sgr,cgr,x0,y0,z0,es,dij,beta,cv)
 IF (kr <= eps) GO TO 40
 sdelx = dxs
 dely = 2.0*es
 ax1  = ax + es*tl
 ay1  = ay + es*cgs
 az1  = az + es*sgs
 ax2  = ax - es*tl
 ay2  = ay - es*cgs
 az2  = az - es*sgs
 CALL incro (ax,ay,az,ax1,ay1,az1,ax2,ay2,az2,sgr,cgr,sgs,cgs,  &
     kr,fl,beta,sdelx,dely,delr,deli)
 40 IF (nb == 0) GO TO 120
 noas = nas(l)
 
!     CHECK FOR ASSOCIATED BODIES
 
 IF (noas == 0) GO TO 120
 dijs  = dij
 delrs = delr
 delis = deli
 diji  = 0.0
 delri = 0.0
 delii = 0.0
 na1   = na2 + 1
 na2   = na2 + noas
 IF (na2 > nb) na2 = nb
 
!     START DO-LOOP FOR THE SUMMATION OF THE WING-IMAGE CONTRIBUTIONS
!     OVER  RANGE(P)
 
 DO  na = na1,na2
   nob  = nasb(na)
   
!     NOB IS THE SEQUENCE NUMBER OF THE CURRENT BODY ASSOCIATED WITH
!     PANEL  L  IN WHICH THE SENDING POINT  J  LIES
   
   nobi = nob
   da   = avr(nob)
   dar  = arb(nob)
   dxle = xle(nob)
   dxte = xte(nob)
   SELECT CASE ( igo )
     CASE (    1)
       GO TO 50
     CASE (    2)
       GO TO 60
     CASE (    3)
       GO TO 70
     CASE (    4)
       GO TO 80
   END SELECT
   50 CONTINUE
   dzb  = zb(nob)
   dyb  = yb(nob)
   deta = ys(ls)
   dzeta= zs(ls)
   GO TO  90
   60 CONTINUE
   dzb  = zb(nob)
   dyb  =-yb(nob)
   deta =-ys(ls)
   dzeta= zs(ls)
   GO TO  90
   70 CONTINUE
   dzb  =-zb(nob)
   dyb  = yb(nob)
   deta = ys(ls)
   dzeta=-zs(ls)
   GO TO  90
   80 CONTINUE
   dzb  =-zb(nob)
   dyb  =-yb(nob)
   deta =-ys(ls)
   dzeta=-zs(ls)
   90 CONTINUE
   dcgam= cgs
   dsgam= sgs
   dee  = es
   dxi  = xic(j)
   IF (dxi < dxle .OR. dxi > dxte) CYCLE
   CALL subi (da,dzb,dyb,dar,deta,dzeta,dcgam,dsgam,dee,dxi,tl,detai,  &
       dzetai,dcgami,dsgami,deei,dtlami,dmuy,dmuz,infl,ioutfl)
   dij  = 0.0
   IF (infl /= 0 .OR. ioutfl == 0) GO TO 100
   dtl  = dtlami
   dsqrtl = SQRT(1.0+dtl**2)
   dsl  = dtl/dsqrtl
   dcl  = 1.0/dsqrtl
   x0i  = x0
   y0i  = yrec - detai
   z0i  = zrec - dzetai
   CALL snpdf (dsl,dcl,dtl,dsgami,dcgami,sgr,cgr,x0i,y0i,z0i,deei,  &
       dij,beta,cv)
   diji = diji + dij
   IF  (kr <= eps) GO TO 100
   delr = 0.0
   deli = 0.0
   ayi  = y0i
   azi  = z0i
   ay1i = ayi - deei*dcgami
   az1i = azi - deei*dsgami
   ay2i = ayi + deei*dcgami
   az2i = azi + deei*dsgami
   deei2 = 2.0*deei
   CALL incro (ax,ayi,azi,ax1,ay1i,az1i,ax2,ay2i,az2i,sgr,cgr,dsgami,  &
       dcgami,kr,fl,beta,sdelx,deei2,delr,deli)
   delri = delri + delr
   delii = delii + deli
   CYCLE
   100 CONTINUE
   delri = 0.0
   delii = 0.0
 END DO
 dij  = dijs
 delr = delrs
 deli = delis
 120 CONTINUE
 dp = CMPLX(((dij+diji)-(delr+delri)),(-deli-delii))
 SELECT CASE ( igo )
   CASE (    1)
     GO TO 140
   CASE (    2)
     GO TO 150
   CASE (    3)
     GO TO 170
   CASE (    4)
     GO TO 180
 END SELECT
 140 CONTINUE
 dpur = dp
 IF (nd == 0) GO TO 160
 
!     UPPER LEFT  SENDING POINT
 
 igo = 2
 sgs =-sgs
 tl  =-tl
 sl  =-sl
 y0  = yrec + ys(ls)
 ay  = y0
 GO TO 30
 150 CONTINUE
 dpul = dp
 160 CONTINUE
 IF (NE == 0) GO TO 190
 
!     LOWER RIGHT SENDING POINT
 
 igo = 3
 tl  = xlam(j)
 sl  = tl/(SQRT(1.0+tl*tl))
 y0  = yrec - ys(ls)
 z0  = zrec + zs(ls)
 ay  = y0
 az  = z0
 sgs =-sg(ls)
 GO TO 30
 170 CONTINUE
 dplr = dp
 IF (nd == 0) GO TO 190
 
!     LOWER LEFT  SENDING POINT
 
 igo = 4
 sgs = sg(ls)
 tl  =-xlam(j)
 sl  = tl/(SQRT(1.0+tl*tl))
 y0  = yrec + ys(ls)
 ay  = y0
 GO TO 30
 180 CONTINUE
 dpll = dp
 190 CONTINUE
 sum  = dpur + flnd*dpul + flne*dplr + flnd*flne*dpll
 RETURN
END SUBROUTINE subpb
