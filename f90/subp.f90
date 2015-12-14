SUBROUTINE subp (i,l,ls,j,sgr,cgr,yrec,zrec,sum,xic,delx,ee,xlam,  &
        sg,cg,ys,zs)
     
!     COMPUTES ELEMENTS OF THE SUBMATRICES  DPP, DPZ  AND  DPY
!     USING  SUBROUTINES  SNPDF,  INCRO  AND SUBI
 
 
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
 REAL :: kr,m
 COMPLEX :: dpur,dpul,dplr,dpll,dp
 
 COMMON /amgmn/ mcb(7),nrow,nd,NE,refc,fmach,kr
 COMMON /dlcom/ dum(3),f
 
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
 
!     UPPER RIGHT SENDING POINT
 
 igo  = 1
 tl   = xlam(j)
 sqtl = SQRT(1.0+tl**2)
 sl   = tl/sqtl
 cl   = 1.0/sqtl
 x    = xic(i) + f*delx(i)
 x0   = x - xic(j)
 y0   = yrec - ys(ls)
 z0   = zrec - zs(ls)
 es   = ee(ls)
 dxs  = delx(j)
 ax   = x0
 ay   = y0
 az   = z0
 cv   = dxs
 
 30 nobi = 1
 CALL snpdf (sl,cl,tl,sgs,cgs,sgr,cgr,x0,y0,z0,es,dij,beta,cv)
 IF (kr <= eps) GO TO  40
 sdelx= dxs
 dely = 2.0*es
 ax1  = ax + es*tl
 ay1  = ay + es*cgs
 az1  = az + es*sgs
 ax2  = ax - es*tl
 ay2  = ay - es*cgs
 az2  = az - es*sgs
 CALL incro (ax,ay,az,ax1,ay1,az1,ax2,ay2,az2,sgr,cgr,sgs,cgs,  &
     kr,fl,beta,sdelx,dely,delr,deli)
 40 CONTINUE
 dp = CMPLX(((dij+diji)-(delr +delri)),(-deli-delii))
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
 
!     TEST FOR ABS(YS(LS)) .LE..001 TAKEN OUT
 
 IF (nd == 0) GO TO 160
 
!     UPPER LEFT  SENDING POINT
 
 igo  = 2
 sgs  =-sgs
 tl   =-tl
 sl   =-sl
 y0   = yrec + ys(ls)
 ay   = y0
 GO TO  30
 150 CONTINUE
 dpul = dp
 160 CONTINUE
 IF (NE == 0) GO TO 190
 
!     LOWER RIGHT SENDING POINT
 
 igo  = 3
 tl   = xlam(j)
 sl   = tl/(SQRT(1.0+tl*tl))
 y0   = yrec - ys(ls)
 z0   = zrec + zs(ls)
 ay   = y0
 az   = z0
 sgs  =-sg(ls)
 GO TO  30
 170 CONTINUE
 dplr = dp
 IF (nd == 0) GO TO 190
 
!     LOWER LEFT  SENDING POINT
 
 igo  = 4
 sgs  = sg(ls)
 tl   =-xlam(j)
 sl   = tl/(SQRT(1.0+tl*tl))
 y0   = yrec + ys(ls)
 ay   = y0
 GO TO  30
 180 CONTINUE
 dpll = dp
 190 CONTINUE
 sum  = dpur + flnd*dpul + flne*dplr + flnd*flne*dpll
 RETURN
END SUBROUTINE subp
