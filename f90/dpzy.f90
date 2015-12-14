SUBROUTINE dpzy(   kb,iz,i,j1,j2,ifirst,ilast,yb,zb,  &
        avr,arb,th1a,th2a,nt121,nt122,nbaray,ncaray,  &
        nzykb,dpz,dpy)
!   ***   GENERATES ROWS OF THE SUBMATRICES  DPZ  AND DPY  USING
!         SUBROUTINE  SUBP
 
 INTEGER, INTENT(IN OUT)                  :: kb
 INTEGER, INTENT(OUT)                     :: iz
 INTEGER, INTENT(IN OUT)                  :: i
 INTEGER, INTENT(IN)                      :: j1
 INTEGER, INTENT(IN)                      :: j2
 INTEGER, INTENT(IN OUT)                  :: ifirst
 INTEGER, INTENT(IN OUT)                  :: ilast
 REAL, INTENT(IN)                         :: yb(1)
 REAL, INTENT(IN)                         :: zb(1)
 REAL, INTENT(IN)                         :: avr(1)
 REAL, INTENT(IN)                         :: arb(1)
 REAL, INTENT(IN)                         :: th1a(1)
 REAL, INTENT(IN)                         :: th2a(1)
 INTEGER, INTENT(IN)                      :: nt121(1)
 INTEGER, INTENT(IN)                      :: nt122(1)
 INTEGER, INTENT(IN)                      :: nbaray(1)
 INTEGER, INTENT(IN)                      :: ncaray(1)
 INTEGER, INTENT(IN OUT)                  :: nzykb
 COMPLEX, INTENT(OUT)                     :: dpz(1)
 COMPLEX, INTENT(OUT)                     :: dpy(1)
 INTEGER :: z
 COMPLEX :: sum
 
 
 COMMON /dlbdy/ nj1,nk1,np,nb,ntp,nbz,nby,ntz,nty,nt0,ntzs,ntys,  &
     inc,ins,inb,inas,izin,iyin,inbea1,inbea2,insbea,izb,iyb,  &
     iavr,iarb,infl,ixle,ixte,int121,int122,izs,iys,ics,iee,isg,  &
     icg,ixij,ia,idelx,ixic,ixlam,ia0,ixis1,ixis2,ia0p,iria  &
     ,inasb,ifla1,ifla2,ith1a,ith2a, ecore,next,scr1,scr2,scr3,scr4,scr5
 COMMON /zzzzzz / z(1)
 
 pi   = 3.1415926
 ix1    = 1
 iz   = iz+1
!  IZ  IS THE BODY-ELEMENT NUMBER FOR BODY  KB  --  IZ RUNS FROM  1
!  THROUGH  NBE-SUB-KB
 ix2 = nt122(kb)
 IF (iz >= ifirst.AND.iz <= ilast)  ix2=nt121(kb)
 DO    ix=ix1,ix2
   l    = 1
   ksp = 0
!  L IS THE PANEL NUMBER ASSOCIATED WITH SENDING   POINT  J
   ls   = 1
!  LS IS THE STRIP NUMBER ASSOCIATED WITH SENDING   POINT  J
   nbxs = nbaray(l)
   nc1  = ncaray(l)
   nbcum= nc1
   ixp1 = ix+1
   IF (ixp1 > ix2)  ixp1=ix1
   ixm1 = ix-1
   IF (ixm1 == 0)   ixm1=ix2
   IF (iz >= ifirst.AND.iz <= ilast)  GO TO  30
   theta= th2a(ix)
   thp1 = th2a(ixp1)
   thm1 = th2a(ixm1)
   GO TO  40
   30 CONTINUE
   theta = th1a(ix)
   thp1 = th1a(ixp1)
   thm1 = th1a(ixm1)
   40 CONTINUE
   IF (ix == ix1)  thm1=thm1-2.0*pi
   IF (ix == ix2)  thp1=thp1+2.0*pi
   delth= 0.5*(thp1 - thm1)
   yrec = yb(kb)+avr(kb)*COS(theta)
   zrec = zb(kb)+avr(kb)*arb(kb)*SIN(theta)
   rho  = SQRT(1.0+(arb(kb)**2 - 1.0) * (COS(theta))**2)
   sgr  = -arb(kb)*COS(theta)/rho
   cgr  = SIN(theta)/rho
   smult= SIN(theta) * rho / pi
   cmult= COS(theta) * rho / pi
   DO    j=j1,j2
     CALL subpb(i,l,ls,j,sgr,cgr,yrec,zrec,sum,z(ixic),z(idelx),z(iee)  &
         ,z(ixlam),z(isg),z(icg),z(iys),z(izs),z(inas),z(inasb+ksp),  &
         z(iavr),z(izb),z(iyb),z(iarb),z(ixle),z(ixte),z(ia),nb)
     SELECT CASE ( nzykb )
       CASE (    1)
         GO TO 50
       CASE (    2)
         GO TO 50
       CASE (    3)
         GO TO 60
     END SELECT
     50 CONTINUE
     dpz(j) = dpz(j) + sum * smult * delth
     IF (nzykb == 1)  GO TO 70
     60 CONTINUE
     dpy(j) = dpy(j) + sum * cmult * delth
     70 CONTINUE
     IF (j == j2)  CYCLE
     IF (j < nbxs)   GO TO  80
     ksp = ksp + z(inas+l-1)
     l    = l+1
     nc1  = ncaray(l)
     nbxs = nbaray(l)
     80 CONTINUE
     IF (j < nbcum)  CYCLE
     ls   = ls+1
     nbcum= nbcum+nc1
   END DO
 END DO
 RETURN
END SUBROUTINE dpzy
