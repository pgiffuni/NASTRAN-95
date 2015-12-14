SUBROUTINE fvrs1c (z,w1,omega,nf)
!----------------------------------------------------------------------
 
 COMPLEX, INTENT(OUT)                     :: z(3,nf)
 REAL, INTENT(IN)                         :: w1(nf)
 REAL, INTENT(IN)                         :: omega
 INTEGER, INTENT(IN)                      :: nf
 COMPLEX :: p,z1,py1a,py1b,pz1a,pz1b,py2a,py2b,pz2a,pz2b
 
 
 
 COMMON /condas/ pi,twopi,radeg,degra,s4piq
 COMMON /BLANK / dum(5),ixt,ixp,iyt,iyp,izt,izp,dum1(3)
!----------------------------------------------------------------------
 ll=1
 DO  kkk=1,nf
   IF(w1(kkk) == 0.0)GO TO 30
   a=1.0
   IF(w1(kkk)-omega < 0.0)a=-1.0
   b=1.0
   IF(w1(kkk)+omega < 0.0)b=-1.0
   
!     COMPUTE BASE(FI)(3X3)  IF W.NE.0--MODFRL=TRUE
   
!     ZERO OUT MATRIX
   
   kk=ll+2
   DO  i=1,3
     DO  j=ll,kk
       z(i,j)=(0.0,0.0)
     END DO
   END DO
   f=w1(kkk)/twopi
   IF(ixt == -1) GO TO 9
   CALL tab(ixt,f,xo)
   IF(ixp == -1)GO TO 4
   CALL tab(ixp,f,phi)
   rad=phi*degra
   z1=CMPLX(0.0,rad)
   p=CEXP(z1)
   GO TO 5
   4    p=(1.0,0.0)
   5    z(1,ll+1)=xo*p
   9    IF(iyt == -1)GO TO 11
   CALL tab(iyt,f,yo)
   GO TO 12
   11   yo=0.0
   12   IF(izt == -1)GO TO 13
   CALL tab(izt,f,zo)
   GO TO 14
   13   zo=0.0
   14   IF(iyp == -1)GO TO 15
   CALL tab(iyp,f,phi)
   GO TO 16
   15   phi=0.0
   16   rad=phi*degra
   z1=CMPLX(0.0,rad)
   py1a=CEXP(a*z1)
   py1b=CEXP(b*z1)
   z1=CMPLX(0.0,rad-0.5*pi*a)
   py2a=CEXP(a*z1)
   z1=CMPLX(0.0,rad-0.5*pi*b)
   py2b=CEXP(b*z1)
   IF(izp == -1)GO TO 17
   CALL tab(izp,f,phi)
   GO TO 18
   17   phi=0.0
   18   rad=phi*degra
   z1=CMPLX(0.0,rad)
   pz1a=CEXP(a*z1)
   pz1b=CEXP(b*z1)
   z1=CMPLX(0.0,rad-0.5*pi*a)
   pz2a=CEXP(a*z1)
   z1=CMPLX(0.0,rad-0.5*pi*b)
   pz2b=CEXP(b*z1)
   z(2,ll)=(yo*py1a-a*zo*pz2a)*0.5
   z(3,ll)=(a*yo*py2a+zo*pz1a)*0.5
   z(2,ll+2)=(yo*py1b+b*zo*pz2b)*0.5
   z(3,ll+2)=(-b*yo*py2b+zo*pz1b)*0.5
   ll=ll+3
   CYCLE
   30   CONTINUE
   
!     COMPUTE BASE(FI)(3X2) IF W1=0.0, FOR MODFRL=TRUE
   
   a=1.0
   IF(omega < 0.0)a=-1.0
!------ZERO OUT MATRIX(3X2)
   kk=ll+1
   DO  i=1,3
     DO  j=ll,kk
       z(i,j)=(0.0,0.0)
     END DO
   END DO
   f=w1(kkk)/twopi
   IF(ixt == -1)GO TO 90
   CALL tab(ixt,f,xo)
   IF(ixp == -1)GO TO 40
   CALL tab(ixp,f,phi)
   rad=phi*degra
   z1=CMPLX(0.0,rad)
   p=CEXP(z1)
   GO TO 50
   40   p=(1.0,0.0)
   50   z(1,ll)=xo*p
   GO TO 100
   90   z(1,ll)=(0.0,0.0)
   100  IF(iyt == -1)GO TO 190
   CALL tab(iyt,f,yo)
   IF(iyp == -1) GO TO 140
   CALL tab(iyp,f,phi)
   rad=phi*degra
   cy=COS(rad)
   GO TO 150
   140  cy=1.0
   150  yy=yo*cy
   GO TO 200
   190  yy=0.0
   200  IF(izt == -1)GO TO 290
   CALL tab(izt,f,zo)
   IF(izp == -1) GO TO 240
   CALL tab(izp,f,phi)
   rad=phi*degra
   cz=COS(rad)
   GO TO 250
   240  cz=1.0
   250  zz=zo*cz
   GO TO 300
   290  zz=0.0
   300  z(2,kk)=yy-a*CMPLX(0.0,zz)
   z(3,kk)=zz+a*CMPLX(0.0,yy)
   ll=ll+2
 END DO
 RETURN
END SUBROUTINE fvrs1c
