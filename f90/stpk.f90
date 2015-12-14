SUBROUTINE stpk(ek,n,nstop,nopen,nsted,tsr,pm,cr,ci,im,j1)
!     COMPUTES K MATRIX FOR STRIP NUMBER N
!     EK= LOCAL REDUCED FREQUENCY
!     NSTOP  =2 FOR NO CONTROL SURFACE
!     NOPEN  =1 FOR OPEN GAP
!     TSR = GAP/SEMICHORD RATIO  (FOR CLOSED STAGE ONLY)
!     NSTED =1 FOR STEADY CASE
 
 REAL, INTENT(IN)                         :: ek
 INTEGER, INTENT(IN OUT)                  :: n
 INTEGER, INTENT(IN OUT)                  :: nstop
 INTEGER, INTENT(IN OUT)                  :: nopen
 INTEGER, INTENT(IN OUT)                  :: nsted
 REAL, INTENT(IN OUT)                     :: tsr
 REAL, INTENT(IN)                         :: pm(1)
 REAL, INTENT(OUT)                        :: cr
 REAL, INTENT(OUT)                        :: ci
 INTEGER, INTENT(IN OUT)                  :: im
 INTEGER, INTENT(IN OUT)                  :: j1
 DIMENSION p(37)
 COMPLEX :: ekm,w,t,r,v1,v2,UNIT,cmp0
 COMPLEX :: w2
 COMMON /stripc/nns,bref,clam,fm,ncirc,nncirc,ekr(1),  &
     dum,       bb(4),beta(4),ekm(4,4)
 DATA nhek,nhtr,nhti,nhsize /4HEK  ,4HTR  ,4HTI  ,4HSIZE/
 
 UNIT=CMPLX(1.0,0.0)
 cmp0=CMPLX(0.0,0.0)
 t=2.0*UNIT
 DO  i=1,37
   p(i)=pm(i)
 END DO
 DO  i=1,4
   DO  j=1,4
     ekm(i,j)= cmp0
   END DO
 END DO
 a1=0.318310
 a2=0.101321
 IF(nsted /= 1) GO TO 50
!     STEADY CASE
 e1k = 1.e20
 ekm(1,2)=2.0
 IF(nstop == 2) GO TO 100
 ekm(1,3)=a1*2.0*p(1)
 ekm(2,3)=a1*p(5)
 ekm(3,2)=a1*2.0*p(31)
 ekm(3,3)=a2*(2.0*p(1)*p(31) + p(35) )
 ekm(4,2)=a1*p(8)
 ekm(4,3)=a2*(p(1)*p(8) + p(10) )
 IF(nopen == 1) GO TO 100
!     CLOSED STAGE
 ekm(1,4)=a1*2.0*p(13)
 ekm(2,4)=a1*p(15)
 tst=AMAX1(0.01  ,tsr)
 ekm(3,4)=a2*(2.0*p(13)*p(31) +2.0*ALOG(tst) +p(21))
 ekm(4,4)=a2*(p(13)*p(8) + p(18))
 50 IF(nsted == 1) GO TO 100
 
!     UNSTEADY CASE, EM(1,1)=(K SUB A)/EK**2, ETC.
 e1k  = 1./ek
 t=cmp0
 v1=cmp0
 v2=cmp0
 IF(ek > 1000.0) GO TO 71
 IF ( ncirc > 0 ) GO TO 73
 CALL stpbs0(ek,1,bj0,by0)
 CALL stpbs1(ek,1,bj1,by1)
 denom=(bj1+by0)**2 + (by1-bj0)**2
 cr= (bj1*(bj1+by0) + by1*(by1-bj0) )/denom
 ci=-(by1*by0 + bj1*bj0)/denom
!     (CR + I*CI  = THEODORSEN FUNCTION)
 GO TO 72
!  NEXT 8 STATEMENTS ARE FOR GENERATION OF WAGNER FUNCTIONS
 73 cr = bb(1)
 ci = 0.0
 DO  nn = 2,nncirc
   beoek = beta(nn)/ek
   fcr = bb(nn)/(1.0 + beoek*beoek)
   cr = cr + fcr
   ci = ci + fcr*beoek
 END DO
 72    t=2.0*CMPLX(cr,ci) - UNIT
 w=CMPLX(0.0,ek)
 v1=UNIT/w
 v2=v1*v1
 60 r=t+UNIT
 w2 = -w*w
 ekm(1,1)=-( r*v1 +1. )
 ekm(1,1) = ekm(1,1) * w2
 ekm(1,2)=-( r*(v2+v1) + v1 + 0.5 )
 ekm(1,2) = ekm(1,2) * w2
 ekm(2,1)=-( 0.5 )
 ekm(2,1) = ekm(2,1) * w2
 ekm(2,2)=-( v1 + 0.375 )
 ekm(2,2) = ekm(2,2) * w2
 IF(nstop == 2) GO TO 100
 ekm(1,3)=-a1*( r*(v2*p(1)+0.5*v1*p(2)) + v1*p(3) + 0.5 *p(4) )
 ekm(1,3) = ekm(1,3) * w2
 ekm(2,3)=-a1*( v2*p(5) + 0.5*v1*p(6) + 0.25*p(7) )
 ekm(2,3) = ekm(2,3) * w2
 ekm(3,1)=-a1*( r*v1*p(31) + p(3) )
 ekm(3,1) = ekm(3,1) * w2
 ekm(3,2)=-a1*( r*(v2+v1)*p(31) + v1*p(32) + 0.25*p(6) )
 ekm(3,2) = ekm(3,2) * w2
 ekm(3,3)=-a2*( r*(v2*p(1)+0.5*v1*p(2))*p(31) +  &
     v2*p(35) + v1*p(36) + 0.5*p(37) )
 ekm(3,3) = ekm(3,3) * w2
 ekm(4,1)=-a1*0.5*( r*v1*p(8) + p(4) )
 ekm(4,1) = ekm(4,1) * w2
 ekm(4,2)=-a1*0.5*( r*(v2+v1)*p(8) + v1*p(9) + 0.5*p(7) )
 ekm(4,2) = ekm(4,2) * w2
 ekm(4,3)=-a2*( r*(v2*p(1)+0.5*v1*p(2))*0.5*p(8) +  &
     v2*p(10) + 0.5*v1*p(11) + 0.25*p(12) )
 ekm(4,3) = ekm(4,3) * w2
 IF(nopen /= 1) GO TO 70
!     OPEN STAGE
 ekm(1,4)=-a1*( r*v1*p(1) + p(3) )
 ekm(1,4) = ekm(1,4) * w2
 ekm(2,4)=-a1*( v1*p(5) + 0.25*p(6) )
 ekm(2,4) = ekm(2,4) * w2
 ekm(3,4)=-a2*( r*v1*p(1)*p(31) + v1*p(35) + p(17) )
 ekm(3,4) = ekm(3,4) * w2
 ekm(4,4)=-a2*( r*0.5*v1*p(1)*p(8) + v1*p(10) + 0.5*p(37) )
 ekm(4,4) = ekm(4,4) * w2
 GO TO 100
 70 CONTINUE
!     CLOSED STAGE
 ekm(1,4)=-a1*( r*(v2*p(13)+ v1*p(1) ) + v1*p(14) + p(3) )
 ekm(1,4) = ekm(1,4) * w2
 ekm(2,4)=-a1*( v2*p(15) + 2.0*v1*p(5) + 0.25*p(6) )
 ekm(2,4) = ekm(2,4) * w2
 tst=AMAX1(0.01  ,tsr)
 ekm(3,4)=-a2*( r*(v2*p(13)+v1*p(1))*p(31) +  &
     v2*(2.0*ALOG(tst) + p(21)) + v1*p(16) + p(17) )
 ekm(3,4) = ekm(3,4) * w2
 ekm(4,4)=-a2*( r*(v2*p(13)+v1*p(1))*0.5*p(8) +  &
     v2*p(18) + v1*p(19) + 0.5*p(37) )
 ekm(4,4) = ekm(4,4) * w2
 100 CONTINUE
 CALL bug(nhek  ,100,ek,1)
 CALL bug(nhtr  ,100,cr,1)
 CALL bug(nhti  ,100,ci,1)
 CALL bug(nhsize,100,n,1)
 RETURN
 71    cr = .5
 ci = 0.
 w = CMPLX(0.0,ek)
 GO TO 60
END SUBROUTINE stpk
