SUBROUTINE alg06(r1,r2,x1,x2,h,s,vm,tb1,tb2,w,xk,sclfac,speed,&
                 spdfac,g,ej,hmin,nstrms,pi)
     
 
 REAL, INTENT(IN)                         :: r1(1)
 REAL, INTENT(IN)                         :: r2(1)
 REAL, INTENT(IN)                         :: x1(1)
 REAL, INTENT(IN OUT)                     :: x2(1)
 REAL, INTENT(IN)                         :: h(1)
 REAL, INTENT(IN)                         :: s(1)
 REAL, INTENT(IN)                         :: vm(1)
 REAL, INTENT(IN OUT)                     :: tb1(1)
 REAL, INTENT(IN)                         :: tb2(1)
 REAL, INTENT(OUT)                        :: w( 1)
 REAL, INTENT(IN)                         :: xk
 REAL, INTENT(IN OUT)                     :: sclfac
 REAL, INTENT(IN)                         :: speed
 REAL, INTENT(IN OUT)                     :: spdfac
 REAL, INTENT(IN)                         :: g
 REAL, INTENT(IN OUT)                     :: ej
 REAL, INTENT(IN)                         :: hmin
 INTEGER, INTENT(IN OUT)                  :: nstrms
 REAL, INTENT(IN)                         :: pi
 
 DIMENSION r(150),w2d(150),w3d(150),xx1(150),xx2(150),xx3(150),xx5( 9,9),b(150)
 
 EQUIVALENCE (xx2(1),xx5(1,1))
 
 ntub=nstrms-1
 DO  j=1,nstrms
   q1=h(j)-vm(j)**2*(1.0+(tb2(j)+r2(j)*speed*spdfac*pi&
      /(sclfac*30.0*vm(j)))**2)/(2.0*g*ej)
   IF(q1 < hmin)q1=hmin
   xx1(j)=alg4(q1,s(j))
   xx2(j)=alg5(q1,s(j))
 END DO
 CALL alg01(r2,xx1,nstrms,r2,q1,xx3,nstrms,0,1)
 DO  j=1,nstrms
   xx1(j)=xx3(j)*g/xx2(j)
 END DO
 q1=(r2(nstrms)-r2(1))/149.0
 r(1)=r2(1)
 DO  j=2,150
   r(j)=r(j-1)+q1
 END DO
 CALL alg01(r2,xx1,nstrms,r,xx2,q1,150,0,0)
 DO  j=1,nstrms
   xx3(j)=((r2(j)-r1(j))**2+(x2(j)-x1(j))**2)*(1.0+((tb1(j)+tb2(j))*0.5)**2)
 END DO
 CALL alg01(r2,xx3,nstrms,r,xx1,q1,150,0,0)
 DO  j=1,nstrms
   w2d(j)=vm(j)**2*(1.0+tb2(j)**2)
 END DO
 CALL alg01(r2,w2d,nstrms,r,xx3,q1,150,0,0)
 CALL alg01(r2,w  ,nstrms,r,w2d,q1,150,0,0)
 nkeep=nstrms
 nstrms=150
 ntub=149
 q2=(speed*spdfac*pi/(30.0*sclfac))**2
 DO  j=1,nstrms
   w3d(j)=0.0
 END DO
 b(1)=(r(2)-r(1))/2.0
 b(nstrms)=(r(nstrms)-r(ntub))/2.0
 DO  j=2,ntub
   b(j)=(r(j+1)-r(j-1))/2.0
 END DO
 DO  j=1,nstrms
   dr=xk*xx1(j)/xx3(j)*(q2*r(j)-xx2(j))
   IF(dr < 0.0) THEN
     GO TO   130
   ELSE IF (dr == 0.0) THEN
     GO TO   120
   ELSE
     GO TO   200
   END IF
   120   w3d(j)=w3d(j)+w2d(j)
   CYCLE
   130   IF(j == 1)GO TO 120
   IF(r(j)+dr <= r(1))GO TO 180
   DO  jj=2,j
     jjj=j-jj+1
     IF(r(j)+dr >= r(jjj))EXIT
   END DO
   150   jjj=jjj+1
   q1=w2d(j)*b(j)/(b(j)-dr)
   DO  jj=jjj,j
     w3d(jj)=w3d(jj)+q1
   END DO
   CYCLE
   180   a=b(j)*w2d(j)/(r(nstrms)-r(1))
   IF(j /= nstrms)a=b(j)*w2d(j)/((r(j+1)+r(j))*0.5-r(1))
   DO  jj=1,j
     w3d(jj)=w3d(jj)+a
   END DO
   CYCLE
   200   IF(j == nstrms)GO TO 120
   IF(r(j)+dr >= r(nstrms))GO TO 250
   DO  jj=j,nstrms
     IF(r(j)+dr < r(jj))EXIT
   END DO
   220   jj=jj-1
   q1=w2d(j)*b(j)/(b(j)+dr)
   DO  jjj=j,jj
     w3d(jjj)=w3d(jjj)+q1
   END DO
   CYCLE
   250   a=b(j)*w2d(j)/(r(nstrms)-r(1))
   IF(j /= 1)a=b(j)*w2d(j)/(r(nstrms)-(r(j)+r(j-1))*0.5)
   DO  jj=j,nstrms
     w3d(jj)=w3d(jj)+a
   END DO
 END DO
 nstrms=nkeep
 xx1(1)=0.0
 DO  ll=1,150
   xx1(1)=xx1(1)+w3d(ll)
 END DO
 DO  l=2,9
   xx1(l)=0.0
   DO  ll=1,150
     xx1(l)=xx1(l)+r(ll)**(l-1)*w3d(ll)
   END DO
 END DO
 DO  l=1,9
   DO  j=l,9
     IF(j == 1)GO TO 310
     xx5(l,j)=0.0
     DO  ll=1,150
       xx5(l,j)=xx5(l,j)+r(ll)**(l+j-2)
     END DO
     GO TO 320
     310   xx5(1,1)=150
     320   xx5(j,l)=xx5(l,j)
   END DO
 END DO
 CALL alg30(xx5,xx1)
 DO  j=1,nstrms
   w(j)=(((((((xx1(9)*r2(j)+xx1(8))*r2(j)+xx1(7))*r2(j)+xx1(6))*r2(j)  &
       +xx1(5))*r2(j)+xx1(4))*r2(j)+xx1(3))*r2(j)+xx1(2))*r2(j)+xx1(1)
 END DO
 
 RETURN
END SUBROUTINE alg06
