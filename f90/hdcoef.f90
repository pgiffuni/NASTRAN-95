SUBROUTINE hdcoef(x,y,z,xxx,jxx,ns,ccc,lz)
     
 
!     THIS SUBROUTINE DETERMINES EQUATION OF LINES AND PLANES.
 
 
 
 REAL, INTENT(IN)                         :: x(1)
 REAL, INTENT(IN)                         :: y(1)
 REAL, INTENT(IN)                         :: z(1)
 REAL, INTENT(OUT)                        :: xxx(1)
 INTEGER, INTENT(IN OUT)                  :: jxx
 INTEGER, INTENT(IN OUT)                  :: ns
 REAL, INTENT(OUT)                        :: ccc(1)
 INTEGER, INTENT(IN)                      :: lz
 INTEGER :: zcoef1,zcoef,ibcoef(5)
 DIMENSION  coe(8)
 COMMON/zzzzzz/rz(1)
 COMMON/hdptrs/xdum,xcc,xasolv,yasolv,zasolv,x1skt,y1skt,z1skt, zcoef1,zcoef
 COMMON/go3/l0,l1,l00,l01,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13
 DATA epsi / 1.0E-5 /
 
 LE=0
 ja=l13+(jxx-1)*lz
 jf=l12+(jxx-1)*5
 i=0
 j=1
 10 CONTINUE
 
 
!     SEARCH FOR MATCHING COORDINATES.
 
 
 i=i+1
 t=x(i+1)-x(i)
 s=y(i+1)-y(i)
 u=z(i+1)-z(i)
 IF (ABS(t) > epsi)  GO TO 20
 IF (ABS(s) > epsi)  GO TO 20
 IF (ABS(u) > epsi)  GO TO 20
 
 
!     MATCH FOUND.....PROCEED IF LIST IS NOT EXHAUSTED.
 
 
 i=i+2
 20 CONTINUE
 IF(i > ns)GO TO 70
 
 
!     DETERMINE EQUATION OF LINE-SEGMENTS.
 
 
 t=x(i+1)-x(i)
 t1=y(i+1)-y(i)
 IF ((ABS(t1) < epsi) .AND. (ABS(t) < epsi))  GO TO 10
 IF (ABS(t) > epsi)  GO TO 30
 29 CONTINUE
 ccc(j+ja)=0
 ccc(j+1+ja)=1
 ccc(j+2+ja)=-x(i)
 GO TO 40
 30 CONTINUE
 ccc(j+ja)=1
 e=(y(i+1)-y(i))/(x(i+1)-x(i))
 IF(ABS(e) > 100000.)GO TO 29
 f=(e*x(i))-y(i)
 ccc(j+1+ja)=-e
 ccc(j+2+ja)=f
 40 CONTINUE
 IF (ABS(ccc(j+ja)) > epsi)  GO TO 50
 ccc(j+3+ja)=y(i)
 ccc(j+4+ja)=y(i+1)
 GO TO 60
 50 CONTINUE
 ccc(j+3+ja)=x(i)
 ccc(j+4+ja)=x(i+1)
 60 CONTINUE
 j=j+5
 rz(zcoef1+LE)=z(i)
 rz(zcoef+LE)=z(i+1)
 LE=LE+1
 IF(LE > 3)GO TO 10
 ibcoef(LE)=i
 GO TO 10
 70 CONTINUE
 
 
!     DETERMINE EQUATION OF PLANE.
 
 j=(j-1)/5
 xxx(jf+5)=j
 IF(ns <= 3)GO TO 120
 k1=1
 k2=2
 k3=3
 a1=x(k3)-x(k1)
 b1=y(k3)-y(k1)
 c1=z(k3)-z(k1)
 a2=x(k2)-x(k1)
 b2=y(k2)-y(k1)
 c2=z(k2)-z(k1)
 coe(1)=b1*c2-b2*c1
 coe(2)=c1*a2-c2*a1
 coe(3)=a1*b2-a2*b1
 coe(4)=coe(1)*x(1)+coe(2)*y(1)+coe(3)*z(1)
 coe(4)=-coe(4)
 DO  j=1,4
   xxx(jf+j)=coe(j)
 END DO
 IF (ABS(coe(3)) > epsi)  GO TO 140
 j=1
 DO  k=1,LE
   japj=ja+j
   ccc(japj)=rz(zcoef1-1+k)
   ccc(japj+1)=rz(zcoef-1+k)
   j=j+5
 END DO
 IF (ABS(coe(1)) > epsi)  i=1
 IF (ABS(coe(2)) > epsi)  i=2
 p=coe(i)
 IF (ABS(p) < epsi)  p = epsi
 DO  k=1,4
   jfpk=jf+k
   xxx(jfpk)=xxx(jfpk)/p
 END DO
 GO TO 140
 120 CONTINUE
 xxx(jf+5)=1
 DO  ix=1,2
   xxx(jf+ix)=z(ix)
 END DO
 xxx(jf+3)=0
 140 CONTINUE
 RETURN
END SUBROUTINE hdcoef
