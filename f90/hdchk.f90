SUBROUTINE hdchk(xxx,ccc,nno,ii,xi,yi,ngx,zm,zmi,  &
        rv,rvi,tgm,tgi,zi,lz,xcc)
     
 
 
!     THIS SUBROUTINE SOLVES FOR THE POINTS OF INTERSECTION ON THE
!     LINES OF THE JTH ELEMENT WITH OTHER LINES AND PLANES(RELEVANT)
 
 
 
 REAL, INTENT(IN)                         :: xxx(1)
 REAL, INTENT(IN)                         :: ccc(1)
 INTEGER, INTENT(IN OUT)                  :: nno(1)
 INTEGER, INTENT(IN)                      :: ii
 REAL, INTENT(OUT)                        :: xi(1)
 REAL, INTENT(OUT)                        :: yi(1)
 INTEGER, INTENT(OUT)                     :: ngx(1)
 REAL, INTENT(IN)                         :: zm(1)
 REAL, INTENT(IN)                         :: zmi(1)
 REAL, INTENT(IN)                         :: rv(1)
 REAL, INTENT(IN)                         :: rvi(1)
 REAL, INTENT(IN)                         :: tgm(1)
 REAL, INTENT(IN)                         :: tgi(1)
 REAL, INTENT(OUT)                        :: zi(1)
 INTEGER, INTENT(IN)                      :: lz
 REAL, INTENT(IN)                         :: xcc(1)
 
 
 COMMON/hedg/js,m,jt,vx,vx1,vx2,vx3,nn
 COMMON/go3/l0,l1,l00,l01,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13
 
 jm=1
 eex=.015
 EXP=.005
 ngx(1)=0
 IF(ii == 0)GO TO 190
 IF(nn == 1) GO TO 5
 IF(vx3 /= 0.)GO TO 5
 a=xxx(jt+2)
 b=xxx(jt+1)
 c=xxx(jt+4)
 z1=xcc(js)
 z2=xcc(js+1)
 IF(a == 0.) GO TO 1
 y1=-xcc(js+3)*b-c
 y2=-xcc(js+4)*b-c
 x1=xcc(js+3)
 x2=xcc(js+4)
 GO TO 50
 1 CONTINUE
 y1=xcc(js+3)
 y2=xcc(js+4)
 x1=-c
 x2=x1
 GO TO 50
 5 CONTINUE
 a=xcc(js)
 b=xcc(js+1)
 c=xcc(js+2)
 IF(a == 0.)GO TO 20
 y1=-xcc(js+3)*xcc(js+1)-xcc(js+2)
 y2=-xcc(js+4)*xcc(js+1)-xcc(js+2)
 x1=xcc(js+3)
 x2=xcc(js+4)
 GO TO 30
 20 CONTINUE
 y1=xcc(js+3)
 y2=xcc(js+4)
 x1=-xcc(js+2)
 x2=x1
 30 CONTINUE
 IF(nn /= 1)GO TO 40
 z1=xxx(1+jt)
 z2=xxx(2+jt)
 GO TO 50
 40 CONTINUE
 z1=-(vx+vx1*y1+vx2*x1)/vx3
 z2=-(vx+vx1*y2+vx2*x2)/vx3
 50 CONTINUE
 al=x2-x1
 bl=y2-y1
 cl=z2-z1
 eg=AMIN1(z1,z2)
 egx=AMAX1(x1,x2)
 egx1=AMIN1(x1,x2)
 egy=AMAX1(y1,y2)
 egy1=AMIN1(y1,y2)
 
 
!     THIS CODE DETERMINES THE POINTS OF INTERSECTIONS ON THE LINES OF
!     JTH ELEMENT RESULTING FROM THE INTERSECTION OF THE PLANES WITH
!     THESE LINES.
 
 
 DO  jr=1,ii
   lg=nno(l4+jr)
   nno(l4+jr)=IABS(nno(jr+l4))
   LE=nno(l4+jr)
   je=l13+lz*(LE-1)
   ju=l12+5*(LE-1)
   nk=xxx(5+ju)
   jv=1
   ac=xxx(1+ju)
   bc=xxx(2+ju)
   cc=xxx(3+ju)
   d=xxx(4+ju)
   IF(egx < tgm(l5+LE))CYCLE
   IF(egx1 > tgi(l6+LE))CYCLE
   IF(egy < rvi(l8+LE))CYCLE
   IF(egy1 > rv(l7+LE))CYCLE
   IF(eg > zm(l2+LE))CYCLE
   IF(lg < 0)GO TO 80
   IF((al == 0.).AND.(bl == 0.))GO TO 80
   IF(al == 0.)GO TO 60
   xp=((bc*bl)/al)*x1+(cc*cl/al)*x1-d
   xp=xp-bc*y1-cc*z1
   vu=ac+(bc*bl/al)+(cc*cl/al)
   IF(vu == 0.)GO TO 80
   xp=xp/vu
   t=(xp-x1)/al
   yp=t*bl+y1
   zp=t*cl+z1
   GO TO 70
   60 CONTINUE
   yp=(ac*al/bl)*y1+(cc*cl/bl)*y1-d
   yp=yp-cc*z1-ac*x1
   vu=bc+(ac*al/bl)+(cc*cl/bl)
   IF(vu == 0.)GO TO 80
   yp=yp/vu
   t=(yp-y1)/bl
   xp=t*al+x1
   zp=t*cl+z1
   70 CONTINUE
   s=zp-zm(l2+LE)
   s1=zp-zmi(l3+LE)
   IF((ABS(s) < eex).OR.(ABS(s1) < eex))GO TO 56
   IF(s*s1 > 0.)GO TO 80
   56 CONTINUE
   s=xp-tgm(l5+LE)
   s1=xp-tgi(l6+LE)
   IF(s*s1 > 0.)GO TO 80
   s=yp-rv(l7+LE)
   s1=yp-rvi(l8+LE)
   IF(s*s1 > 0.)GO TO 80
   t=xp
   IF(a == 0.)t=yp
   s=t-xcc(js+3)
   s1=t-xcc(js+4)
   IF(s*s1 >= 0.)GO TO 80
   m=m+1
   
!     STORES INTERSECTIONS.
   
   xi(m+1)=xp
   yi(m+1)=yp
   zi(m+1)=zp
   80 CONTINUE
   
!     THIS CODE DETERMINES INTERSECTION POINTS OF LINES WITH LINES.
   
   DO  jc=1,nk
     b1=ccc(jv+1+je)
     a1=ccc(jv+je)
     c1=ccc(jv+2+je)
     t=a1*b-b1*a
     IF(t == 0.)GO TO 160
     xo=(c1*a-c*a1)/t
     IF((ABS(b) <= 50.).AND.(a /= 0.))GO TO 90
     yo=-c1-b1*xo
     GO TO 100
     90 CONTINUE
     yo=-c-b*xo
     100 CONTINUE
     t=xo
     IF(a == 0.)t=yo
     s=t-xcc(js+3)
     s1=t-xcc(js+4)
     IF(s*s1 >= 0.)GO TO 160
     t=xo
     IF(a1 == 0.)t=yo
     s1=t-ccc(jv+4+je)
     s=t-ccc(jv+3+je)
     IF((ABS(s) <= eex).OR.(ABS(s1) <= eex))GO TO 110
     IF(s*s1 > 0.)GO TO 160
     110 CONTINUE
     IF(cc == 0.)GO TO 160
     zx=-(ac*xo+bc*yo+d)/cc
     IF(nn /= 1 .AND. vx3 /= 0.)GO TO 130
     tsz=z2-z1
     tsx=x2-x1
     vt=xo-x1
     IF(tsx /= 0.)GO TO 120
     vt=yo-y1
     tsx=y2-y1
     120 CONTINUE
     zx1=(tsz/tsx)*vt+z1
     GO TO 140
     130 CONTINUE
     zx1=-(vx+vx1*yo+vx2*xo)/vx3
     140 CONTINUE
     IF(ABS(zx-zx1) <= EXP)GO TO 150
     IF(zx1 > zx)GO TO 160
     150 CONTINUE
     m=m+1
     
!     STORES INTERSECTIONS.
     
     xi(m+1)=xo
     yi(m+1)=yo
     zi(m+1)=zx1
     160 jv=jv+5
   END DO
 END DO
 ngx(1)=m
 190 RETURN
END SUBROUTINE hdchk
