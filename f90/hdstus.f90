SUBROUTINE hdstus (oj,tmj,xxx,tgm,rv,rvi,tgi,zm,nno,ii,h,im,jxt,  &
        zj,nc,zmi,ccc,lz)
     
!     THIS SUBROUTINE DETERMINES THE VISIBILITY OF AN ARBITRARY POINT
!     BY DRAWING A LINE FROM THE POINT IN QUESTION TO INFINITY AND
!     COUNTING THE NUMBER OF TIMES IT CROSSES THE BOUNDARIES OF A
!     RELEVANT ELEMENT.
 
 
 REAL, INTENT(IN)                         :: oj
 REAL, INTENT(IN)                         :: tmj
 REAL, INTENT(IN)                         :: xxx(1)
 REAL, INTENT(IN)                         :: tgm(1)
 REAL, INTENT(IN)                         :: rv(1)
 REAL, INTENT(IN)                         :: rvi(1)
 REAL, INTENT(IN)                         :: tgi(1)
 REAL, INTENT(IN)                         :: zm(1)
 INTEGER, INTENT(IN)                      :: nno(1)
 INTEGER, INTENT(IN)                      :: ii
 REAL, INTENT(IN OUT)                     :: h(8)
 INTEGER, INTENT(OUT)                     :: im
 INTEGER, INTENT(IN OUT)                  :: jxt
 REAL, INTENT(IN OUT)                     :: zj
 INTEGER, INTENT(IN OUT)                  :: nc
 REAL, INTENT(IN OUT)                     :: zmi(1)
 REAL, INTENT(IN)                         :: ccc(1)
 INTEGER, INTENT(IN)                      :: lz
 
 COMMON /go3/ l0,l1,l00,l01,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13
 
 ggk = .015
 ei  = 0
 im  = 0
 10 CONTINUE
 IF (ei >= 1.) GO TO 70
 ei  = ei + .2
 d   = ei*oj - tmj
 loop60:  DO  jo = 1,ii
   ggk = .015
   i   = 0
   jg  = nno(l4+jo)
   js  = l13 + (jg-1)*lz
   jt  = l12 + (jg-1)*5
   
!     PRELIMINARY CHECK TO SEE IF THE POINT IS OUTSIDE THE BOUNDARY
!     BOXES IN THE X,Y,Z DIMENSIONS.
   
   IF (tmj >= rv(l7+jg) .OR. tmj <= rvi(l8+jg)) CYCLE loop60
   IF (oj >= tgi(l6+jg) .OR. oj <= tgm(l5+jg) ) CYCLE loop60
   IF (zj >= zm(l2+jg)) CYCLE loop60
   vx  = xxx(4+jt)
   vx1 = xxx(2+jt)*tmj
   vx2 = xxx(1+jt)*oj
   zs  =-(vx+vx1+vx2)/xxx(3+jt)
   IF (ABS(zj-zs) < ggk) CYCLE loop60
   IF (zj >= zs) CYCLE loop60
   ns  = xxx(5+jt)
   ib  = ns*5
   IF (h(8) == 1.) GO TO 25
   DO  j = 1,ib,5
     ggk = .015
     nsub= j + 1 + js
     IF (ABS(ccc(nsub)) >= 100.) ggk = ALOG10(ABS(ccc(nsub)))
     ve  = oj
     IF (ccc(j+js) == 0.) ve = tmj
     s   = ve - ccc(j+3+js)
     s1  = ve - ccc(j+4+js)
     yg  = tmj
     IF (ccc(j+js) /= 0.) GO TO 15
     dy  =-ccc(j+2+js)/ccc(j+1+js)
     yg  = oj
     GO TO 16
     15 CONTINUE
     dy  =-ccc(j+2+js) - ccc(j+1+js)*oj
     16 CONTINUE
     IF (ABS(yg-dy) < ggk .AND. s*s1 <= 0.) CYCLE loop60
   END DO
   25 CONTINUE
   
!     THE FOLLOWING CODE COUNTS THE INTERSECTIONS OF BOUNDARIES
!     OF A GIVEN ELEMENT WITH THE INFINITE LINE AND CHECKS,IF INSIDE
!     OF THE BOUNDARY, WHETHER OR NOT THE POINT IS BEHIND OR IN FRONT
!     OF THE ELEMENT.
   
   DO  j = 1,ib,5
     t  =-ccc(j+2+js)  + ccc(j+js)*d
     r  = ei*ccc(j+js) + ccc(j+1+js)
     IF (r == 0.) CYCLE
     t  = t/r
     IF (t < oj) CYCLE
     IF (ccc(j+js) /= 0.) GO TO 30
     t  = ei*t - d
     30 CONTINUE
     s  = t - ccc(j+3+js)
     s1 = t - ccc(j+4+js)
     IF (s == 0. .OR. s1 == 0.) GO TO 10
     IF (s*s1 >= 0.) CYCLE
     i  = i + 1
   END DO
   IF (MOD(i,2) == 0) CYCLE loop60
   im = 1
   GO TO 70
 END DO loop60
 im = 0
 70 CONTINUE
 RETURN
END SUBROUTINE hdstus
