SUBROUTINE ptintr (a,aa,b,bb,s,k,eps)
     
!     RETURNS DOUBLE PRECISION VALUES OF X,Y COORDINATES (S) OF
!         POINT OF INTERSECTION (IF ANY) OF LINE SEGMENTS
!         FROM A TO AA AND B TO BB
!           A .NE. AA  AND  B .NE. BB
!     K IS CONDITION FLAG RETURNED --
!         K = 1   LINES INTERSECT AT S
!         K = 0   LINES INTERSECT AT S, AN ENDPOINT OF ONE LINE SEGMENT
!         K =-1   LINES DO NOT INTERSECT
 
 
 DOUBLE PRECISION, INTENT(IN)             :: a(2)
 DOUBLE PRECISION, INTENT(IN)             :: aa(2)
 DOUBLE PRECISION, INTENT(IN)             :: b(2)
 DOUBLE PRECISION, INTENT(IN)             :: bb(2)
 DOUBLE PRECISION, INTENT(OUT)            :: s(2)
 INTEGER, INTENT(OUT)                     :: k
 DOUBLE PRECISION, INTENT(IN OUT)         :: eps(2)
 DOUBLE PRECISION :: p(2)
 DOUBLE PRECISION :: ax,ay,bx,by,aaa,pa,paa,bbb,pb,pbb
 DOUBLE PRECISION :: d
 DOUBLE PRECISION :: dist,x,y,u,v
 
!     EPS ARRAY FOR SIGNIFICANCE TESTING
!         EPS(1) IS AREA, ANGLE LIMIT
!         EPS(2) IS LENGTH LIMIT
 
 
!     DOUBLE PRECISION FUNCTION FOR DISTANCE BETWEEN 2 POINTS
 
 dist(x,y,u,v) = (x-u)**2 +(y-v)**2
 
 x    = 0.d0
 y    = x
 u    = x
 v    = x
 p(1) = 0.d0
 p(2) = 0.d0
 s(1) = 0.d0
 s(2) = 0.d0
 
 k  =-1
 
 ax = aa(1) - a(1)
 ay = aa(2) - a(2)
 bx = bb(1) - b(1)
 by = bb(2) - b(2)
 
 aaa= ax**2 + ay**2
 bbb= bx**2 + by**2
 d  = bx*ay - ax*by
 
!     IS EITHER LINE TOO SHORT?
 
 IF (aaa <= eps(1) .OR. bbb <= eps(1)) RETURN
 
!     ARE A AND B PARALLEL?
 
 IF (DABS(d) > eps(1)) GO TO 80
 
!     A AND B ARE PARALLEL -- ARE THEY SAME LINE?
 
 p(1) = b(1)
 p(2) = b(2)
 IF (dist(b(1),b(2), a(1), a(2)) <= eps(1) .OR.  &
     dist(b(1),b(2),aa(1),aa(2)) <= eps(1)) GO TO 100
 p(1) = bb(1)
 p(2) = bb(2)
 IF (dist(bb(1),bb(2), a(1), a(2)) <= eps(1) .OR.  &
     dist(bb(1),bb(2),aa(1),aa(2)) <= eps(1)) GO TO 100
 
!     A PARALLEL TO B AND NOT SAME LINE
 
 RETURN
 
!     IS A PARALLEL TO Y AXIS?
 
 80 IF (DABS(ax) > eps(2)) GO TO 90
 p(1) = a(1)
 p(2) = b(2) + (p(1)-b(1))*by/bx
 GO TO 100
 90 p(1) = ((b(2)-a(2))*ax*bx + a(1)*ay*bx-b(1)*ax*by)/d
 p(2) =  a(2) + (p(1)-a(1))*ay/ax
 
 100 aaa = aaa + eps(1)
 bbb = bbb + eps(1)
 pa  = dist(p(1),p(2), a(1), a(2))
 pb  = dist(p(1),p(2), b(1), b(2))
 paa = dist(p(1),p(2),aa(1),aa(2))
 pbb = dist(p(1),p(2),bb(1),bb(2))
 
!     POINT OF INTERSECTION NOT ON EITHER SEGMENT
 
 IF (pa > aaa .OR. paa > aaa .OR. pb > bbb .OR. pbb > bbb) RETURN
 
!     LINES INTERSECT AT P
 
 k    = 1
 s(1) = p(1)
 s(2) = p(2)
 
!     LINES INTERSECT AT P, AN ENDPOINT OF ONE SEGMENT
 
 IF ((pa < eps(2) .OR. paa < eps(2)) .OR.  &
     (pb < eps(2) .OR. pbb < eps(2))) k= 0
 RETURN
END SUBROUTINE ptintr
