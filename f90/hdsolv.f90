SUBROUTINE hdsolv (ixr,j,xxx,ccc,ii,nno,nit,x21,y21,z21,iia,nc,  &
        zm,zmi,lz)
     
!     THIS SUBROUTINE SOLVES FOR THE LINES OF INTERSECTION RESULTING
!     FROM THE INTERSECTIONS OF THE JTH ELEMENT WITH THE OTHER
!     RELEVANT ELEMENTS.
 
 
 
 
 INTEGER, INTENT(IN OUT)                  :: ixr
 INTEGER, INTENT(IN)                      :: j
 REAL, INTENT(IN OUT)                     :: xxx(1)
 REAL, INTENT(IN)                         :: ccc(1)
 INTEGER, INTENT(IN)                      :: ii
 INTEGER, INTENT(IN)                      :: nno(1)
 INTEGER, INTENT(IN)                      :: nit
 REAL, INTENT(IN OUT)                     :: x21(1)
 REAL, INTENT(IN OUT)                     :: y21(1)
 REAL, INTENT(IN OUT)                     :: z21(1)
 INTEGER, INTENT(IN OUT)                  :: iia(1)
 INTEGER, INTENT(IN OUT)                  :: nc
 REAL, INTENT(IN OUT)                     :: zm(1)
 REAL, INTENT(IN)                         :: zmi(1)
 INTEGER, INTENT(IN)                      :: lz
 INTEGER :: xcc,xasolv,yasolv,zasolv
 DIMENSION  iv(2)
 COMMON /hdptrs/ xdum,xcc,xasolv,yasolv,zasolv
 COMMON /zzzzzz/ rz(1)
 COMMON /go3   / l0,l1,l00,l01,l2,l3,l4,l5,l6,l7,l8,l9,l10, l11,l12,l13
 
 ers = .015
 er  =  ers
 exx = .015
 EXP = .015
 jt  = l12 + (j-1)*5
 jb  = l13 + (j-1)*lz
 IF (ii == 0) GO TO 80
 a3 = xxx(1+jt)
 b3 = xxx(2+jt)
 c3 = xxx(3+jt)
 d3 = xxx(4+jt)
 IF (xxx(jt+3) == 0.) GO TO 80
 loop70:  DO  l = 1,ii
   k  = nno(l4+l)
   
!     CHECKS TO SEE IF THIS RELEVANT ELEMENT IS TO BE CONSIDERED FOR
!     INTERSECTION
   
   IF (k > 0) GO TO 9
   CYCLE loop70
   9 CONTINUE
   IF (k < j) CYCLE loop70
   jx = l12 + (k-1)*5
   IF (zm(l2+j) < zmi(l3+k)) CYCLE loop70
   IF (ABS(xxx(3+jx)) < ers) CYCLE loop70
   mt = 0
   a4 = xxx(1+jx)
   b4 = xxx(2+jx)
   c4 = xxx(3+jx)
   d4 = xxx(4+jx)
   
!     DETERMINES THE EQUATION OF LINE OF INTERSECTION.
   
   b = a3*c4 - a4*c3
   a = b3*c4 - b4*c3
   c = d3*c4 - d4*c3
   IF (a == 0. .AND. b == 0.) CYCLE loop70
   IF (a /= 0.) GO TO 10
   a = 0
   c = c/b
   b = 1
   GO TO 20
   10 CONTINUE
   b = b/a
   c = c/a
   a = 1
   20 CONTINUE
   iv(1) = j
   iv(2) = k
   DO  m = 1,2
     jv = 1
     i  = iv(m)
     jj = l13 + (i-1)*lz
     ig = l12 + (i-1)*5 + 5
     nk = xxx(ig)
     DO  ix = 1,nk
       a1 = ccc(jv+  jj)
       b1 = ccc(jv+1+jj)
       c1 = ccc(jv+2+jj)
       
!     CHECK TO BE SURE LINE OF INTERSECTION IS NOT BOUNDARY LINE
!     OF THE JTH SET.
       
       s  = a1 + b1 + c1
       s1 = a  + b  + c
       e  = ABS(s-s1)
       s  = a1*50 + b1*50 + c1
       s1 = a *50 + b *50 + c
       f  = ABS(s-s1)
       IF (f < EXP .AND. e < EXP) CYCLE loop70
       
       
!     DETERMINES THE POINTS OF INTERSECTIONS OF THE LINE OF INTERSECTION
!     WITH OTHER LINES OF RELEVANT ELEMENTS.
       
       
       t = a1*b - b1*a
       IF (ABS(t) < er) GO TO 50
       xo = (c1*a-c*a1)/t
       IF (a /= 0.) GO TO 30
       yo = -c1 - b1*xo
       GO TO 40
       30 CONTINUE
       yo = -c - b*xo
       40 CONTINUE
       t = xo
       IF (a1 == 0.) t = yo
       s  = t - ccc(jv+4+jj)
       s1 = t - ccc(jv+3+jj)
       IF (s*s1 > 0.) GO TO 50
       mt = mt + 1
       
!     STORE THE PTS OF INTERSECTIONS.
       
       rz(xasolv-1+mt) = xo
       rz(yasolv-1+mt) = yo
       rz(zasolv-1+mt) =-(d3+a3*xo+b3*yo)/c3
       zt = -(d4+a4*xo+b4*yo)/c4
       IF (ABS(zt-rz(zasolv-1+mt)) > exx) CYCLE loop70
       50 jv = jv + 5
     END DO
   END DO
   CALL hdstat (mt,nit,ixr,x21,y21,z21,iia,iv,a,b,c,j,  &
       rz(xasolv),rz(yasolv),rz(zasolv),ccc,xxx,lz)
 END DO loop70
 80 CONTINUE
 nr = 5*xxx(5+jt)
 DO  is = 1,nr
   rz(xcc-1+is) = ccc(is+jb)
 END DO
 xxx(5+jt) = xxx(5+jt) + nit
 RETURN
END SUBROUTINE hdsolv
