SUBROUTINE sdhtff
     
!     THIS ROUTINE CALCULATES THE PHASE 1 FLUX-TEMPERATURE RELATIONSHIPS
 
 INTEGER :: sub,nels(18),ip(32),smap(52),strspt,sig
 REAL :: c(12),k,kq(9),dr(3,4),mato,el,zi(3),vec(3),vvec(3)
 COMMON /condas/ consts(5)
 COMMON /sdr2x4/ dumx(109),strspt
 COMMON /sdr2x5/ est(100),ide,sig(32),nq,nsil,NAME(2),k(9),ce(96), dshpb(3,32)
 COMMON /sdr2x6/ sub,imat,af,theta,r(3,32)
 COMMON /hmtout/ mato(6)
 EQUIVALENCE     (consts(1),pi)
 DATA    nels  / 1,1,4,1,4,1,3,5,10,1,1,1,1,4,1,1,1,1 /
 DATA    smap  / 1        ,2        ,3        ,6      ,  &
     1        ,2        ,6        ,5      ,  &
     1        ,4        ,5        ,6      ,  &
     1        ,2        ,3        ,6      ,  &
     1        ,3        ,4        ,8      ,  &
     1        ,3        ,8        ,6      ,  &
     1        ,5        ,6        ,8      ,  &
     3        ,6        ,7        ,8      ,  &
     2        ,3        ,4        ,7      ,  &
     1        ,2        ,4        ,5      ,  &
     2        ,4        ,5        ,7      ,  &
     2        ,5        ,6        ,7      , 4        ,5        ,7        ,8      /
 
 SELECT CASE ( sub )
   CASE (    1)
     GO TO 30
   CASE (    2)
     GO TO 40
   CASE (    3)
     GO TO 40
   CASE (    4)
     GO TO 40
   CASE (    5)
     GO TO 40
   CASE (    6)
     GO TO 50
   CASE (    7)
     GO TO 50
   CASE (    8)
     GO TO 50
   CASE (    9)
     GO TO 50
   CASE (   10)
     GO TO 30
   CASE (   11)
     GO TO 30
   CASE (   12)
     GO TO 30
   CASE (   13)
     GO TO 30
   CASE (   14)
     GO TO 30
   CASE (   15)
     GO TO 30
   CASE (   16)
     GO TO 50
   CASE (   17)
     GO TO 50
   CASE (   18)
     GO TO 30
 END SELECT
 30 k(1) = mato(1)
 nq   = 1
 GO TO  60
 40 k(1) = mato(1)
 k(2) = mato(2)
 k(3) = k(2)
 k(4) = mato(3)
 nq   = 2
 GO TO  60
 50 k(1) = mato(1)
 k(2) = mato(2)
 k(3) = mato(3)
 k(4) = k(2)
 k(5) = mato(4)
 k(6) = mato(5)
 k(7) = k(3)
 k(8) = k(6)
 k(9) = mato(6)
 nq   = 3
 60 CONTINUE
 ip(1)= 1
 ip(2)= 2
 ip(3)= 3
 IF (sub == 17) GO TO 111
 IF (sub /= 3 .AND. sub /= 5) GO TO 100
 
!     MOVE  QUADS TO ELEMENT COORDINATES
!     (CQUAD4? APPEARENTLY UP TO ELEMENT TYPE 52 ONLY)
 
 DO  i = 1,3
   dr(i,1) = r(i,2) - r(i,1)
   dr(i,3) = r(i,3) - r(i,1)
   dr(i,2) = r(i,4) - r(i,2)
 END DO
 CALL saxb (dr(1,3),dr(1,2),dr(1,4))
 
 el   = SQRT(dr(1,1)**2 + dr(2,1)**2 + dr(3,1)**2)
 area = SQRT(dr(1,4)**2 + dr(2,4)**2 + dr(3,4)**2)
 
 DO  i = 1,3
   dr(i,1) = dr(i,1)/el
   dr(i,4) = dr(i,4)/area
 END DO
 
 CALL saxb (dr(1,4),dr(1,1),dr(1,2))
 DO  i = 1,3
   dr(i,4) = r(i,4) - r(i,1)
 END DO
 CALL gmmats (dr(1,1),2,3,0,dr(1,3),2,3,1,kq)
 dr(1,3) = kq(1)
 dr(1,4) = kq(2)
 dr(2,3) = kq(3)
 dr(2,4) = kq(4)
 dr(1,2) = el
 dr(1,1) = 0.0
 dr(2,1) = 0.0
 dr(2,2) = 0.0
 GO TO 120
 100 IF (sub /= 2 .AND. sub /= 4) GO TO 120
 
!     MOVE  TRIANGLES TO ELEMENT COORDINATES
!     (CTRIA3?)
 
 DO  i = 1,3
   dr(i,1) = r(i,2) - r(i,1)
   dr(i,2) = r(i,3) - r(i,1)
 END DO
 
 el   = dr(1,1)**2 + dr(2,1)**2 + dr(3,1)**2
 el   = SQRT(el)
 area = sadotb(dr(1,1),dr(1,2))/el
 CALL saxb (dr(1,1),dr(1,2),dr(1,3))
 dr(2,3) = SQRT(dr(1,3)**2 + dr(2,3)**2 + dr(3,3)**2)/el
 dr(1,3) = area
 dr(1,1) = 0.0
 dr(1,2) = el
 dr(2,1) = 0.0
 dr(2,2) = 0.0
 GO TO 120
 
!     IS2D8-CENTROID ONLY-WE NEED TO CONVERT ONLY GRIDS 5-8 TO LOCAL
!     COORDS
 
 111 DO  i = 1,3
   zi(i) = r(i,2) - r(i,1)
 END DO
 zlen  = SQRT(zi(1)**2 + zi(2)**2 + zi(3)**2)
 DO  i = 1,3
   zi(i) = zi(i)/zlen
 END DO
 DO  i = 5,8
   DO  j = 1,3
     vec(j) = r(j,i) - r(j,1)
   END DO
   dr(1,i-4) = sadotb(vec,zi)
   CALL saxb (zi,vec,vvec)
   dr(2,i-4) = SQRT(vvec(1)**2 + vvec(2)**2 + vvec(3)**2)
 END DO
 120 CONTINUE
 
!     LOOP  ON  SUBELEMENTS  (ONE FOR MOST)
 
 fact = 0.0
 nel  = nels(sub)
 xels = FLOAT(nel)
 DO  iel = 1,nel
   
   SELECT CASE ( sub )
     CASE (    1)
       GO TO 130
     CASE (    2)
       GO TO 160
     CASE (    3)
       GO TO 160
     CASE (    4)
       GO TO 140
     CASE (    5)
       GO TO 140
     CASE (    6)
       GO TO 200
     CASE (    7)
       GO TO 220
     CASE (    8)
       GO TO 240
     CASE (    9)
       GO TO 240
     CASE (   10)
       GO TO 330
     CASE (   11)
       GO TO 330
     CASE (   12)
       GO TO 330
     CASE (   13)
       GO TO 330
     CASE (   14)
       GO TO 330
     CASE (   15)
       GO TO 330
     CASE (   16)
       GO TO 285
     CASE (   17)
       GO TO 291
     CASE (   18)
       GO TO 330
   END SELECT
   
!     RODS,BARS, ETC.
   
   130 el = 0.0
   DO  i = 1, 3
     el = el + (r(i,1)-r(i,2))**2
   END DO
   el = SQRT(el)
   c(1) = -1.0/el
   c(2) =  1.0/el
   np = 2
   GO TO 300
   
!     RING ELEMENTS, TRIANGLES AND QUADRILATERALS
   
   140 af = 1.0
   160 DO  i = 1,3
     ig = i + iel - 1
     IF (ig > 4) ig = ig - 4
     ip(i) = ig
   END DO
   i1 = ip(1)
   i2 = ip(2)
   i3 = ip(3)
   area = dr(1,i1)*(dr(2,i2)-dr(2,i3)) + dr(1,i2)*(dr(2,i3)-dr(2,i1))  &
       + dr(1,i3)*(dr(2,i1)-dr(2,i2))
   c(1) = (dr(2,i2) - dr(2,i3))/area
   c(2) = (dr(2,i3) - dr(2,i1))/area
   c(3) = (dr(2,i1) - dr(2,i2))/area
   c(4) = (dr(1,i3) - dr(1,i2))/area
   c(5) = (dr(1,i1) - dr(1,i3))/area
   c(6) = (dr(1,i2) - dr(1,i1))/area
   
   np   = 3
   GO TO 300
   
!     SOLID ELEMENTS
   
   200 DO  i = 1,4
     ip(i) = i
   END DO
   GO TO 260
   
!     WEDGE
   
   220 lrow = 4*iel - 4
   DO  i = 1,4
     i1 = lrow + i
     ip(i) = smap(i1)
   END DO
   GO TO 260
   
!     HEXA1 AND HEXA2 ELEMENTS
   
   240 lrow = 4*iel + 8
   DO  i = 1,4
     i1 = lrow + i
     ip(i) = smap(i1)
   END DO
   260 i1 = ip(1)
   DO  i = 1,3
     ig = ip(i+1)
     DO  j = 1,3
       dr(j,i) = r(j,ig) - r(j,i1)
     END DO
   END DO
   
!     NO NEED TO COMPUTE DETERMINANT SINCE IT IS NOT USED SUBSEQUENTLY.
   
   ising = -1
   CALL invers (3,dr,3,c,0,determ,ising,c(4))
   DO  i = 1,3
     ig = 4*i - 4
     c(ig+1) = 0.0
     DO  j = 2,4
       i1 = ig + j
       c(i1  ) = dr(j-1,i)
       c(ig+1) = c(ig+1) - c(i1)
     END DO
   END DO
   np = 4
   GO TO 300
   
!     ISOPARAMETRIC SOLIDS
   
   285 ig = 0
   DO  i = 1,3
     DO  j = 1,nsil
       ig = ig + 1
       ce(ig) = dshpb(i,j)
     END DO
   END DO
   CYCLE
   
!     IS2D8- SINCE CENTROID ONLY, WE CAN EASILY COMPUTE SHAPE FUNCTIONS
!     DERIVATIVES, JACOBIAN,ETC.. THE FINAL RESULT OF DNDX,DNDY=DNL IS
!     GIVEN
   
   291 x68 = dr(1,2) - dr(1,4)
   x57 = dr(1,1) - dr(1,3)
   y68 = dr(2,2) - dr(2,4)
   y57 = dr(2,1) - dr(2,3)
   denom = -x68*y57 + x57*y68
   DO  i = 1,24
     ce(i)  = 0.
   END DO
   ce( 5) = y68/denom
   ce( 6) =-y57/denom
   ce( 7) =-y68/denom
   ce( 8) = y57/denom
   ce(13) =-x68/denom
   ce(14) = x57/denom
   ce(15) = x68/denom
   ce(16) =-x57/denom
   CYCLE
   
!     SUPERIMPOSE C MATRICES ONTO CE MATRICES OF THE WHOLE ELEMENT
   
   300 DO  i = 1,np
     DO  j = 1,nq
       i1 = np*(j-1) + i
       ig = nsil*(j-1) + ip(i)
       ce(ig) = ce(ig) + c(i1)/xels
     END DO
   END DO
   CYCLE
   
!     BOUNDARY HEAT CONVECTION ELEMENTS
   
   330 itype = sub - 9
   IF (itype > 7) RETURN
   SELECT CASE ( itype )
     CASE (    1)
       GO TO 340
     CASE (    2)
       GO TO 350
     CASE (    3)
       GO TO 370
     CASE (    4)
       GO TO 380
     CASE (    5)
       GO TO 380
     CASE (    6)
       GO TO 350
     CASE (    7)
       GO TO 350
   END SELECT
   340 np = 1
   c(1) = 1.0
   fact = af*k(1)
   GO TO 410
   350 np = 2
   c(1) = 0.5
   c(2) = 0.5
   el = SQRT((r(1,1)-r(1,2))**2 + (r(2,1)-r(2,2))**2 + (r(3,1)-r(3,2))**2)
   fact = af*el*k(1)
   GO TO 410
   
!     RING SURFACE
   
   370 el   = ((r(1,2)-r(1,1))**2 + (r(3,2)-r(3,1))**2)
   fact = 3.0*(r(1,1) + r(1,2))
   c(1) = (2.0*r(1,1) + r(1,2))/fact
   c(2) = (r(1,1) + 2.0*r(1,2))/fact
   fact = (r(1,1) + r(1,2))*pi*SQRT(el)*k(1)
   np   = 2
   GO TO 410
   
!     TRIANGLES  (ALSO FOR SUBELEMENT OF QUAD)
   
   380 DO  i = 1,3
     ig = i + iel - 1
     IF (ig > 4) ig = ig - 4
     ip(i) = ig
   END DO
   i1 = ip(1)
   i2 = ip(2)
   i3 = ip(3)
   DO  i = 1,3
     dr(i,1) = r(i,i2) - r(i,i1)
     dr(i,2) = r(i,i3) - r(i,i1)
   END DO
   CALL saxb (dr(1,1),dr(1,2),dr(1,3))
   area = (SQRT(dr(1,3)**2 + dr(2,3)**2 +dr(3,3)**2))/2.0
   IF (itype == 5) area = area/2.0
   fact = fact + area*mato(1)
   c(1) = 1.0/3.0
   c(2) = c(1)
   c(3) = c(1)
   np   = 3
   
!     SUPERIMPOSE C MATRIX INTO CE MATRIX
   
   410 DO  i = 1,np
     ig = ip(i)
     ce(ig) = ce(ig) + c(i)/xels
     ig = ip(i) + 4
     ce(ig) = ce(ig) - c(i)/xels
   END DO
   k(1) = fact
 END DO
 RETURN
END SUBROUTINE sdhtff
