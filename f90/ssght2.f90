SUBROUTINE ssght2 (FILE,delta,uni)
     
!     THIS ROUTINE USES THE TEMPERATURE VECTOR DATA TO CALCULATE
!     LOAD VECTER TERMS WITH THE EQUATION-
 
!     DELTAP = (K(TI) - K(TO))*T1
!          WHERE       TO IS THE INITIAL TEMPERATURE
!                      TI IS THE  NEW  TEMPERATURE VECTOR
!                      K  IS THE TEMPERATURE DEPENDENT CONDUCTIVITY
!                          MATRIX
!                      DELTAP  IS THE NONLINEAR LOAD
 
 
 INTEGER, INTENT(IN OUT)                  :: FILE
 REAL, INTENT(OUT)                        :: delta(1)
 REAL, INTENT(IN)                         :: uni(1)
 INTEGER :: elid,sub,sil,npts(15),nels(15),ip(4),smap(52),  &
      flag,sindx(4),subr(2)
 REAL :: mato,matout,temp
 DOUBLE PRECISION :: constd,drtemp(3,4),drt(4,4),c(12),k(9),kq(9),  &
     dr(3,4),t1(8),el,area,rbar,pi,fact,determ,dadotb
 COMMON/ condad/  constd(5)
 COMMON/ matin /  matid,inflag,temp,dum,sinth,costh
 COMMON/ hmtout/  matout(6)
 COMMON/ estout/  elid,sub,NAME(2),sil(8),imat,af,theta,r(3,8), mato(6)
 EQUIVALENCE      (constd(1),pi)
 DATA    npts  /  2,3,4,3,4,4,6,8, 8,1,2,2,3,4,2 /
 DATA    nels  /  1,1,4,1,4,1,3,5,10,1,1,1,1,4,1 /
 DATA    subr  /  4HSSGH ,4HT2   /
 DATA    smap  /  1      ,2      ,3      ,6      ,  &
     1      ,2      ,6      ,5      , 1      ,4      ,5      ,6      ,  &
     1      ,2      ,3      ,6      , 1      ,3      ,4      ,8      ,  &
     1      ,3      ,8      ,6      , 1      ,5      ,6      ,8      ,  &
     3      ,6      ,7      ,8      , 2      ,3      ,4      ,7      ,  &
     1      ,2      ,4      ,5      , 2      ,4      ,5      ,7      ,  &
     2      ,5      ,6      ,7      , 4      ,5      ,7      ,8      /
 
!     READ DATA, 45 WORDS PER ELEMENT.
 
 10 CALL READ (*480,*470,FILE,elid,45,0,flag)
 
!     CALCULATE AVERAGE ELEMENT TEMPERATURE
 
 np    = npts(sub)
 xpts  = FLOAT(np)
 IF (sub > 9) xpts = xpts*2.0
 
 temp  = 0.0
 DO  i = 1,np
   ltemp = sil(i)
   temp  = temp + uni(ltemp)
   IF (sub <= 9) CYCLE
   IF (sil(i+4) == 0) CYCLE
   ltemp = sil(i+4)
   temp  = temp + uni(ltemp)
 END DO
 temp  = temp/xpts
 
!     SET UP CALL TO MATERIAL SUBROUTINE
 
 inflag = 1
 IF (sub >= 2 .AND. sub <= 5) inflag = 2
 IF (sub >= 6 .AND. sub <= 9) inflag = 3
 sinth = 0.0
 costh = 1.0
 IF (theta == 0.0 .OR. inflag /= 2) GO TO 30
 sinth = SIN(theta)
 costh = COS(theta)
 30 matid = imat
 CALL hmat (elid)
 
!     SUBTRACT  CONDUCTIVITY AT INITIAL TEMPERATURE AND PLACE IN  MATRIX
 
 IF (inflag == 2) GO TO 40
 IF (inflag == 3) GO TO 50
 k(1) = matout(1) - mato(1)
 GO TO 60
 40 k(1) = matout(1) - mato(1)
 k(2) = matout(2) - mato(2)
 k(3) = k(2)
 k(4) = matout(3) - mato(3)
 GO TO 60
 50 k(1) = matout(1) - mato(1)
 k(2) = matout(2) - mato(2)
 k(3) = matout(3) - mato(3)
 k(4) = k(2)
 k(5) = matout(4) - mato(4)
 k(6) = matout(5) - mato(5)
 k(7) = k(3)
 k(8) = k(6)
 k(9) = matout(6) - mato(6)
 60 CONTINUE
 ip(1) = 1
 ip(2) = 2
 ip(3) = 3
 IF (sub /= 3 .AND. sub /= 5) GO TO 100
 
!     MOVE  QUADS TO ELEMENT COORDINATES
 
 DO  j = 1,2
   l  = 1
   m  = 2
   i1 = 1
   i2 = 2
   i3 = 3
   i4 = 4
   IF (j == 1) GO TO 65
   l  = 3
   m  = 4
   i1 = 3
   i2 = 4
   i3 = 1
   i4 = 2
   65 CONTINUE
   DO  i = 1,3
     dr(i,1) = r(i,i2) - r(i,i1)
     dr(i,3) = r(i,i3) - r(i,i1)
     dr(i,2) = r(i,i4) - r(i,i2)
   END DO
   CALL daxb (dr(1,3),dr(1,2),dr(1,4))
   
   area = DSQRT(dr(1,4)**2 + dr(2,4)**2 + dr(3,4)**2)
   
   DO  i = 1,3
     dr(i,4) = dr(i,4)/area
   END DO
   el = dr(1,4)*dr(1,1) + dr(2,4)*dr(2,1) + dr(3,4)*dr(3,1)
   DO  i = 1,3
     dr(i,1) = dr(i,1) - el*dr(i,4)
   END DO
   el = DSQRT(dr(1,1)**2 + dr(2,1)**2 + dr(3,1)**2)
   DO  i = 1,3
     dr(i,1) = dr(i,1)/el
   END DO
   
   CALL daxb (dr(1,4),dr(1,1),dr(1,2))
   DO  i = 1,3
     dr(i,4) = r(i,i4) - r(i,i1)
   END DO
   CALL gmmatd (dr(1,1),2,3,0, dr(1,3),2,3,1, kq)
   drt(l,3) = kq(1)
   drt(l,4) = kq(2)
   drt(m,3) = kq(3)
   drt(m,4) = kq(4)
   drt(l,1) = 0.0D0
   drt(l,2) = el
   drt(m,1) = 0.0D0
   drt(m,2) = 0.0D0
 END DO
 GO TO 120
 100 IF (sub /= 2 .AND. sub /= 4) GO TO 120
 
!     MOVE  TRIANGLES TO ELEMENT COORDINATES
 
 DO  i = 1,3
   dr(i,1) = r(i,2) - r(i,1)
   dr(i,2) = r(i,3) - r(i,1)
 END DO
 
 el   = dr(1,1)**2 + dr(2,1)**2 + dr(3,1)**2
 el   = DSQRT(el)
 area = dadotb( dr(1,1),dr(1,2))/el
 CALL daxb (dr(1,1), dr(1,2), dr(1,3))
 dr(2,3) = DSQRT(dr(1,3)**2 + dr(2,3)**2 + dr(3,3)**2)/el
 dr(1,3) = area
 dr(1,1) = 0.0D0
 dr(1,2) = el
 dr(2,1) = 0.0D0
 dr(2,2) = 0.0D0
 120 CONTINUE
 
!     LOOP  ON  SUBELEMENTS  (ONE FOR MOST)
 
 nel = nels(sub)
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
       GO TO  200
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
   END SELECT
   
!     RODS,BARS, ETC.
   
   130 c(1) = 1.0D0
   c(2) =-1.0D0
   el = (r(1,2)-r(1,1))**2 + (r(2,2)-r(2,1))**2 + (r(3,2)-r(3,1))**2
   el = DSQRT(el)
   kq(1) = af*k(1)/el
   ip(1) = 1
   ip(2) = 2
   np = 2
   nq = 1
   GO TO 300
   
!     RING ELEMENTS, TRIANGLES AND QUADRILATERALS
   
   140 rbar = 0.0
   DO  i = 1,3
     ig = i + iel - 1
     IF (ig > 4) ig = ig - 4
     rbar  = rbar + r(1,ig)
     ip(i) = ig
   END DO
   af = rbar/3.0*pi
   IF (sub == 5) GO TO 160
   i1 = ip(1)
   i2 = ip(2)
   i3 = ip(3)
   GO TO 180
   160 j  = 1
   i1 = 1
   i2 = 2
   i3 = 3
   IF (iel == 2 .OR. iel == 4) i3 = 4
   ip(1) = 1
   ip(2) = 2
   ip(3) = 3
   IF (iel == 1) GO TO 165
   ip(3) = 4
   IF (iel == 2) GO TO 165
   j     = 3
   ip(1) = 3
   ip(2) = 4
   ip(3) = 1
   IF (iel == 3) GO TO 165
   ip(3) = 2
   165 DO  i = 1,4
     dr(1,i) = drt(j,i)
     dr(2,i) = drt(j+1,i)
   END DO
   180 CONTINUE
   area = dr(1,i1)*(dr(2,i2) - dr(2,i3)) + dr(1,i2)*(dr(2,i3) - dr(2,i1))  &
       + dr(1,i3)*(dr(2,i1) - dr(2,i2))
   
   c(1) = (dr(2,i2) - dr(2,i3))/area
   c(2) = (dr(2,i3) - dr(2,i1))/area
   c(3) = (dr(2,i1) - dr(2,i2))/area
   
   c(4) = (dr(1,i3) - dr(1,i2))/area
   c(5) = (dr(1,i1) - dr(1,i3))/area
   c(6) = (dr(1,i2) - dr(1,i1))/area
   
   IF (sub == 3 .OR. sub == 5) area = area/2.0D0
   DO  i = 1,4
     kq(i) = k(i)*area*af/2.0D0
   END DO
   
   np = 3
   nq = 2
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
     i1 = lrow +i
     ip(i) = smap(i1)
   END DO
   i1 = ip(1)
   260 DO  i = 1,3
     ig = ip(i+1)
     DO  j = 1,3
       dr(j,i) = r(j,ig) - r(j,i1)
     END DO
   END DO
   
!     COMPUTE INVERSE AND BRING ALONG THE DETERMINANT FROM INVERD.
   
   ising = 0
   CALL inverd (3, dr(1,1),3,c(1), 0, determ,ising,c(5))
   DO  i = 1,3
     ig = 4*i - 4
     c(ig+1) = 0.0D0
     DO  j = 2,4
       i1 = ig + j
       c(i1  ) = dr(j-1,i)
       c(ig+1) = c(ig+1) - c(i1)
     END DO
   END DO
   fact = determ/6.0D0
   IF (sub == 9)  fact = fact/2.0D0
   DO  i = 1,9
     kq(i) = k(i)*fact
   END DO
   np = 4
   nq = 3
   
!     PERFORM  MATRIX MULTPLIES FOR EACH SUBELEMENT
!                               T
!                       DP  =  C  K  C * T1
   
   300 DO  i = 1,np
     ig = ip(i)
     ltemp = sil(ig)
     t1(i) = uni(ltemp)
     sindx(i) = sil(ig)
   END DO
   CALL gmmatd (c,nq,np,0, t1,np,1,0, drtemp)
   CALL gmmatd (kq,nq,nq,0, drtemp,nq,1,0, drtemp(1,3))
   CALL gmmatd (c,nq,np,1, drtemp(1,3),nq,1,0, kq(1))
   DO  i = 1,np
     ig = sindx(i)
     delta(ig)  = delta(ig) + kq(i)
   END DO
   CYCLE
   
!     BOUNDARY HEAT CONVECTION ELEMENTS
   
   330 itype = sub - 9
   IF (itype > 7) GO TO 10
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
   kq(1) = af*k(1)
   GO TO 410
   350 np = 2
   el = (r(1,2)-r(1,1))**2 + (r(2,2)-r(2,1))**2 + (r(3,2)-r(3,1))**2
   kq(1) = af*k(1)*DSQRT(el)/3.0D0
   kq(2) = kq(1)/2.0D0
   kq(3) = kq(2)
   kq(4) = kq(1)
   GO TO 410
   
!     RING SURFACE
   
   370 el    = ((r(1,2)-r(1,1))**2 + (r(3,2)-r(3,1))**2)
   c(1)  = pi*k(1)*DSQRT(el)/6.0D0
   kq(1) = c(1)*(3.0D0*r(1,1) + r(1,2))
   kq(2) = c(1)*(      r(1,1) + r(1,2))
   kq(3) = kq(2)
   kq(4) = c(1)*(      r(1,1) + 3.0D0*r(1,2))
   np= 2
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
   CALL daxb (dr(i,1),dr(i,2),dr(i,3))
   area = DSQRT(dr(1,3)**2 + dr(2,3)**2 + dr(3,3)**2)/12.0D0
   IF (itype == 5) area = area/2.0D0
   kq(1) = area*k(1)
   kq(2) = kq(1)/2.0D0
   kq(3) = kq(2)
   kq(4) = kq(2)
   kq(5) = kq(1)
   kq(6) = kq(2)
   kq(7) = kq(2)
   kq(8) = kq(2)
   kq(9) = kq(1)
   np    = 3
   
!     PERFORM MATRIX MULTIPLY, FIRST GET TEMPERATURE VECTOR
   
   410 DO  i = 1,np
     ig = ip(i)
     ltemp  = sil(ig)
     t1(i)  = uni(ltemp)
     IF (sil(ig+4) /= 0) GO TO 420
     t1(i+4) = 0.0D0
     CYCLE
     420 ltemp   = sil(ig+4)
     t1(i+4) = uni(ltemp)
   END DO
   CALL gmmatd (kq(1),np,np,0, t1(1),np,1,0, c)
   CALL gmmatd (kq(1),np,np,0, t1(5),np,1,0, c(5))
   DO  i = 1,np
     ig  = ip(i)
     ipg = sil(ig)
     delta(ipg) = delta(ipg) + c(i) - c(i+4)
     ig = ig + 4
     IF (sil(ig) > 0.0) THEN
       GO TO   440
     ELSE
       GO TO   450
     END IF
     440 ipg = sil(ig)
     delta(ipg) = delta(ipg) + c(i+4) - c(i)
     450 CONTINUE
   END DO
   
 END DO
 GO TO 10
 470 RETURN
 
 480 j = -2
 CALL mesage (FILE,j,subr)
 RETURN
END SUBROUTINE ssght2
