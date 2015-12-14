SUBROUTINE saxif1 (iopt)
     
!     THIS ROUTINE GENERATES MATRICES WHICH RELATE PRESSURE TO VELOCITY
!     IN A FLUID. IOPT DETERMINES THE ELEMENT TYPE
 
!        IOPT     TYPE
!         0      CAXIF2
!         1      CAXIF3
!         2      CAXIF4
 
 
 INTEGER, INTENT(IN)                      :: iopt
 INTEGER :: nest(100),sil
 DIMENSION       a(9)
 COMMON /sdr2x5/ est(100),nid,sil(4),sv(95)
 COMMON /sdr2x6/ hm(9),r(4),z(4),am(9),coef,en,el,rbar,zbar,  &
     r1n,r2n,rbn1,dr,dz,i1,i2,i3,iret,ij,ik,kj
 EQUIVALENCE     (est(1),nest(1)),(am(1),a(1))
 
 DO  i = 1,44
   sv(i) = 0.0
 END DO
 nid = nest(1)
 IF (iopt-1 < 0) THEN
   GO TO    20
 ELSE IF (iopt-1 == 0) THEN
   GO TO    50
 ELSE
   GO TO    70
 END IF
 
!     CAXIF2 ELEMENTS
 
 20 IF (nest(6) >= 1) GO TO 30
 coef = est(4)*(est(13)-est(9))
 IF (coef == 0.0) RETURN
 sv(3) = 1.0/coef
 sv(4) = -sv(3)
 GO TO 40
 30 IF (nest(6) > 1) GO TO 40
 coef  = est(4)*(est(8)+est(12))
 IF (coef == 0.0) RETURN
 sv(1) = -1.0/coef
 sv(2) =  sv(1)
 40 CONTINUE
 en    = FLOAT(nest(6))
 rbar  = (est(8)+est(12))/2.0
 zbar  = (est(9)+est(13))/2.0
 dr    = est(12) - est(8)
 dz    = est(13) - est(9)
 r1n   = est( 8)**nest(6)
 r2n   = est(12)**nest(6)
 rbn1  = rbar**(nest(6)-1)
 hm(1) = est(13)/(r1n*dz)
 hm(2) =-est( 9)/(r2n*dz)
 hm(3) =-1.0/(r1n*dz)
 hm(4) = 1.0/(r2n*dz)
 el    = SQRT(dz**2 +dr**2)
 coef  = rbn1/(est(4)*el)
 am(1) = en*dr*coef
 am(2) = (en*dr*zbar + rbar*dz)*coef
 am(3) = en*el*coef
 am(4) = en*zbar*el*coef
 sv(5) = am(1)*hm(1) + am(2)*hm(3)
 sv(6) = am(1)*hm(2) + am(2)*hm(4)
 sv(7) = am(3)*hm(1) + am(4)*hm(3)
 sv(8) = am(3)*hm(2) + am(4)*hm(4)
 sil(1)= nest(2)
 sil(2)= nest(3)
 RETURN
 
!     CAXIF3 ELEMENT
 
 50 n    = nest(7)
 en   = FLOAT(n)
 rho  = est(5)
 DO  i = 1,3
   sil(i) = nest(i+1)
   ir   = 4*(i-1) + 9
   r(i) = est(ir  )
   z(i) = est(ir+1)
 END DO
 i1   = 1
 i2   = 2
 i3   = 3
 rbar = (r(i1)+r(i2)+r(i3))/3.0
 zbar = (z(i1)+z(i2)+z(i3))/3.0
 iret = 4
 GO TO 120
 
!     CAXIF4 ELEMENT
 
 70 n    = nest(8)
 en   = FLOAT(n)
 rho  = est(6)*4.0
 DO  i = 1,4
   sil(i) = nest(i+1)
   ir   = 4*(i-1) + 10
   r(i) = est(ir  )
   z(i) = est(ir+1)
 END DO
 rbar = (r(1)+r(2)+r(3)+r(4))/4.0
 zbar = (z(1)+z(2)+z(3)+z(4))/4.0
 i1   = 1
 i2   = 2
 i3   = 3
 iret = 1
 GO TO  120
 90 i3   = 4
 iret = 2
 GO TO  120
 100 i2   = 3
 iret = 3
 GO TO  120
 110 i1   = 2
 iret = 4
 
!     ACTUAL SUBTRIANGLE CALCULATION
 
 120 IF (rho == 0.0) RETURN
 a(1) = 0.0
 a(2) =-1.0/rho
 a(3) = 0.0
 a(5) = a(2)*en
 a(4) = a(5)/rbar
 a(6) = a(4)*zbar
 a(7) = 0.0
 a(8) = 0.0
 a(9) = a(2)
 
 coef = (r(i2)-r(i1))*(z(i3)-z(i1)) - (r(i3)-r(i1))*(z(i2)-z(i1))
 IF (coef == 0.0) RETURN
 hm(1) = (r(i2)*z(i3)-r(i3)*z(i2))/coef
 hm(2) = (r(i3)*z(i1)-r(i1)*z(i3))/coef
 hm(3) = (r(i1)*z(i2)-r(i2)*z(i1))/coef
 hm(4) = (z(i2)-z(i3))            /coef
 hm(5) = (z(i3)-z(i1))            /coef
 hm(6) = (z(i1)-z(i2))            /coef
 hm(7) = (r(i3)-r(i2))            /coef
 hm(8) = (r(i1)-r(i3))            /coef
 hm(9) = (r(i2)-r(i1))            /coef
 DO  j = 1,3
   jcol  = i1
   IF (j == 2) jcol = i2
   IF (j == 3) jcol = i3
   DO  i = 1,3
     ij = (2+iopt)*(i-1) + jcol
     DO  k = 1,3
       ik = 3*(i-1) + k
       kj = 3*(k-1) + j
       sv(ij) = sv(ij) + a(ik)*hm(kj)
     END DO
   END DO
 END DO
 SELECT CASE ( iret )
   CASE (    1)
     GO TO 90
   CASE (    2)
     GO TO 100
   CASE (    3)
     GO TO 110
   CASE (    4)
     GO TO 160
 END SELECT
 
!     THE CENTROID  CALCULATIONS ARE COMPLETE.
 
 160 nsta = 3*(iopt+2)
 ncol = iopt + 2
 IF (iopt == 2) rho = est(6)
 DO  i = 1,ncol
   j = i + 1
   IF (j > ncol) j = j - ncol
   el = SQRT((r(j)-r(i))**2  + (z(j)-z(i))**2)*rho
   
   ik = nsta + 2*ncol*(i-1) + i
   ij = ik + j - i
   sv(ik) = -1.0/el
   sv(ij) = -sv(ik)
   coef   = -en/((r(i)+r(j))*rho)
   ik     = ik + ncol
   ij     = ij + ncol
   sv(ik) = coef
   sv(ij) = coef
 END DO
 RETURN
END SUBROUTINE saxif1
