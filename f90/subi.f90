SUBROUTINE subi (da,dzb,dyb,dar,deta,dzeta,dcgam,dsgam,dee,dxi,tl,  &
        detai,dzetai,dcgami,dsgami,deei,dtlami,dmuy,dmuz,  &
        infl,ioutfl)
     
!     COMPUTES THE IMAGE POINT COORDINATES INSIDE ASSOCIATED BODIES AND
!     THE  MU-Z  MU-Y  ELEMENTS USED IN SUBROUTINE FWMW
 
 eps   = 0.1*dee
 dmuy  = 0.0
 dmuz  = 0.0
 igo   = 1
 psqr  = SQRT(((deta-dyb)*dar)**2 + (dzeta-dzb)**2)
 costh = (deta -dyb)*dar/psqr
 sinth = (dzeta-dzb)/psqr
 ct2   = costh*costh
 st2   = sinth*sinth
 ct3   = costh*ct2
 st3   = sinth*st2
 ycbar = da*(1.0-dar*dar)*ct3 + dyb
 zcbar = da*(dar*dar-1.0)*st3/dar + dzb
 paren = st2 + dar*dar*ct2
 par3  = paren*paren**2
 abar  = da*SQRT(par3)/dar
 abar2 = abar*abar
 IF (infl /= 0)  GO TO  300
 eta1  = deta  - dee*dcgam
 eta2  = deta  + dee*dcgam
 zeta1 = dzeta - dee*dsgam
 zeta2 = dzeta + dee*dsgam
 rho12 = (eta1 - ycbar)**2 + (zeta1-zcbar)**2
 rho22 = (eta2 - ycbar)**2 + (zeta2-zcbar)**2
 etai1 = ycbar + (eta1-ycbar)*abar2/rho12
 etai2 = ycbar + (eta2-ycbar)*abar2/rho22
 zeti1 = zcbar + (zeta1-zcbar)*abar2/rho12
 zeti2 = zcbar + (zeta2-zcbar)*abar2/rho22
 deei  = SQRT((etai2-etai1)**2 + (zeti2-zeti1)**2)/2.0
 detai = (etai1 + etai2)/2.0
 dzetai= (zeti1 + zeti2)/2.0
 dcgami=-(etai2 - etai1)/(2.0*deei)
 dsgami=-(zeti2 - zeti1)/(2.0*deei)
 dxi1  = dxi  - dee*tl
 dxi2  = dxi  + dee*tl
 delxi = dxi1 - dxi2
 dtlami= delxi/(2.0*deei)
 IF (ABS(dar-1.0) <= 0.0001)  GO TO  420
 GO TO 301
 300 CONTINUE
 rho2  = (deta-ycbar)**2 + (dzeta-zcbar)**2
 rho4  = rho2*rho2
 detai = ycbar + (deta -ycbar)*abar2/rho2
 dzetai= zcbar + (dzeta-zcbar)*abar2/rho2
 301 CONTINUE
 SELECT CASE ( igo )
   CASE (    1)
     GO TO 302
   CASE (    2)
     GO TO 303
   CASE (    3)
     GO TO 304
 END SELECT
 302 CONTINUE
 xetai = detai
 xzetai= dzetai
 GO TO 307
 303 CONTINUE
 xetai = etai1
 xzetai= zeti1
 GO TO  307
 304 CONTINUE
 xetai = etai2
 xzetai= zeti2
 307 CONTINUE
 IF (dar < 1.0) GO TO  310
 dybm  = dyb - eps
 dybp  = dyb + eps
 IF (deta >= dyb .AND. xetai < dybm) GO TO  325
 IF (deta <= dyb .AND. xetai > dybp) GO TO  325
 GO TO  320
 310 CONTINUE
 dzbm  = dzb - eps
 dzbp  = dzb + eps
 IF (dzeta >= dzb .AND. xzetai < dzbm) GO TO  325
 IF (dzeta <= dzb .AND. xzetai > dzbp) GO TO  325
 320 CONTINUE
 part1 = ((xetai  - dyb)/da)**2
 part2 = ((xzetai - dzb)/(da*dar))**2
 tedif = part1 + part2 - 1.0
 IF (infl  ==   0) GO TO  400
 IF (tedif <= eps) GO TO  330
 325 CONTINUE
 ioutfl = 0
 GO TO  500
 330 CONTINUE
 ioutfl = 1
 trm1 = (deta-ycbar)**2 - (dzeta-zcbar)**2
 trm2 = 2.0*(deta-ycbar)*(dzeta-zcbar)
 dmuy = -(-dsgam*trm1 + dcgam*trm2)*abar2/rho4
 dmuz = -(-dsgam*trm2 - dcgam*trm1)*abar2/rho4
 GO TO 500
 400 CONTINUE
 IF (tedif > eps) GO TO 325
 IF (igo   ==   3) GO TO 420
 igo  = igo + 1
 GO TO 301
 420 CONTINUE
 ioutfl = 1
 500 CONTINUE
 RETURN
END SUBROUTINE subi
