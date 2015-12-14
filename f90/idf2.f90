SUBROUTINE idf2(ee,e2,eta01,zet01,a2r,a2i,b2r,b2i,c2r,c2i,  &
                r1sqx,diijr,diiji)

!   ***   INTEGRATES THE NONPLANAR PARTS OF THE INCREMENTAL
!         OSCILLATORY KERNELS FOR UNSTEADY CASES

 eps  = 0.0001
 azet = ABS(zet01)
 deno = r1sqx-e2
 parn = eta01**2 + zet01**2
 facr = parn*a2r + eta01*b2r + c2r
 faci = parn*a2i + eta01*b2i + c2i
 eta02=eta01**2
 zet02= zet01**2
 IF  ((azet/ee)  <=  0.001)  GO TO  120
 test0= ABS((r1sqx-e2)/(2.0*ee*azet))
 IF (test0 > 0.1)  GO TO 120
 den2 = (eta01+ee)**2+zet02
 den3 = (eta01-ee)**2+zet02
 fac2a= r1sqx*eta01+(eta02-zet02)*ee
 fac3a= r1sqx*eta01-(eta02-zet02)*ee
 fac2b= r1sqx+eta01*ee
 fac3b= r1sqx-eta01*ee
 trm2r= (fac2a*a2r+fac2b*b2r+(eta01+ee)*c2r)/den2
 trm2i= (fac2a*a2i+fac2b*b2i+(eta01+ee)*c2i)/den2
 trm3r=-(fac3a*a2r+fac3b*b2r+(eta01-ee)*c2r)/den3
 trm3i=-(fac3a*a2i+fac3b*b2i+(eta01-ee)*c2i)/den3
 IF (test0 <= 0.0001)  GO TO 110
 coef = (2.0*ee)/(r1sqx-e2)
 arga = coef*zet01
 test = ABS(arga)
 IF  (test > 0.3)  GO TO 90
 s    = arga**2
 ser  = 1./3.+s*(-1./5.+s*(1./7.+s*(-1./9.+s*(1./11.-s/13.))))
 alpha= e2*(coef**2)*ser
 funct= coef*(1.0-alpha*(zet01**2)/e2)
 GO TO 100
 90 CONTINUE
 argt = coef*azet
 atana= ATAN(argt)
 funct= atana/azet
 100 CONTINUE
 trm1r= facr*funct
 trm1i= faci*funct
 diijr= (trm1r + trm2r + trm3r)/(2.0*zet02)
 diiji= (trm1i + trm2i + trm3i)/(2.0*zet02)
 GO TO 170
 110 CONTINUE
 funct= 0.0
 GO TO 100
 120 CONTINUE
 dena = (eta01+ee)**2 + zet01**2
 denb = (eta01-ee)**2 + zet01**2
 up1r = 2.0*(e2*a2r + c2r)
 up1i = 2.0*(e2*a2i + c2i)
 up2r = 4.0*e2*eta01*b2r
 up2i = 4.0*e2*eta01*b2i
 trm1r= (up1r *(r1sqx+e2) + up2r )/(dena*denb)
 trm1i= (up1i *(r1sqx+e2) + up2i )/(dena*denb)
 IF  ((azet/ee)  <=  0.001)  GO TO  130
 coef = (2.0*ee)/(r1sqx-e2)
 arga = coef*zet01
 test = ABS(arga)
 IF  (test > 0.3)  GO TO 125
 s    = arga**2
 ser  = 1./3.+s*(-1./5.+s*(1./7.+s*(-1./9.+s*(1./11.-s/13.))))
 alpha= e2*(coef**2)*ser
 funct= coef*(1.0-alpha*(zet01**2)/e2)
 GO TO 140
 125 CONTINUE
 argt= coef*azet
 atana= ATAN(argt)
 funct= atana/azet
 alpha= (e2/zet02)*(1.0-funct*(deno/(2.0*ee)))
 GO TO 140
 130 CONTINUE
 alpha= ((2.0*e2)/(eta02-e2))**2
 140 CONTINUE
 trm2r= -alpha*facr/e2
 trm2i= -alpha*faci/e2
 diijr= ee*(trm1r + trm2r)/deno
 diiji= ee*(trm1i + trm2i)/deno
 170 CONTINUE
 
 RETURN
END SUBROUTINE idf2
