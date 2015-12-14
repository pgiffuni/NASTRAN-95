SUBROUTINE idf1 (ee,e2,    eta01,zet01,are,aim,bre,bim,cre,cim,  &
        r1sqx,xiijr,xiiji)
!   ***   INTEGRATES THE PLANAR PARTS OF THE INCREMENTAL
!         OSCILLATORY KERNELS FOR UNSTEADY CASES
 pi   = 3.1415926
 parn = eta01**2 - zet01**2
 facr = parn*are + eta01*bre + cre
 faci = parn*aim + eta01*bim + cim
 parnr= bre/2.0  + eta01*are
 parni= bim/2.0  + eta01*aim
 up   = (eta01-ee)**2 + zet01**2
 down = (eta01+ee)**2 + zet01**2
 arg2 = up/down
 alarg2 = ALOG(arg2)
 trm2r= parnr * ALOG(arg2)
 trm2i= parni * ALOG(arg2)
 trm3r= 2.0*ee* are
 trm3i= 2.0*ee* aim
 azet = ABS(zet01)
 IF  ((azet/ee)  <=  0.001)  GO TO  100
 test0= ABS((r1sqx-e2)/(2.0*ee*azet))
 IF (test0 <= 0.0001)  GO TO 110
 coef = (2.0*ee)/(r1sqx-e2)
 arga = coef*zet01
 test = ABS(arga)
 IF (test <= 0.3)  GO TO 120
 argt = coef*azet
 atana= ATAN(argt)
 funct= atana/azet
 GO TO 170
 100 CONTINUE
 funct= (2.0*ee)/(eta01**2-e2)
 GO TO 170
 110 CONTINUE
 funct= 0.0
 GO TO 170
 120 CONTINUE
 s    = arga**2
 ser  = 1./3.+s*(-1./5.+s*(1./7.+s*(-1./9.+s*(1./11.-s/13.))))
 alpha= e2*(coef**2)*ser
 funct= coef*(1.0-alpha*(zet01**2)/e2)
 170 CONTINUE
 trm1r= facr * funct
 trm1i= faci * funct
 xiijr= trm1r + trm2r + trm3r
 xiiji= trm1i + trm2i + trm3i
 RETURN
END SUBROUTINE idf1
