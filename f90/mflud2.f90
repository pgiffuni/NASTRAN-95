SUBROUTINE mflud2
     
!     THIS ROUTINE GENERATES THE PSUEDO STIFFNESS MATRIX TERMS
!     FOR THE CENTER PLUG FLUID ELEMENT
 
!     THE ECPT DATA BLOCK CONTAINS THE FOLLOWING DATA
 
!         FIELD    SYMBOL
!           1        ID
!           2        SIL1
!           3        SIL2
!           4        RHO
!           5        BULK
!           6        N
!           7        CSF
!           8        R1
!           9        Z1
!           10       -
!           11       CSF
!           12       R2
!           13       Z2
!           14       -
!           15       -
 
 INTEGER :: necpt(100)
 DOUBLE PRECISION :: constd,dpi,r1,z1,r2,z2,z1p,z2p,z1p1,z2p1,rk,ri,  &
     kfact,f0,a,b,i2n0,i2n1,i2n2,i2np2,dz,hpq,pirho, twopr,kh,k1,k2
 COMMON /condad/  constd(5)
 COMMON /sma2et/  ecpt(100)
 COMMON /sma2io/  dum1(10),ifmgg
 COMMON /sma2cl/  dun2(2),npvt
 COMMON /sma2dp/  z1p,z2p,rk,ri,kfact,f0,a,b,i2n0,i2n1,i2n2,i2np2,  &
     dz,hpq(4),pirho,twopr,kh(4),k1,k2
 EQUIVALENCE      (constd(1),dpi),(ecpt(1),necpt(1))
 
 
 IF (ecpt(13) - ecpt(9) < 0.0) THEN
   GO TO     5
 ELSE
   GO TO    10
 END IF
 5 r1 = ecpt(12)
 r2 = ecpt(8)
 z1 = ecpt(13)
 z2 = ecpt(9)
 i  = necpt(3)
 necpt(3) = necpt(2)
 necpt(2) = i
 GO TO 15
 10 r1 = ecpt(8)
 z1 = ecpt(9)
 r2 = ecpt(12)
 z2 = ecpt(13)
 15 IF (ecpt(5) <= 0.0) RETURN
 IF (r1 == 0.0 .OR. r2 == 0.0) GO TO 350
 IF (z1 == z2) GO TO 350
 
!     CALCULATE THE INTEGRAL PARAMETERS I2N0,I2N1,I2N2,AND I2NP2
 
 k  = 2*necpt(6) + 2
 rk = k
 b   = (r2-r1)/(z2-z1)
 dum = DABS(b)
 IF (dum > 1.0E-6) GO TO 30
 z1p  = ((r1+r2)/2.0D0)**k
 i2n0 = (z1p/rk)*(z2-z1)
 i2n1 = i2n0*(z2+z1)/2.0D0
 i2n2 = i2n0*(z2**2+z2*z1+z1**2)/3.0D0
 i2np2= i2n0*rk/(rk+2.0D0)*r1**2
 GO TO 300
 
 30 z1p  = r1**(k+1)
 z2p  = r2**(k+1)
 z1p1 = z1p*r1
 z2p1 = z2p*r2
 a    = 1.0D0/b
 i2n0 = a/(rk*(rk+1.0D0))*(z2p-z1p)
 i2n1 = a/(rk*(rk+1.0D0))*(z2p*z2-z1p*z1 -a/(rk+2.0D0)*(z2p1-z1p1))
 i2n2 = a/(rk*(rk+1.0D0))*(z2p*z2**2 -z1p*z1**2 -a/(rk+2.0D0)*2.0D0  &
     * (z2p1*z2 -z1p1*z1 -a/(rk+3.0D0)*(z2p1*r2-z1p1*r1)))
 i2np2= a/((rk+2.0D0)*(rk+3.0D0))*(z2p1*r2-z1p1*r1)
 
 300 dz   = z2 - z1
 n    = necpt(6)
 z1p  = r1**n
 z2p  = r2**n
 hpq(1) = z2/(dz*z1p)
 hpq(2) =-z1/(dz*z2p)
 hpq(3) =-1.0D0/(dz*z1p)
 hpq(4) = 1.0D0/(dz*z2p)
 lp = 1
 IF (npvt == necpt(2)) GO TO 320
 IF (npvt == necpt(3)) GO TO 310
 GO TO 350
 310 lp = 2
 320 pirho = dpi/DBLE(ecpt(5))
 IF (necpt(6) == 0) pirho = 2.0D0*pirho
 kh(1) = pirho*(i2n0*hpq(lp)+i2n1*hpq(lp+2))
 kh(2) = pirho*(i2n1*hpq(lp)+i2n2*hpq(lp+2))
 k1 = kh(1)*hpq(1) + kh(2)*hpq(3)
 k2 = kh(1)*hpq(2) + kh(2)*hpq(4)
 ifile = ifmgg
 i  = npvt
 j  = necpt(2)
 CALL sma2b (k1,j,i,ifile,0.0D0)
 j  = necpt(3)
 CALL sma2b (k2,j,i,ifile,0.0D0)
 350 RETURN
END SUBROUTINE mflud2
