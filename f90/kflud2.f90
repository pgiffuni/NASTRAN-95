SUBROUTINE kflud2
     
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
 
 LOGICAL :: nogo
 INTEGER :: out,eltype,necpt(100)
 DOUBLE PRECISION :: r1,z1,r2,z2,constd,dpi,  &
     z1p,z2p,z1p1,z2p1,rk,ri,kfact,f0,a,b,i2n0,i2n1,  &
     i2n2,i2np2,dz,hpq,pirho,twopr,kh,k1,k2
 CHARACTER (LEN=23) :: ufm
 COMMON  /xmssg / ufm
 COMMON  /condad/ constd(5)
 COMMON  /system/ sysbuf,out,nogo
 COMMON  /emgdic/ eltype
 COMMON  /sma1io/ dum1(10),ifkgg
 COMMON  /sma1cl/ iopt4,k4ggsw,npvt
 COMMON  /sma1dp/ z1p,z2p,rk,ri,kfact,f0,a,b,i2n0,i2n1,i2n2,i2np2,  &
     dz,hpq(4),pirho,twopr,kh(4),k1,k2
 COMMON  /sma1et/ ecpt(100)
 EQUIVALENCE      (constd(1),dpi),(ecpt(1),necpt(1))
 
 
 IF (ecpt(13)-ecpt(9) < 0.0) THEN
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
 15 IF (r1 == 0.0D0 .OR. r2 == 0.0D0) GO TO 5000
 IF (z1 == z2) RETURN
 
!     CALCULATE THE INTEGRAL PARAMETERS I2N0,I2N1,I2N2,AND I2NP2
 
 k  = 2*necpt(6)
 rk = k
 IF (k > 0) GO TO 20
 
 i2n0 = 0.0
 i2n1 = 0.0
 i2n2 = 0.0
 i2np2= (z2-z1)*(r2**2 + r2*r1 + r1**2)/6.0D0
 
 GO TO 300
 
 20 b    = (r2-r1)/(z2-z1)
 dum  = DABS(b)
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
 i2n1 = a/(rk*(rk+1.0D0))*(z2p*z2-z1p*z1-a/(rk+2.0D0)*(z2p1-z1p1))
 i2n2 = a/(rk*(rk+1.0D0))*(z2p*z2**2-z1p*z1**2 -a/(rk+2.0D0)*2.0D0  &
     * (z2p1*z2-z1p1*z1-a/(rk+3.0D0)*(z2p1*r2-z1p1*r1)))
 i2np2= a/((rk+2.0D0)*(rk+3.0D0))*(z2p1*r2-z1p1*r1)
 300 dz   = z2 - z1
 n    = necpt(6)
 z1p  = r1**n
 z2p  = r2**n
 hpq(1) = z2/(dz*z1p)
 hpq(2) =-z1/(dz*z2p)
 hpq(3) =-1.0D0/(dz*z1p)
 hpq(4) = 1.0D0/(dz*z2p)
 lp   = 1
 IF (npvt == necpt(2)) GO TO 320
 IF (npvt == necpt(3)) GO TO 310
 RETURN
 
 310 lp = 2
 320 IF (ecpt(4) == 0.0) RETURN
 pirho  = dpi/DBLE(ecpt(4))
 IF (n == 0) pirho = pirho*2.0D0
 rk = n
 twopr = 2.0*pirho*rk**2
 kh(1) = twopr*(i2n0*hpq(lp)+i2n1*hpq(lp+2))
 kh(2) = twopr*(i2n1*hpq(lp)+i2n2*hpq(lp+2)) +pirho*i2np2*hpq(lp+2)
 k1    = kh(1)*hpq(1) + kh(2)*hpq(3)
 k2    = kh(1)*hpq(2) + kh(2)*hpq(4)
 ifile = ifkgg
 i     = npvt
 j     = necpt(2)
 CALL sma1b (k1,j,i,ifile,0.0D0)
 j     = necpt(3)
 CALL sma1b (k2,j,i,ifile,0.0D0)
 RETURN
 
 5000 n = necpt(1)
 IF (eltype == 43) n = n/1000
 WRITE  (out,6000) ufm,n
 6000 FORMAT (a23,' 5000, NEGATIVE OR ZERO RADIUS DETECTED FOR ',  &
     'CFLUID2/CAXIF2 ELEMENT ID',i9)
 nogo = .true.
 RETURN
END SUBROUTINE kflud2
