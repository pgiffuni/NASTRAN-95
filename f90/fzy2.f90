SUBROUTINE fzy2 (xij, x1, x2,eta,zeta, yb, zb, a, beta2,cbar, k,  &
        fzzr, fzzi, fzyr, fzyi, fyzr, fyzi, fyyr, fyyi)
!   ***   THIS SUBROUTINE IS AN ALTERNATIVE TO SUBROUTINE  FMZY   ---
!         IT IS USED WHENEVER THE OPTION FLAG   IBFS  =  1
!   ***
 
 REAL, INTENT(IN)                         :: xij
 REAL, INTENT(IN)                         :: x1
 REAL, INTENT(IN)                         :: x2
 REAL, INTENT(IN OUT)                     :: eta
 REAL, INTENT(IN OUT)                     :: zeta
 REAL, INTENT(IN OUT)                     :: yb
 REAL, INTENT(IN OUT)                     :: zb
 REAL, INTENT(IN)                         :: a
 REAL, INTENT(IN)                         :: beta2
 REAL, INTENT(IN)                         :: cbar
 REAL, INTENT(IN)                         :: k
 REAL, INTENT(OUT)                        :: fzzr
 REAL, INTENT(OUT)                        :: fzzi
 REAL, INTENT(OUT)                        :: fzyr
 REAL, INTENT(OUT)                        :: fzyi
 REAL, INTENT(OUT)                        :: fyzr
 REAL, INTENT(OUT)                        :: fyzi
 REAL, INTENT(OUT)                        :: fyyr
 REAL, INTENT(OUT)                        :: fyyi
 REAL :: m, kbar,kbar2,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,kbar3
 DATA lastbr /0/
 DATA    test1,test2,cth,sth,raij,raij2 /0.142857, 0.5, 1.0,3*0.0/
 DATA    capdr,capdi,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11/ 13*0.0 /
 
 m      = SQRT(1.0 - beta2)
 IF  (k <= 0.0001 .AND. m <= 0.0001)  GO TO  110
 kbar   = 2.0 *k *m *a / cbar
 kbar2  = kbar*kbar
 GO TO  120
 110 CONTINUE
 kbar   = 0.0
 kbar2  = 0.0
 120 xa     = 0.5 * (x1 + x2)
 dx     = x2 - x1
 a2     = a * a
 eps    = 0.001 * a2
 IF  (eta == yb .AND. zeta == zb)  GO TO  130
 raij2  = (eta-yb)**2 + (zeta-zb)**2
 raij   = SQRT(raij2)
 cth    = (eta- yb) / raij
 sth    = (zeta-zb) / raij
 IF  (raij2 > a2)  GO TO  150
 GO TO  140
 130 CONTINUE
 raij   = 0.0
 raij2  = 0.0
 cth    = 1.0
 sth    = 0.0
 140 rwig2  = a2
 GO TO 160
 150 rwig2  = raij2
 160 raa    = SQRT((xa -xij)**2 + beta2*rwig2)
 ct2    = cth*cth
 st2    = 0.0
 IF  (ABS(sth) > 0.0001)  st2=sth*sth
 rwig   = SQRT(rwig2)
 raa2   = raa * raa
 raa3   = raa * raa2
 raa4   = raa * raa3
 capa   = m - (xa- xij) / raa
 delta  = dx / raa
 delta2 = delta * delta
 earg   = 0.0
 IF  (kbar <= 0.0001)  GO TO  180
 earg   = kbar * (m * (xa-xij) - raa) / (beta2 * a)
 qr     = COS(earg) / (4.0 * dx)
 qi     = SIN(earg) / (4.0 * dx)
 GO TO  190
 180 qr     = 1.0 / (4.0 * dx)
 qi     = 0.0
 190 CONTINUE
 IF  (delta < test1)  GO TO  240
 i1     = delta / raa2
 trm1   = beta2 * a * i1
 fthr   = a * qr * trm1
 fthi   = a * qi * trm1
 IF  (kbar <= 0.0001)  GO TO  210
 i4     = delta / raa
 trm2   = kbar * i4
 fthr   = fthr - a * qi * trm2
 fthi   = fthi + a * qr * trm2
 210 CONTINUE
 IF  (raij2 < (a2+eps))  GO TO  220
 frr    = fthr
 fri    = fthi
 GO TO  370
 220 i6     = delta / raa4
 trm1   = -3.0 * a2 * beta2*beta2 * i6
 capdr  = raij2 * qr * trm1
 capdi  = raij2 * qi * trm1
 IF  (kbar <= 0.0001)  GO TO  230
 i9     = delta / raa3
 trm1   = trm1 + kbar2 * i1
 trm2   = -3.0 * a * beta2 * kbar * i9
 capdr  = raij2 * (qr * trm1 - qi * trm2)
 capdi  = raij2 * (qr * trm2 + qi * trm1)
 230 frr    = fthr + capdr
 fri    = fthi + capdi
 GO TO  370
 240 CONTINUE
 IF  (delta > test2)  GO TO  320
 lastbr = 0
 tau    = (xa - xij) / raa
 tau2   = tau * tau
 i1     = delta * (1.0 - (-1.0+5.0*tau2)*delta2/8.0) / raa2
 250 trm1   = a * beta2 * i1
 fthr   = a * qr * trm1
 fthi   = a * qi * trm1
 IF  (kbar <= 0.0001)  GO TO  270
 IF  (lastbr /= 0)  GO TO  350
 delta3 = delta * delta2
 i2     = -(tau * delta3) / (4.0 * raa)
 i3     = delta3 / 12.0
 i4     = delta * (1.0 + (-1.0+3.0*tau2)*delta2/12.0) / raa
 i5     = -(tau * delta3) / 6.0
 260 trm1   = trm1 - (kbar2 * capa * i5) / (a * beta2)
 trm2   = kbar * (capa * i2 + i4 - i3*beta2*rwig2/(2.0*raa3) )
 fthr   = a * (qr * trm1 - qi * trm2)
 fthi   = a * (qr * trm2 + qi * trm1)
 270 IF  (raij2 > (a2+eps))  GO TO  280
 frr    = fthr
 fri    = fthi
 GO TO  370
 280 CONTINUE
 kbar3  = kbar*kbar2
 IF  (lastbr /= 0)  GO TO  340
 i6     = delta * (1.0 + 5.0*(-1.0+7.0*tau2)*delta2/24.0) / raa4
 290 trm1   = -3.0 * a2 * beta2*beta2 * i6
 capdr  = raij2 * qr * trm1
 capdi  = raij2 * qi * trm1
 IF  (kbar <= 0.0001)  GO TO  310
 IF  (lastbr /= 0)  GO TO  360
 i7     = -5.0 * tau * delta3 / (12.0 * raa3)
 i8     = delta3 / (12.0 * raa2)
 i9     = delta * (1.0 + (-1.0+6.0*tau2)*delta2/6.0) / raa3
 i10    = -delta3 * tau / (3.0 * raa2)
 300 trm1   = trm1 + kbar2 * (i1 + 3.0 * capa * i10)
 trm2   = 3.0*a*beta2 * kbar * (-capa*i7 +i8*beta2*rwig2/(2.0*raa3)  &
     -i9) + kbar3 * capa * i2 / (a * beta2)
 capdr  = raij2 * (qr * trm1 - qi * trm2)
 capdi  = raij2 * (qr * trm2 + qi * trm1)
 310 frr    = fthr + capdr
 fri    = fthi + capdi
 GO TO  370
 320 CONTINUE
 lastbr = 1
 rwig   = SQRT(rwig2)
 ra12   = (x1 - xij)**2 + beta2 * rwig2
 ra22   = (x2 - xij)**2 + beta2 * rwig2
 ra1    = SQRT(ra12)
 ra2    = SQRT(ra22)
 i1     = ((x2-xij)/ra2 - (x1-xij)/ra1) / (beta2*rwig2)
 GO TO  250
 340 CONTINUE
 ra13   = ra1 * ra12
 ra23   = ra2 * ra22
 i6     = ((x2-xij)/ra23-(x1-xij)/ra13 + 2.0*i1)/(3.0*beta2*rwig2)
 GO TO  290
 350 part1  = 0.5 * dx * (xa - xij)
 i2     = -((part1+raa2)/ra2 + (part1-raa2)/ra1)/ (beta2*rwig2)
 denom  = x1 - xij + ra1
 i11    = ALOG(ABS((x2 - xij + ra2) / denom))
 i3     = i11 - 2.0*(xa - xij)*i2 - raa2 * i1
 deno4  = SQRT(beta2) * rwig
 arg1   = (x2 - xij) / deno4
 arg2   = (x1 - xij) / deno4
 i4     = (ATAN(arg1) - ATAN(arg2)) / deno4
 i5     = 0.5 * ALOG(ra22 / ra12) - (xa - xij) * i4
 GO TO  260
 360 CONTINUE
 i7     = -(1.0/ra23 - 1.0/ra13) / 3.0 - (xa - xij) * i6
 i8     = i1 - 2.0 * (xa - xij) * i7 - raa2 * i6
 i9     = ((x2-xij)/ra22-(x1-xij)/ra12 + i4) / (2.0*beta2*rwig2)
 i10    = -((part1 + raa2)/ra22 + (part1 - raa2)/ra12 +  &
     (xa - xij) * i4) / (2.0 * beta2 * rwig2)
 GO TO 300
 370 CONTINUE
 fzzr   = ct2 * fthr + st2 * frr
 fzzi   = ct2 * fthi + st2 * fri
 fyyr   = st2 * fthr + ct2 * frr
 fyyi   = st2 * fthi + ct2 * fri
 IF  (cth == 0.0 .OR. sth == 0.0)  GO TO  400
 fzyr   = cth * sth * (frr - fthr)
 fzyi   = cth * sth * (fri - fthi)
 GO TO  410
 400 fzyr   = 0.0
 fzyi   = 0.0
 410 CONTINUE
 fyzr = fzyr
 fyzi = fzyi
 RETURN
END SUBROUTINE fzy2
