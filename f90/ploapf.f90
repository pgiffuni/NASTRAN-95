SUBROUTINE ploapf (ecpt,iecpt,l,pa,pb)
     
!     THIS ROUTINE IS CALLED ONLY BY PLOAD1 FOR HANDLING PIN FLAGS OF
!     THE CBAR
 
 
 REAL, INTENT(IN)                         :: ecpt(33)
 INTEGER, INTENT(IN)                      :: iecpt(9)
 REAL, INTENT(IN)                         :: l
 REAL, INTENT(IN OUT)                     :: pa(1)
 REAL, INTENT(IN OUT)                     :: pb(1)
 REAL :: i1,i2,i12,k1,k2, l2,l3,ke,kep,lb,lr1,lr2, l2b3,l2b6
 DIMENSION  pe(12),pep(12), ipin(10),ke(144),kep(144)
 COMMON /matout/ e,g
 
 ka = iecpt(8)
 kb = iecpt(9)
 IF (ka == 0 .AND. kb == 0) GO TO 200
 DO  i = 1,6
   pe(i  ) = pa(i)
   pe(i+6) = pb(i)
 END DO
 l2  = l**2
 l3  = l2*l
 a   = ecpt(17)
 i1  = ecpt(18)
 i2  = ecpt(19)
 fj  = ecpt(20)
 k1  = ecpt(31)
 k2  = ecpt(32)
 i12 = ecpt(33)
 ei1 = e*i1
 ei2 = e*i2
 r1  = 12.0*ei1/l3
 r2  = 12.0*ei2/l3
 IF (k1 == 0.0 .OR. i12 /= 0.0) GO TO 20
 gak = g*a*k1
 r1  = (12.0*ei1*gak)/(gak*l3 + 12.0*l*ei1)
 20 IF (k2 == 0.0 .OR. i12 /= 0.0) GO TO 30
 gak = g*a*k2
 r2  = (12.0*ei2*gak)/(gak*l3 + 12.0*l*ei2)
 
!     COMPUTE THE -SMALL-K-S. SK1, SK2, SK3 AND SK4
 
 30 sk1 = 0.25*r1*l2 + ei1/l
 sk2 = 0.25*r2*l2 + ei2/l
 sk3 = 0.25*r1*l2 - ei1/l
 sk4 = 0.25*r2*l2 - ei2/l
 
!     COMPUTE THE 12 X 12 MATRIX KE
 
 ael = a*e /l
 lr1 = l*r1/2.0
 lr2 = l*r2/2.0
 gjl = g*fj/l
 
 DO  i = 1,144
   ke(i) = 0.0
 END DO
 ke(  1) = ael
 ke(  7) =-ael
 ke( 14) = r1
 ke( 18) = lr1
 ke( 20) =-r1
 ke( 24) = lr1
 ke( 27) = r2
 ke( 29) =-lr2
 ke( 33) =-r2
 ke( 35) =-lr2
 ke( 40) = gjl
 ke( 46) =-gjl
 ke( 51) =-lr2
 ke( 53) = sk2
 ke( 57) = lr2
 ke( 59) = sk4
 ke( 62) = lr1
 ke( 66) = sk1
 ke( 68) =-lr1
 ke( 72) = sk3
 ke( 73) =-ael
 ke( 79) = ael
 ke( 86) =-r1
 ke( 90) =-lr1
 ke( 92) = r1
 ke( 96) =-lr1
 ke( 99) =-r2
 ke(101) = lr2
 ke(105) = r2
 ke(107) = lr2
 ke(112) =-gjl
 ke(118) = gjl
 ke(123) =-lr2
 ke(125) = sk4
 ke(129) = lr2
 ke(131) = sk2
 ke(134) = lr1
 ke(138) = sk3
 ke(140) =-lr1
 ke(144) = sk1
 IF (i12 == 0.0) GO TO 50
 beta    =-12.0*e*i12/l3
 lb      = l *beta/2.0
 l2b3    = l2*beta/3.0
 l2b6    = l2*beta/6.0
 ke( 15) =-beta
 ke( 17) = lb
 ke( 21) = beta
 ke( 23) = lb
 ke( 26) =-beta
 ke( 30) =-lb
 ke( 32) = beta
 ke( 36) =-lb
 ke( 50) = lb
 ke( 54) = l2b3
 ke( 56) =-lb
 ke( 60) = l2b6
 ke( 63) =-lb
 ke( 65) = l2b3
 ke( 69) = lb
 ke( 71) = l2b6
 ke( 87) = beta
 ke( 89) =-lb
 ke( 93) =-beta
 ke( 95) =-lb
 ke( 98) = beta
 ke(102) = lb
 ke(104) =-beta
 ke(108) = lb
 ke(122) = lb
 ke(126) = l2b6
 ke(128) =-lb
 ke(132) = l2b3
 ke(135) =-lb
 ke(137) = l2b6
 ke(141) = lb
 ke(143) = l2b3
 
!     SET UP THE IPIN ARRAY
 
 50 DO  i = 1,5
   ipin(i  ) = MOD(ka,10)
   ipin(i+5) = MOD(kb,10) + 6
   IF (ipin(i+5) == 6) ipin(i+5) = 0
   ka = ka/10
   kb = kb/10
 END DO
 
!     ALTER KE MATRIX DUE TO PIN FLAGS
 
 DO  i = 1,10
   ip = ipin(i)
   IF (ip == 0) CYCLE
   ii = ip*13 - 12
   IF (ke(ii) /= 0.0) GO TO 80
   il = ip
   ii = ii - il
   DO  j = 1,12
     ii = ii + 1
     ke(ii) = 0.0
     ke(il) = 0.0
     il = il + 12
   END DO
   CYCLE
   80 ip12 = (ip-1)*12
   DO  j = 1,12
     j12 = (j-1)*12
     ji  = j12  + ip
     ij  = ip12 + j
     DO  ll = 1,12
       jll = j12  + ll
       ill = ip12 + ll
       kep(jll) = ke(jll) - (ke(ill)/ke(ii))*ke(ji)
     END DO
     pep(j  ) = pe(j  ) - (ke(ji )/ke(ii))*pe(ip)
     kep(ij ) = 0.0
     kep(ji ) = 0.0
   END DO
   pep(ip ) = 0.0
   DO  k = 1,144
     ke(k) = kep(k)
   END DO
   DO  k = 1,12
     pe(k) = pep(k)
   END DO
 END DO
 
 DO  i = 1,10
   ip = ipin(i)
   IF (ip == 0) CYCLE
   pe(ip) = 0.0
 END DO
 DO  i = 1,6
   pa(i) = pe(i  )
   pb(i) = pe(i+6)
 END DO
 
 200 RETURN
END SUBROUTINE ploapf
