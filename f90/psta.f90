SUBROUTINE psta(deltay,bi,ca,alph,thi,ajjl)
     
 REAL, INTENT(IN)                         :: deltay(1)
 REAL, INTENT(IN)                         :: bi(1)
 REAL, INTENT(IN OUT)                     :: ca(1)
 REAL, INTENT(IN)                         :: alph(1)
 REAL, INTENT(IN)                         :: thi(13)
 REAL, INTENT(IN OUT)                     :: ajjl
 DIMENSION a(3,3),ai(6),aj(6),h(3,3),ek(6),g(3,3),gi(3,3),q(3,3)
 
 DIMENSION p(3,6),qi(3,3)
 COMPLEX :: pc(3)
 COMMON /packx/ iti,it0,ii,nn,incr
 COMMON /amgmn/ mcb(7),nrow,nd,NE,refc,emach,rfk
 COMMON /pstonc/ njj,nmach,nthry,nthick,nalpha,nxis,ntaus,nstrip, seclam
 COMMON / condas / pi,twopi,radg,degra
 DATA a /9*0.0/ , h /9*0.0/
 
 bref = refc * .5
 rfc = rfk/bref
 ii = nrow +1
 nn = nrow
 
!     BUILD AJJL FOR EACH STRIP
 
 DO  i=1,nstrip
   b = bi(i)
   const = 8.0 * deltay(i) * (rfc * b)**2
   a(1,1) = -1.0
   a(2,1) = -.5 * b
   a(2,2) = b
   a(3,3) = b
   h(1,1) = -1.0
   h(1,2) = a(2,1)
   h(2,2) = b
   h(3,3) = b
   alpha = alph(1)
   IF(nalpha /= 1 ) alpha = alph(i)
   alpha = alpha * degra
   alpha2 = alpha*alpha
   n = 2
   IF(ca(i) /= 0.0 ) n = 3
   IF( nthick == 0 ) GO TO 20
   DO  j=1,6
     ai(j) = thi(j)
     aj(j) = 0.0
     IF(n == 3 ) aj(j) = thi(j+6)
   END DO
   zetah = 1.0
   IF( nxis == 1 ) zetah = thi(13)
   IF( nxis > 1) zetah = thi(i+12)
   GO TO 70
   20 IF( ntaus /= 1 ) GO TO 30
   tau = thi(1)
   tauh = thi(2)
   taut = thi(3)
   IF( n == 2 ) taut = 0.
   t = tauh-taut
   zetam = thi(4)
   zetah = thi(5)
   GO TO 50
   30 k = (i-1) * 3+1
   tau = thi(k)
   tauh = thi(k+1)
   taut = thi(k+2)
   IF( n == 2 ) taut = 0.
   t = tauh - taut
   k = (i-1)*2+1 + 3*nstrip
   zetam = thi(k)
   zetah = thi(k+1)
   DO  j=1,6
     aj(j) = 0.
   END DO
   50 IF(n == 2 ) zetah = 1.0
   IF( n == 2 ) GO TO 60
   aj(1) = -.5 * t
   aj(2) = -.25*t*(1.0+zetah)
   aj(3) = -(1./6.)*t*(1.+zetah+zetah*zetah)
   aj(4) = .25*t*t/ (1.-zetah)
   aj(5) = .125*t*t*(1.0+zetah) / (1.-zetah)
   aj(6) = (1./12.)*t*t*(1.+zetah+zetah*zetah) / (1.-zetah)
   60 ts = tau-tauh*(tau-tauh)
   ai(1) = tauh*.5 + aj(1)
   ai(2) = -(tau/3.)*zetah + (tauh/6.) * (2.*zetah+zetam) + aj(2)
   ai(3) = -(tau/12.)*zetah *(3.*zetah+2.*zetam) + (tauh/12.) *  &
       (3. *zetah*zetah + 2.*zetah*zetam+zetam*zetam) + aj(3)
   ai(4) = (tau*tau/(3. *zetam)) + (1./3.)*ts*(zetah-zetam) + aj(4)
   ai(5) = (tau*tau/12.) + (1./12.)*ts*(3.*zetah + zetam) /  &
       (zetah-zetam) + aj(5)
   ai(6) = (tau*tau/30.) * zetam + (1./30.)*ts*(6.*zetah*zetah +  &
       3.*zetah*zetam + zetam*zetam) / (zetah-zetam) + aj(6)
   70 ems = emach*emach
   secs = seclam*seclam
   IF( nthry /= 0 ) GO TO 80
   cbar1 = 1.
   cbar2 = (1.4+1.)/4.
   GO TO 90
   80 cbar1 = emach / SQRT(ems-secs)
   cbar2 = (ems*ems*(1.4+1.)- 4.*secs*(ems-secs)) /(4.*(ems-secs)**2)
   90 cbar3 = (1.4+1.) / 12.
   ek(1) = (1./emach ) *(cbar1+2.*cbar2*emach*ai(1)  &
       + 3.*cbar3*ems*(ai(4)+alpha2))
   ek(2) = (1./emach) * (cbar1+4.*cbar2*emach * ai(2)  &
       + 3.*cbar3*ems*(2.*ai(5)+alpha2))
   ek(3) = (4./(3.*emach)) * (cbar1+6.*cbar2*emach*ai(3)  &
       + 3.*cbar3*ems*(3.*ai(6)+alpha2))
   IF( n == 3 ) GO TO 100
   ek(4) = (1./emach) * (cbar1*(1.-zetah) + 2.*cbar2*emach*aj(1)  &
       +3.*cbar3*ems*aj(4) + alpha2*(1.-zetah))
   ek(5) = (1./emach) * (cbar1*(1.-zetah*zetah) + 4.*cbar2*emach*  &
       aj(2) + 3.*cbar3*ems*(2.*aj(5)+alpha2*(1.-zetah*zetah)))
   ek(6) = (4./(3.*emach))*(cbar1*(1.-zetah**3) + 6.*cbar2*emach*  &
       aj(3) + 3.*cbar3*ems*(3.*aj(6)+ alpha2*(1.-zetah**3)))
   e1k = 1.0/(rfc *b)
   e1ks = e1k*e1k
   g(1,1) = 0.
   g(1,2) = -ek(1) * e1ks
   g(2,1) = 0.
   g(2,2) = -ek(2) * e1ks
   gi(1,1) = -ek(1) * e1k
   gi(1,2) = -ek(2) * e1k
   gi(2,1) = gi(1,2)
   gi(2,2) = -ek(3) * e1k
   IF( n == 3 ) GO TO 100
   g(1,3) = -ek(4) * e1ks
   g(2,3) = -ek(5) * e1ks
   g(3,1) = 0.
   g(3,2) =-(ek(5)-2.*ek(4)*zetah) * e1ks
   g(3,3) = g(3,2)
   gi(1,3) = -(ek(5)-2.*ek(4)*zetah) * e1k
   gi(2,3) = -(ek(6)-2.*ek(5)*zetah) * e1k
   gi(3,1) = gi(1,3)
   gi(3,2) = -(ek(6) -2.*ek(5)*zetah) * e1k
   gi(3,3) = -(ek(6)-4.*ek(5)*zetah+4.*ek(4)*zetah*zetah) * e1k
   
!     MATRICES BUILT TIME TO MULTIPLY
   
   100 DO  k=1,n
     DO  l=1,n
       q(k,l) = 0.
       qi(k,l) = 0.
       DO  m1=1,n
         q(k,l) = q(k,l) + a(k,m1) * g(m1,l)
         qi(k,l) = qi(k,l) + a(k,m1) * gi(m1,l)
       END DO
     END DO
   END DO
   n2 = 2*n
   DO  k=1,n
     DO  l=1,n2,2
       it = l/2+1
       p(k,l) = 0.
       p(k,l+1) = 0.
       DO  m1=1,n
         p(k,l) = p(k,l) + q(k,m1) * h(m1,it)
         p(k,l+1) = p(k,l+1) + qi(k,m1)*h(m1,it)
       END DO
       p(k,l) = p(k,l) * const
       p(k,l+1) = p(k,l+1) * const
     END DO
   END DO
   
!     PACK OUT
   
   nn = nn+n
   DO  j=1,n2,2
     DO  k=1,n
       pc(k) = CMPLX(p(k,j),p(k,j+1))
     END DO
     CALL pack(pc,ajjl,mcb)
   END DO
   ii = ii+n
 END DO
 RETURN
END SUBROUTINE psta
