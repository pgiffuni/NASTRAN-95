SUBROUTINE amgt1b (q,nstns2,c1sbar,c2sbar)
     
!      SUBSONIC RAO (CASCADES) CODE FOR SWEPT TURBOPROPS.
 
 
 COMPLEX, INTENT(OUT)                     :: q(nstns2,nstns2)
 INTEGER, INTENT(IN)                      :: nstns2
 REAL, INTENT(IN)                         :: c1sbar
 REAL, INTENT(IN)                         :: c2sbar
 INTEGER :: sln
 REAL :: m,kappa,mu,mus,lamda,lamdm,nu,x(20),disp(20,10), w(8),ww(8),m2sbar
 COMPLEX :: loads(21),stt(20),sum1,sum2,stti,  &
     an(401),ab(401),fk(401),cn(401),cb(401),pd(401),  &
     so(100),s1(100),p(50),a(20,40),ff,st,stp,fg,fs,fo, slope
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /tamg1l/ iref,minmac,maxmac,nlines,nstns,refstg,refcrd,  &
     refmac,refden,refvel,refswp,sln,nstnsx,stag,  &
     chord,dcbdzb,bspace,mach,den,vel,sweep,amach, redf,blspc,amachr,tsonic
 COMMON /amgmn / mcb(7),nrow,dum(2),refc,sigma,rfreq
 COMMON /system/ ibuf,iout
 DATA    w     / 1.48283,  .89414, .83521,  .66721,  &
     .64172,  .55519, .54026,  .48547 /
 DATA    ww    / 1.13333,-0.00036,0.18796,-0.00027,  &
     0.08469,-0.00022,0.05049,-0.00019 /
 
!     THEORY DEPENDENT RESTRICTION OF NO MORE THAN 10 COMPUTING
!     STATIONS PER STREAMLINE IS REFLECTED IN CODING.
 
 IF (nstns > 10) GO TO 1000
 
 m      = amach
 omega  = redf
 ss     = 2*blspc
 deltm  =-sigma
 xlam   = stag
 nm     = nstns
 nnm    = nstns2
 csbar  = .25*(den*vel**2*chord**2)/(refden*refvel**2)
 csbar1 = 2.0/chord
 m2sbar =-dcbdzb/chord
 c2ssch = csbar1*c2sbar
 csblsb = csbar*csbar1
 csbm2s = csbar*m2sbar
 n      = 20
 pi     = 3.141593
 pi2    = pi*2
 con    = 1.0E-5
 nnn    = 100
 kkk    = 2*nnn + 1
 deltm  = deltm/360
 xl     = xlam*pi/180
 b      = 1.0/n
 b2     = 2*b
 d      = ss*SIN(xl)
 hh     = ss*COS(xl)
 beta   = SQRT(1. - m**2)
 h      = hh*beta
 zer    = 0.0
 s      = SQRT(h**2 + d**2)
 lamdm  = ATAN(d/h)
 cl     = COS(lamdm)
 sl     = SIN(lamdm)
 nu     = omega/beta**2
 kappa  = m*nu
 lamda  = m*kappa
 delta  = deltm + lamda*d/pi2
 mu     = kappa*s/pi2
 mus    = mu**2
 ff     = (0.0,1.0)
 fg     = CMPLX(zer,nu*s)
 l      = 1
 cc     = delta**2 - mus
 IF (cc == 0.0) GO TO 200
 IF (cc < 0.0) fk(l) = SQRT(-cc)*ff
 IF (cc > 0.0) fk(l) = SQRT(cc)
 an(l)  = fk(l)*cl + ff*delta*sl
 ab(l)  = fk(l)*cl - ff*delta*sl
 pd(l)  = fk(l)*(pi2*ab(l) + fg)
 ck     = pi2*b/s
 cn(l)  = CEXP(-an(l)*ck)
 cb(l)  = CEXP(-ab(l)*ck)
 DO  i = 1,nnn
   l      = l + 1
   cc     = (delta+i)**2 - mus
   IF (cc == 0.0) GO TO 200
   IF (cc < 0.0) fk(l) = SQRT(-cc)*ff
   IF (cc > 0.0) fk(l) = SQRT(cc)
   an(l)  = fk(l)*cl + (delta+i)*ff*sl
   ab(l)  = fk(l)*cl - (delta+i)*ff*sl
   pd(l)  = fk(l)*(pi2*ab(l) + fg)
   cn(l)  = CEXP(-an(l)*ck)
   cb(l)  = CEXP(-ab(l)*ck)
   l      = l + 1
   cc     = (delta-i)**2 - mus
   IF (cc == 0.0) GO TO 200
   IF (cc > 0.0) fk(l) = SQRT(cc)
   IF (cc < 0.0) fk(l) = SQRT(-cc)*ff
   an(l)  = fk(l)*cl + (delta-i)*ff*sl
   ab(l)  = fk(l)*cl - (delta-i)*ff*sl
   pd(l)  = fk(l)*(pi2*ab(l) + fg)
   cn(l)  = CEXP(-an(l)*ck)
   cb(l)  = CEXP(-ab(l)*ck)
 END DO
 stp    = 0.0
 l      = 1
 st     = ((1-cn(l))/an(l) + (1-cb(l))/ab(l))/fk(l)
 DO  i = 2,kkk,2
   l      = i
   st     = ((1-cn(l))/an(l) + (1-cb(l))/ab(l))/fk(l) + st
   l      = l + 1
   st     = ((1-cn(l))/an(l) + (1-cb(l))/ab(l))/fk(l) + st
   IF (cabs(st-stp) < con) EXIT
   stp    = st
 END DO
 30 CONTINUE
 so(1)  =-st*s/(2*pi2*b2)
 DO  j = 2,n
   jk     = 2*(j-1)
   l      = 1
   stp    = 0.0
   st     = cn(l)**jk/fk(l)
   DO  i = 2,kkk,2
     l      = l + 1
     st     = cn(l)**jk/fk(l) + st
     l      = l + 1
     st     = cn(l)**jk/fk(l) + st
     IF (cabs(st-stp) < con) GO TO 35
     stp    = st
   END DO
   35 so(j)  =-0.5*st
 END DO
 n1     = n + 1
 n2     = 3*n - 1
 DO  j = n1,n2
   jk     = j - n
   stp    = 0.0
   l      = 1
   st     = cb(l)**jk/fk(l)
   DO  i = 2,kkk,2
     l      = l + 1
     st     = cb(l)**jk/fk(l) + st
     l      = l + 1
     st     = cb(l)**jk/fk(l) + st
     IF (cabs(st-stp) < con) GO TO 45
     stp    = st
   END DO
   45 so(j)  =-0.5*st
 END DO
 DO  j = 1,n
   jk     = (j-1)*2 + 1
   l      = 1
   stp    = 0.0
   st     = an(l)*cn(l)**jk/fk(l)
   DO  i = 2,kkk,2
     l      = l + 1
     st     = an(l)*cn(l)**jk/fk(l) + st
     l      = l + 1
     st     = an(l)*cn(l)**jk/fk(l) + st
     IF (cabs(st-stp) < con) EXIT
     stp    = st
   END DO
   54 s1(j)  =-pi/s*st
 END DO
 n1     = n + 1
 n2     = 2*n
 DO  j = n1,n2
   jk     = (j-n1)*2 + 1
   l      = 1
   stp    = 0.0
   st     = ab(l)*cb(l)**jk/fk(l)
   DO  i = 2,kkk,2
     l      = l + 1
     st     = ab(l)*cb(l)**jk/fk(l) + st
     l      = l + 1
     st     = ab(l)*cb(l)**jk/fk(l) + st
     IF (cabs(st-stp) < con) EXIT
     stp    = st
   END DO
   59 s1(j)  = pi/s*st
 END DO
 DO  j = 1,n
   jk     = (j-1)*2 + 1
   l      = 1
   stp    = 0.0
   st     = cb(l)**jk/pd(l)
   DO  i = 2,kkk,2
     l      = l + 1
     st     = cb(l)**jk/pd(l) + st
     l      = l + 1
     st     = cb(l)**jk/pd(l) + st
     IF (cabs(st-stp) < con) EXIT
     stp    = st
   END DO
   62 p(j)   =-s/2*st
 END DO
 fg     = CMPLX(zer,-nu*b)
 fg     = 1/(CEXP(fg) + CMPLX(zer,nu*b2))
 fs     = CMPLX(zer,nu)
 cj     = (nu*beta)**2
 l      = 0
 ct     = 2*kappa**2*b
 DO  j = 1,n
   DO  i = 1,n
     l      = l + 1
     nk     = i - j + 1
     nk1    = i - j
     nk2    = nk1 + 1
     IF (i == j) nk1 = n + 1
     IF (i == j) nk2 = 1
     IF (j <= i) GO TO 65
     nk1    = n + j - i + 1
     nk2    = nk1 - 1
     nk     = n + 2*(j-i)
     65 a(i,j) = s1(nk1) - s1(nk2) + ct*so(nk)
     IF (j /= n) CYCLE
     nk     = n + 2*(j-i) + 1
     nk2    = j - i + 1
     a(i,j) = a(i,j) - fg*(s1(nk1) + so(nk)*fs + cj*p(nk2))
   END DO
 END DO
 x(1)   =-1.0 + b
 DO  i = 2,n
   x(i)   = x(i-1) + b2
 END DO
 n1     = n  + nm
 nn     = n1
 nn1    = nn + nm
 n1n    = n  - 1
 n1m    = nm - 1
 n11    = n  + 1
 nn11   = nn + 1
 n22    = n  + 2
 nn22   = nn + 2
 fo     = ff*omega
 tanlam = TAN(sweep*pi/180.0)
 dlsdzb = dcbdzb/2.0
 td     = tanlam*dlsdzb
 DO  i = 1,n
   disp(i,1) =-1.0
   disp(i,2) =-1.0 - x(i)
   stt(i) = CEXP(-ff*lamda*x(i))*pi2/beta
   stti   = stt(i)
   a(i,n11 ) = stti*disp(i,1)*(fo+td)
   a(i,nn11) = stti*disp(i,1)*tanlam
   a(i,n22 ) = stti*(disp(i,2)*(fo+td)-1.)
   a(i,nn22) = stti*disp(i,2)*tanlam
 END DO
 DO  jj = 3,nm
   nf     = n  + jj
   nnf    = nn + jj
   con2   = pi*(jj-2)/2
   DO  i = 1,n
     con    = con2*disp(i,2)
     disp(i,jj) = SIN(con)
     a(i,nf ) = stt(i)*(disp(i,jj)*(fo+td)-con2*COS(con))
     a(i,nnf) = stt(i)*disp(i,jj)*tanlam
   END DO
 END DO
!WKBR SPR93019 10/93 CALL GAUSS (A,N,NN1)
 CALL gauss2 (a,n,nn1)
 DO  j = 1,nnm
   nf     = n + j
   DO  i = 1, n
     loads(i) = a(i,nf)
   END DO
   
   slope   = loads(2)/3./b
   a(1,nf) = 2.*CEXP(lamda*ff*x(1))*(ff*nu*loads(1) + slope)
   
   slope   = (loads(n) - loads(n1n))/b2
   a(n,nf) = 2.*CEXP(lamda*ff*x(n))*(ff*nu*loads(n) + slope)
   
   DO  i = 2,n1n
     slope   = (loads(i+1) - loads(i-1))/4./b
     a(i,nf) = 2.*CEXP(lamda*ff*x(i))*(ff*nu*loads(i) + slope)
   END DO
 END DO
 DO  i = 1,n
   a(i,1) = SQRT((1-x(i))/(1+x(i)))
   DO  j = 2,n1m
     a(i,j) =-disp(i,j+1)
   END DO
   DO  j = nm,n
     con2   =-pi*(j-1)*disp(i,2)/2
     a(i,j) = SIN(con2)
   END DO
 END DO
!WKBR SPR93019 10/93      CALL GAUSS (A,N,NN1)
 CALL gauss2 (a,n,nn1)
 a(1,1) = c2ssch*pi+c1sbar*pi/2.
 a(2,1) = (c2ssch+c1sbar)*pi/2.
 con    = 1.
 conn   = 1.
 DO  j = 1,n1n
   a(1,j+1) = (c2ssch*con+c1sbar*conn)*4./j/pi
   a(2,j+1) = (c2ssch+2.*c1sbar)*conn*4./j/pi-con*c1sbar*32./ (j*pi)**3
   con    = 1. - con
   conn   =-conn
 END DO
 DO  i = 3,nm
   ir     = i - 2
   DO  j = 2,n
     is     = j - 1
     IF (ir == is) GO TO 901
     IF ((ir+is)/2*2 == (ir+is)) GO TO 902
     a(i,j) =-c1sbar*16.*ir*is/(pi*(ir+is)*(ir-is))**2
     CYCLE
     902 a(i,j) = (0.,0.)
     CYCLE
     901 a(i,j) = c2ssch + c1sbar
   END DO
 END DO
 DO  j = 3,nm
   a(j,1) = c2ssch*w(j-2)+c1sbar*ww(j-2)
 END DO
 DO  j = 1,nm
   DO  k = 1,nm
     nf     = n  + k
     nnf    = nn + k
     sum1   = (0.,0.)
     sum2   = (0.,0.)
     do150 i = 1,n
     sum1 = sum1 + a(j,i)*a(i,nf)
     150 sum2 = sum2 + a(j,i)*a(i,nnf)
     q(j,k   ) = csblsb*sum1 + csbm2s*sum2
     q(j,k+nm) = csbar*sum2
     q(j+nm,k) = (0.,0.)
     q(j+nm,k+nm) = (0.,0.)
   END DO
 END DO
 200 RETURN
 
 1000 WRITE  (iout,3001) ufm,sln,nstns
 3001 FORMAT (a23,' - AMG MODULE - NUMBER OF COMPUTING STATIONS ON ',  &
     'STREAMLINE',i8,4H is ,i3,1H., /39X,'SUBSONIC CASCADE ',  &
     'ROUTINE AMGT1B ALLOWS ONLY A MAXIMUM OF 10.')
 CALL mesage (-61,0,0)
 RETURN
END SUBROUTINE amgt1b
