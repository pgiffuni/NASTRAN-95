SUBROUTINE fa1pke (khh,bhh,mhh,bxhh,fsave,nloop,bref,rref,neiw,  &
        eps)
     
!     FA1PKE COMPUTES THE EIGENVALUES FOR THE PK METHOD
 
!     LAST REVISED  2/91, BY J.PETKAS/LOCKHEED
!     ELEMENTS OF INTERPOLATION MATRIX IN D.P. AND LEAST SQUARE FIT
 
 
 INTEGER, INTENT(IN)                      :: khh
 INTEGER, INTENT(IN)                      :: bhh
 INTEGER, INTENT(IN)                      :: mhh
 INTEGER, INTENT(IN OUT)                  :: bxhh
 INTEGER, INTENT(IN)                      :: fsave
 INTEGER, INTENT(IN)                      :: nloop
 REAL, INTENT(IN)                         :: bref
 REAL, INTENT(IN OUT)                     :: rref
 INTEGER, INTENT(OUT)                     :: neiw
 REAL, INTENT(OUT)                        :: eps
 LOGICAL :: eigv
 INTEGER :: sysbuf,NAME(2),trl(7),buf1,floop
 REAL :: kint
 DOUBLE PRECISION :: dx1,dx2,dsum,dz(1)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm,uwm,uim,sfm
 COMMON /system/  sysbuf,nout
 COMMON /unpakx/  iout,inn,nnn,incr1
 COMMON /zzzzzz/  z(1)
 COMMON /fa1pkc/  ncore,nk,imvr,ik,ia,iq,icp,iflag
 COMMON /BLANK /  floop
 COMMON /condas/  pi,twopi
 EQUIVALENCE      (z(1),dz(1))
 DATA    NAME  /  4HFA1P,4HKE  /
 DATA    istart/  0 /
 
!     REINITIALIZE EVERY TIME MACH CHANGES
 
 IF (iflag == 0) GO TO 100
 CALL sswtch (39,l39)
 trl(1) = khh
 CALL rdtrl (trl)
 nrow  = trl(2)
 neiw  = MIN0(neiw,nrow)
 neign = nrow*2
 iout  = 1
 inn   = 1
 incr1 = 1
 nnn   = nrow
 ieigns= ncore  - nrow*5 - 1
 buf1  = ieigns - sysbuf
 nn    = nrow*nrow
 nn2   = nn*2
 imh   = icp
 ibh   = imh + nn
 ikh   = ibh + nn
 iv    = ikh + nn
 ib    = iv  + nn
 ima   = ib  + nn
 IF (MOD(ima,2) == 0) ima = ima + 1
 iop   = ima + nn2*4
 
!     CORE CHECK
 
 IF (iop+sysbuf > ieigns) CALL mesage (-8,0,NAME)
 
!     PUT K B M IN CORE
 
 ifl = khh
 ji  = ikh
 10 CALL gopen (ifl,z(buf1),0)
 DO  i = 1,nrow
   CALL unpack (*15,ifl,z(ji))
   GO TO 20
   15 CALL zeroc (z(ji),nrow)
   20 ji = ji + nrow
 END DO
 CALL CLOSE (ifl,1)
 IF (ifl == mhh) GO TO 40
 IF (ifl == bhh) GO TO 30
 ifl = bhh
 ji  = ibh
 trl(1) = bhh
 CALL rdtrl (trl)
 IF (trl(1) > 0) GO TO 10
 CALL zeroc (z(ji),nn)
 30 ifl = mhh
 ji  = imh
 GO TO 10
 40 CONTINUE
 
!     MODIFICATION FOR LEVEL 17.7 UPDATE
!     REPLACE CALLS TO INVAER WITH CALLS TO INVERS.
!     NO NEED TO COMPUTE DETERMINANT SINCE IT IS NOT USED SUBSEQUENTLY.
 
 ising = -1
 CALL invers (nrow,z(imh),nrow,0,0,det,ising,z(iop))
 IF (ising == 2) CALL mesage (-7,0,NAME)
 
!     START OF LARGE LOOP WITH K = 0.0
 
 100 kint = 0.0
 IF (eps <= 0.0) eps = .001
 kn1 = nk + 1
 iq0 = ieigns - (nn2+1)
 ic0 = iq0 - kn1*2 - 2
 IF (MOD(ic0,2) == 0) ic0 = ic0 - 1
 ip0 = ic0 - kn1*2 - 2
 IF (MOD(ip0,2) == 0) ip0 = ip0 - 1
 
 i   = (floop-1)*3
 eigv= .false.
 IF (z(imvr+i+1) < 0.0) eigv = .true.
 vel = ABS(z(imvr+i+1))
 vels= vel*vel
 rho = (rref*z(imvr+i+2))/2.0
 IF (l39 /= 0) WRITE (nout,105) floop,z(imvr+i),z(imvr+i+1), z(imvr+i+2)
 105 FORMAT ('0 TRACE FOR PK METHOD LOOP',i5,6X,4HMACH,8X,  &
     8HVELOCITY,8X,7HDENSITY,/,30X,1P,e15.5,1P,e15.5,1P,e15.5)
 nit   = 0
 nroot = 0
 
!     INITIALIZE LEAST SQUARE COEFFCIENTS
 
 xav   = 0.
 yav1  = 0.
 x10   = 0.
 x11   = 0.
 x12   = 0.
 y10   = 0.
 y11   = 0.
 
!     BUILD P
 
 110 CONTINUE
 nit = nit + 1
 
!     SUM LEAST SQUARE COEFFICIENTS ASSOCIATED WITH INDEPENDENT
!     VARIABLE STARTING WITH SECOND TRIAL
 
 IF (nit == 1) GO TO 115
 xav = xav + kint
 x10 = x10 + 1.
 x11 = x11 + kint
 x12 = x12 + kint**2
 
 115 ip0d = ip0/2 + 1
 dx1  = kint
 DO  i = 1,nk
   dx2 = z(ik+i-1)
   dz(ip0d+i) = DABS((dx1-dx2)**3) + (dx1+dx2)**3
 END DO
 dz(ip0d+kn1) = 1.d0
 
!     FIND C = A-1  P
 
 iad  = ia/2  + 1
 ic0d = ic0/2 + 1
 l = iad
 DO  i = 1,kn1
   dsum = 0.d+0
   DO  j = 1,kn1
     dsum = dsum + dz(l)*dz(ip0d+j)
     l = l + 1
   END DO
   dz(ic0d+i) = dsum
 END DO
 
!     FIND QR AND QI = Q  C  Q IS COLUMN STORED
 
 l = iq
 DO  i = 1,nn2
   dsum = 0.d+0
   DO  j = 1,nk
     k = l + (j-1)*nn2
     dsum = dsum + z(k)*dz(ic0d+j)
   END DO
   l = l + 1
   z(iq0+i) = dsum
 END DO
 
!     COLUMN STORED M-1  BHH  KNH  QR (Z(IQ0+1)   QI (Z(IQ3+NN+1)
 
!     B  =  -BHH  + RHO*BREF*VEL  QHHI
 
!     K  =  -KHH  + RHO*VELS      QHHR
 
!     BUILD  A
!                  0         I
 
!                   -1       -1
!                 -M K     -M B
 
 nrem = iq0 - iop
 IF (nrem-nn <= 0) CALL mesage (-8,0,NAME)
 it = iop
 IF (MOD(it,2) == 0) it = it + 1
 IF (eigv .AND. it+nn > buf1) CALL mesage (-8,0,NAME)
 bov = bref/vel
 rbv = rho*bref*vel
 iqr = iq0
 iqi = iq0 + nn
 rvs = rho*vels
 
!     BUILD M-1K IN IB AND M-1B IN IT  THEN GMMATS INTO IV AND IB
 
 DO  i = 1,nn
   z(it+i-1) = -z(ibh+i-1) + rbv*z(iqi+i)
   z(ib+i-1) = -z(ikh+i-1) + rvs*z(iqr+i)
 END DO
 CALL gmmats (z(ib),nrow,nrow,0,z(imh),nrow,nrow,0,z(iv))
 CALL gmmats (z(it),nrow,nrow,0,z(imh),nrow,nrow,0,z(ib))
 
!     CALL FA1PKA TO MAKE A MATRIX AND GET EIGENVALUES
 
 CALL fa1pka (z(ima),z(iv),z(ib),z(it),ieigns-it,nrow)
 
!     SORT EIGENVALUES
 
 j = neign*2
 CALL rsort (2,1,z(it),j)
 CALL rsort (2,2,z(it),j)
 IF (kint /= 0.0) GO TO 180
 nlft = neign
 DO  i = 1,j,2
   IF (z(it+i) >= 0.0) EXIT
   nlft = nlft - 1
 END DO
 170 nl = it + (neign-nlft)*2
 nr = 0
 DO  i = 1,j,2
   IF (z(it+i) /= 0.0) CYCLE
   nr = nr + 1
   IF (eigv) CALL fa1pkv (z(ima),z(iv),z(ib),nrow,z(it+i-1),z(ima),  &
       bref,pi,vel,z(buf1))
 END DO
 nrs = nr + 1
 nr  = nr/2
 nra = 0
 180 CONTINUE
 IF (l39 == 0) GO TO 200
 WRITE  (nout,182) kint
 182 FORMAT (1H0,29H estimated reduced frequency ,1P,e15.5, /10X,  &
     11HEIGENVALUES,10X,18H reduced frequency,4X,9HFREQUENCY,  &
     6X,8H damping,/,7X,4HREAL,10X,4HIMAG)
 DO  i = 1,j,2
   er = z(it+i-1)
   ei = z(it+i  )
   IF (ei == 0.0) GO TO 183
   rk = bov*ei
   rf = (1.0/twopi)*ei
   rg = (2.0*er)/ei
   GO TO 185
   183 rk = 0.0
   rf = 0.0
   rg = (bref/(pi*vel))*er
   185 WRITE  (nout,186) er,ei,rk,rf,rg
   186 FORMAT (1H ,1P,e15.5,1P,e15.5,3X,1P,e15.5,1P,e15.5,1P,e15.5)
 END DO
 
!     ROOT ACCEPTANCE AND SAVING
 
 200 j = nlft*2
 l = nroot*2 + 1 + nra*2
 imhere = 200
 IF (l > j) GO TO 360
 
 DO  i = l,j,2
   k = (nroot*5) + 1 + ieigns
   IF (z(nl+i) /= 0.0) GO TO 220
   IF (kint    /= 0.0) GO TO 220
   IF (nrs /= nr) nrs = nrs - 1
   IF (nrs /= nr) CYCLE
   nra = nra + 1
   z(k  ) = z(nl+i-1)
   z(k+1) = z(nl+i  )
   z(k+2) = 0.0
   z(k+3) = 0.0
   z(k+4) = (bref/(.34657*vel))*z(nl+i-1)
   210 nroot  = nroot + 1
   
!     PRINT EIGENVECTORS IF ASKED FOR
   
   nit = 0
   
!     NO. OF ITERATIONS RESET TO ZERO. RE-INITIALIZE LEASE SQUARE COEFF.
   
   xav = 0.
   yav1= 0.
   x10 = 0.
   x11 = 0.
   x12 = 0.
   y10 = 0.
   y11 = 0.
   IF (nroot >= neiw) GO TO 300
   CYCLE
   220 rktst = bov*z(nl+i)
   IF (ABS(rktst-kint) < eps) GO TO 230
   IF (rktst == 0.0) GO TO 230
   
!     SUM LEAST SQUARE COEFFICIENTS ASSOCIATED WITH DEPENDENT VARIABLE
!     STARTING WITH RESULT OF SECOND TIRAL
   
   IF (nit == 1) GO TO 225
   yav1 = yav1 + rktst
   y10  = y10  + rktst
   y11  = y11  + rktst*kint
   225 kint = rktst
   IF (nit == 10) GO TO 240
   GO TO 110
   
!     START LOOP OVER
   
   230 z(k  ) = z(nl+i-1)
   z(k+1) = z(nl+i  )
   z(k+2) = rktst
   z(k+3) = (1.0/twopi)*z(nl+i)
   IF (z(nl+i) /= 0.0) z(k+4) = (2.0*z(nl+i-1))/z(nl+i)
   IF (z(nl+i) == 0.0) z(k+4) = (bref/(.34657*vel))*z(nl+i-1)
   IF (eigv) CALL fa1pkv (z(ima),z(iv),z(ib),nrow,z(k),z(ima),  &
       bref,pi,vel,z(buf1))
   GO TO 210
   
!     FAILURE TO CONVERGE. REPLACE LOOP END WITH LEAST SQUARES FIT
   
   240 nit  = nit + 1
   xav1 = xav/(nit-2)
   xav  = (xav + rktst)/(nit-1)
   yav1 = yav1/(nit-2)
   d1   = x12*x10  - x11*x11
   a11  = (x10*y11 - x11*y10)/d1
   a10  = (x12*y10 - x11*y11)/d1
   rktst= -a10/(a11-1.)
   WRITE  (nout,250) uwm,nit,floop,nroot,neiw
   250 FORMAT (a25,', PK METHOD FIALED TO CONVERGE', /1X,i4,  &
       ' ITERATIONS ON LOOP',i5,',  FOUND',i5,',  ROOTS WANTED',  &
       i5, /5X,'LEAST SQUARES FIT APPROXIMATION IMPLEMENTED.')
   IF (l39 == 1) WRITE (nout,260) xav1,yav1,xav, a11,a10,rktst
   260 FORMAT (/5X,'AVG. TRIAL = ',1P,e12.5,',  AGV. RESLT. = ',1P,e12.5,  &
       ',  NET AVG. = ',1P,e12.5,  //9X,'SLOPE = ',1P,e12.5,  &
       ',    INTERCEPT = ',1P,e12.5,',  VALUE    = ',1P,e12.5)
   GO TO 230
   
 END DO
 
!     LOGIC ERROR
 
 imhere = 270
 GO TO 360
 
!     SAVE EIGENVALUES ON BXHH
 
 300 IF (istart /= 0) GO TO 310
 istart = 1
 CALL gopen (bxhh,z(buf1),1)
 CALL CLOSE (bxhh,2)
 310 CALL gopen (bxhh,z(buf1),3)
 CALL WRITE (bxhh,z(ieigns+1),nroot*5,1)
 IF (floop >= nloop) GO TO 320
 CALL CLOSE (bxhh,3)
 RETURN
 
!     LAST LOOP BUILD FSAVE
 
 320 CALL CLOSE (bxhh,1)
 ibuf2 = buf1 - sysbuf
 CALL gopen (bxhh,z(buf1),0)
 CALL gopen (fsave,z(ibuf2),0)
 CALL skprec (fsave,3)
 CALL CLOSE (fsave,2)
 CALL gopen (fsave,z(ibuf2),3)
 330 CALL READ  (*350,*340,bxhh,z(1),ibuf2,1,nwr)
 340 CALL WRITE (fsave,z(1),nwr,1)
 GO TO 330
 350 CALL CLOSE (bxhh,1)
 CALL CLOSE (fsave,1)
 trl(1) = fsave
 trl(2) = nloop
 trl(7) = neiw
 CALL wrttrl (trl)
 GO TO 400
 
 360 WRITE  (nout,370) sfm,imhere,l,j
 370 FORMAT (a25,'. ERROR IN FA1PKE/@',i3,'  L,J=',2I7)
 CALL mesage (-61,0,0)
 
 400 RETURN
END SUBROUTINE fa1pke
