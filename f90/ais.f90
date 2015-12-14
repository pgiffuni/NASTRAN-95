REAL FUNCTION ais (np,i,l,rr,zz)
     
!     THIS ROUTINE CALCULATES THE SINGLE PRECISION  DELTA(IJ) INTEGRALS
!     FOR AXISYMMETRIC SOLIDS IN SDR2. CALCULATIONS DONE IN DOUBLE
!     PRECISION
 
!     INPUT
!     NP  = NUMBER OF POINTS (3 OR 4.  MORE THAN 4 WILL FAIL FSN-5)
!     I,L = THE INTEGRAL DESIRED (I SERIES STARTS WITH -1)
!     R   = RADIUS ARRAY (NP LONG)
!     Z   = Z-CORD ARRAY (NP LONG)
 
!     OUTPUT
!     DKL = DESIRED INTEGRAL
 
 
 INTEGER, INTENT(IN)                      :: np
 INTEGER, INTENT(IN)                      :: i
 INTEGER, INTENT(IN)                      :: l
 REAL, INTENT(IN)                         :: rr(4)
 REAL, INTENT(IN)                         :: zz(4)
 INTEGER :: nam(2)
 DOUBLE PRECISION :: a,aj,ar,beta,dr,dz,dfact,dzj,dl1,dkl,eps,fact,  &
     gkl,zero,one,two,three,pr,ra,rb,r(4),rak,rbk, za,zb,z(4)
 
 DATA       eps / .01D0 /
 DATA       zero, one,two,three / 0.0D0, 1.0D0, 2.0D0, 3.0D0 /
 DATA       nam / 4HAIS , 1H    /
 
 DO  m = 1,np
   r(m) = rr(m)
   z(m) = zz(m)
 END DO
 dkl  = zero
 l1   = l + 1
 l2   = l + 2
 dl1  = l1
 k    = i + 1
 
!     LOOP ON NUMBER OF POINTS
 
 IF (r(1) <= zero) GO TO 300
 DO  m = 1,np
   j  = m + 1
   IF (m == np) j = 1
   ra = r(m)
   rb = r(j)
   za = z(m)
   zb = z(j)
   dr = rb - ra
   dz = zb - za
   
!     TEST IF RADIUS IS .LE. 0 (DRIVER SHOULD FIND THIS)
   
   IF (rb <= zero) GO TO 300
   gkl = zero
   pr  = ra + rb
   ar  = pr/two
   
!     CHECK FOR APPROXIMATION, DR/AVE(R)
   
   IF ( DABS ( dr/ar ) < eps ) GO TO 70
   
   a = za*dr - ra*dz
   beta = a/dr
   
!     CHECK FOR BETA .EQ. 0 CASE
   
   IF (DABS(beta/ar) > eps) GO TO 10
   
   IF (dz == zero) GO TO 200
   lk  = l + k + 1
   ar  = lk
   gkl = (dz/dr)**l1*(ra**lk-rb**lk)/(dl1*ar)
   GO TO 200
   
!     GENERAL CASE
   
   10 rak = ra**k
   rbk = rb**k
   IF (k < 0) THEN
     GO TO   300
   ELSE IF (k == 0) THEN
     GO TO    20
   ELSE
     GO TO    30
   END IF
   
!     GENERAL CASE, K.EQ.0, CONSTANT TERM
   
   20 gkl = DLOG(ra/rb)/dl1
   GO TO 40
   
!     GENERAL CASE, CONSTANT TERM
   
   30 ar  = k*l1
   gkl = (rak-rbk)/ar
   
!     GENERAL CASE, SUMMATION
   
   40 IF (dz == zero) GO TO 65
   lfact = 1
   
!     CALCULATE FACTORIAL (L+1)
   
   DO  j = 2,l
     lfact = lfact*j
   END DO
   factl = lfact
   jfact = 1
   aj    = one
   dzj   = one
   lmjf  = lfact*l1
   DO  j = 1,l1
     jfact = jfact*j
     
!     CALCULATE (L+1-J) FACTORIAL IN LMJF
     
     lmjf  = lmjf /(l2-j)
     fact  = factl/FLOAT (jfact*lmjf)
     dfact = k + j
     dfact = fact /dfact
     aj    = aj *a
     rak   = rak*ra
     rbk   = rbk*rb
     dzj   = dzj*dz
     gkl   = gkl + (dfact*dzj*(rak-rbk))/aj
   END DO
   
   65 gkl   = gkl*beta**l1
   GO TO 200
   
!     APPROXIMATE CODE
   
   70 CONTINUE
   IF (dr == zero) GO TO 200
   dzj = l1*l2
   rbk = zb**l1
   j   = k - 1
   gkl =-dr*ar**j*rbk/dl1
   
   IF (dz == zero) GO TO 200
   gkl = gkl + (((two*ra+rb)/three)**j*dr*DABS(za**l2-rbk*zb)) /(dzj*dz)
   
   200 dkl = dkl + gkl
 END DO
 
!     ALL DONE
 
 210 ais = dkl
 RETURN
 
!     ERROR
 
 300 CALL mesage (-7,k,nam)
 GO TO 210
END FUNCTION ais
