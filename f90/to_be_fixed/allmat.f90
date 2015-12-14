SUBROUTINE allmat (a,lambda,h,hl,vect,mult,inth,INT,m,ncal,iopt1)
     
!     SUBROUTINE ALLMAT (A,LAMBDA,M,IA,NCAL)
 
!     A ON ENTRY      = MATRIX TO BE ITERATED
!     A    ON RETURN  = EIGENVECTORS  (OPTIONAL)
!     LAMBDA          = EIGENVALUES
!     M               = FIRST DIMENSION OF (A)
!     NCAL ON ENTRY   = FLAG .NE. 0. COMPUTE VECTORS
!                            .EQ. 0. NO VECTORS
!     NCAL ON RETURN  = NUMBER OF EIGENVALUES
 
 
!     PROG. AUTHORS JOHN RINZEL AND R.E.FUNDERLIC, UNION CARBIDE CORP.
!     NUCLEAR DIVISION,CENTRAL DATA PROCESSING FACILITY,
!     OAK RIDGE, TENNESSEE
 
 
 COMPLEX, INTENT(IN OUT)                  :: a(m,m)
 COMPLEX, INTENT(OUT)                     :: lambda(2)
 COMPLEX, INTENT(OUT)                     :: h(m,m)
 COMPLEX, INTENT(OUT)                     :: hl(m,m)
 COMPLEX, INTENT(OUT)                     :: vect(1)
 COMPLEX, INTENT(OUT)                     :: mult(1)
 LOGICAL, INTENT(OUT)                     :: inth(1)
 INTEGER, INTENT(OUT)                     :: INT(1)
 INTEGER, INTENT(IN)                      :: m
 INTEGER, INTENT(IN OUT)                  :: ncal
 INTEGER, INTENT(IN OUT)                  :: iopt1
 LOGICAL :: twice
 INTEGER :: r,rp1,rp2
 COMPLEX :: shift(3),temp,SIN,COS,temp1,temp2
 
 nvec = ncal
 n    = m
 ncal = n
 IF (n /= 1) GO TO 1
 lambda(1) = a(1,1)
 a(1,1) = 1.
 GO TO 62
 1 icount = 0
 shift(1) = 0.
 IF (n /= 2) GO TO 4
 2 temp = (a(1,1)+a(2,2) + CSQRT((a(1,1)+a(2,2))*(a(1,1)+a(2,2)) -  &
     4.*(a(2,2)*a(1,1)-a(2,1)*a(1,2))))/2.
 IF (REAL(temp) /= 0. .OR. AIMAG(temp) /= 0.) GO TO 3
 lambda(m  ) = shift(1)
 lambda(m-1) = a(1,1) + a(2,2) + shift(1)
 GO TO 37
 3 lambda(m  ) = temp + shift(1)
 lambda(m-1) = (a(2,2)*a(1,1)-a(2,1)*a(1,2))/(lambda(m)-shift(1)) + shift(1)
 GO TO 37
 
!     REDUCE MATRIX A TO HESSENBERG FORM
 
 4 nm2 = n - 2
 DO  r = 1,nm2
   rp1 = r + 1
   rp2 = r + 2
   abig= 0.
   INT(r) = rp1
   DO  i = rp1,n
     abssq = REAL(a(i,r))**2 + AIMAG(a(i,r))**2
     IF (abssq <= abig) CYCLE
     INT(r) = i
     abig   = abssq
   END DO
   IF (abig == 0.) CYCLE
   inter = INT(r)
   IF (inter == rp1) GO TO 8
   DO  i = r,n
     temp = a(rp1,i)
     a(rp1,i) = a(inter,i)
     a(inter,i) = temp
   END DO
   DO  i = 1,n
     temp = a(i,rp1)
     a(i,rp1) = a(i,inter)
     a(i,inter) = temp
   END DO
   8 DO  i  = rp2,n
     mult(i) = a(i,r)/a(rp1,r)
     a(i,r)  = mult(i)
   END DO
   DO  i = 1,rp1
     temp = 0.
     DO  j = rp2,n
       temp = temp + a(i,j)*mult(j)
     END DO
     a(i,rp1) = a(i,rp1) + temp
   END DO
   DO  i = rp2,n
     temp = 0.
     DO  j = rp2,n
       temp = temp + a(i,j)*mult(j)
     END DO
     a(i,rp1) = a(i,rp1) + temp - mult(i)*a(rp1,rp1)
   END DO
   DO  i = rp2,n
     DO  j = rp2,n
       a(i,j) = a(i,j) - mult(i)*a(rp1,j)
     END DO
   END DO
 END DO
 
!     CALCULATE EPSILON
 
 eps = 0.
 DO  i = 1,n
   eps = eps + cabs(a(1,i))
 END DO
 DO  i = 2,n
   sum = 0.
   im1 = i - 1
   DO  j = im1,n
     sum = sum + cabs(a(i,j))
   END DO
   IF (sum > eps) eps = sum
 END DO
 eps = SQRT(FLOAT(n))*eps*1.e-12
 IF (eps == 0.) eps = 1.e-12
 DO  i = 1,n
   DO  j = 1,n
     h(i,j) = a(i,j)
   END DO
 END DO
 20 IF (n /= 1) GO TO 21
 lambda(m) = a(1,1) + shift(1)
 GO TO 37
 21 IF (n == 2) GO TO 2
 22 mn1 = m - n + 1
 ard = REAL (a(n,n))
 aid = AIMAG(a(n,n))
 arn = REAL (a(n,n-1))
 ain = AIMAG(a(n,n-1))
 IF (ard == 0.0 .AND. aid == 0.0) GO TO 23
 term1 = ABS(ard*arn + aid*ain)
 term2 = ABS(ard*ain - aid*arn)
 term3 = ard*ard + aid*aid
 IF ((term1+term2) <= 1.0E-9*term3) GO TO 24
 23 IF ((ABS(arn)+ABS(ain)) >= eps) GO TO 25
 24 lambda(mn1) = a(n,n) + shift(1)
 icount = 0
 n = n - 1
 GO TO 21
 
!     DETERMINE SHIFT
 
 25 shift(2) = (a(n-1,n-1)+a(n,n) + CSQRT((a(n-1,n-1)+a(n,n))*  &
     (a(n-1,n-1)+a(n,n)) - 4.*(a(n,n)*a(n-1,n-1)-a(n,n-1)* a(n-1,n))))/2.
 IF (REAL(shift(2)) /= 0. .OR. AIMAG(shift(2)) /= 0.) GO TO 26
 shift(3) = a(n-1,n-1) + a(n,n)
 GO TO 27
 26 shift(3) = (a(n,n)*a(n-1,n-1) - a(n,n-1)*a(n-1,n))/shift(2)
 27 IF (cabs(shift(2)-a(n,n)) < cabs(shift(3)-a(n,n))) GO TO 28
 INDEX = 3
 GO TO 29
 28 INDEX = 2
 29 IF (cabs(a(n-1,n-2)) >= eps) GO TO 30
 lambda(mn1  ) = shift(2) + shift(1)
 lambda(mn1+1) = shift(3) + shift(1)
 icount = 0
 n = n - 2
 GO TO 20
 30 shift(1) = shift(1) + shift(INDEX)
 DO  i = 1,n
   a(i,i) = a(i,i) - shift(INDEX)
 END DO
 
!     PERFORM GIVENS ROTATIONS, QR ITERATES
 
 IF (icount <= 20) GO TO 32
 ncal  = m - n
 GO TO 37
 32 nm1   = n - 1
 temp1 = a(1,1)
 temp2 = a(2,1)
 DO  r = 1,nm1
   rp1   = r + 1
   rho   = SQRT(REAL(temp1)**2 + AIMAG(temp1)**2 +  &
       REAL(:: temp2)**2 + AIMAG(temp2)**2)
   IF (rho == 0.) CYCLE
   COS   = temp1/rho
   SIN   = temp2/rho
   INDEX = MAX0(r-1,1)
   DO  i = INDEX,n
     temp  = CONJG(COS)*a(r,i) + CONJG(SIN)*a(rp1,i)
     a(rp1,i) =-SIN*a(r,i) + COS*a(rp1,i)
     a(r,i) = temp
   END DO
   temp1  = a(rp1,rp1)
   temp2  = a(r+2,r+1)
   DO  i = 1,r
     temp   = COS*a(i,r) + SIN*a(i,rp1)
     a(i,rp1) =-CONJG(SIN)*a(i,r) + CONJG(COS)*a(i,rp1)
     a(i,r) = temp
   END DO
   INDEX  = MIN0(r+2,n)
   DO  i = rp1,INDEX
     a(i,r) = SIN*a(i,rp1)
     a(i,rp1) = CONJG(COS)*a(i,rp1)
   END DO
 END DO
 icount = icount + 1
 GO TO 22
 
!     CALCULATE VECTORS
 
 37 IF (ncal == 0 .OR. nvec == 0) GO TO 62
 n   = m
 nm1 = n - 1
 IF (n /= 2) GO TO 38
 eps = AMAX1(cabs(lambda(1)),cabs(lambda(2)))*1.e-8
 IF (eps == 0.) eps = 1.e-12
 h(1,1) = a(1,1)
 h(1,2) = a(1,2)
 h(2,1) = a(2,1)
 h(2,2) = a(2,2)
 38 DO  l = 1,ncal
   DO  i = 1,n
     DO  j = 1,n
       hl(i,j) = h(i,j)
     END DO
     hl(i,i) = hl(i,i) - lambda(l)
   END DO
   DO  i = 1,nm1
     mult(i) = 0.
     inth(i) = .false.
     ip1 = i + 1
     IF (cabs(hl(i+1,i)) <= cabs(hl(i,i))) GO TO 42
     inth(i) = .true.
     DO  j = i,n
       temp = hl(i+1,j)
       hl(i+1,j) = hl(i,j)
       hl(i,j  ) = temp
     END DO
     42 IF (REAL(hl(i,i)) == 0. .AND. AIMAG(hl(i,i)) == 0.) CYCLE
     mult(i) = -hl(i+1,i)/hl(i,i)
     DO  j = ip1,n
       hl(i+1,j) = hl(i+1,j) + mult(i)*hl(i,j)
     END DO
   END DO
   DO  i = 1,n
     vect(i) = 1.
   END DO
   twice = .false.
   46 IF (REAL(hl(n,n)) == 0. .AND. AIMAG(hl(n,n)) == 0.) hl(n,n) = eps
   vect(n) = vect(n)/hl(n,n)
   DO  i = 1,nm1
     k = n - i
     DO  j = k,nm1
       vect(k) = vect(k) - hl(k,j+1)*vect(j+1)
     END DO
     IF (REAL(hl(k,k)) == 0. .AND. AIMAG(hl(k,k)) == 0.) hl(k,k) = eps
     vect(k) = vect(k)/hl(k,k)
   END DO
   big = 0.
   DO  i = 1,n
     sum = ABS(REAL(vect(i))) + ABS(AIMAG(vect(i)))
     IF (sum > big) big = sum
   END DO
   DO  i = 1,n
     vect(i) = vect(i)/big
   END DO
   IF (twice) GO TO 52
   DO  i = 1,nm1
     IF (.NOT.inth(i)) GO TO 51
     temp = vect(i)
     vect(i  ) = vect(i+1)
     vect(i+1) = temp
     51 vect(i+1) = vect(i+1) + mult(i)*vect(i)
   END DO
   twice = .true.
   GO TO 46
   52 IF (n == 2) GO TO 55
   nm2 = n - 2
   DO  i = 1,nm2
     n1i = n - 1 - i
     ni1 = n - i + 1
     DO  j = ni1,n
       vect(j) = h(j,n1i)*vect(n1i+1) + vect(j)
     END DO
     INDEX = INT(n1i)
     temp  = vect(n1i+1)
     vect(n1i+1) = vect(INDEX)
     vect(INDEX) = temp
   END DO
   55 DO  i = 1,n
     a(i,l) = vect(i)
   END DO
 END DO
 DO  j = 1,ncal
   te = 0.
   DO  i = 1,n
     tem = cabs(a(i,j))
     IF (te > tem) CYCLE
     l = i
     te = tem
   END DO
   temp1 = a(l,j)
   DO  i = 1,n
     a(i,j) = a(i,j)/temp1
   END DO
 END DO
 62 RETURN
END SUBROUTINE allmat
