SUBROUTINE wilvec (d,o,val,vloc,v,f,p,q,r,vec,nx,svec)
     
!     WILKINSON EIGENVECTOR SOLUTION FOR LARGE SYM MATRICES
 
 
 DOUBLE PRECISION, INTENT(IN)             :: d(1)
 DOUBLE PRECISION, INTENT(IN)             :: o(1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: val(1)
 INTEGER, INTENT(IN OUT)                  :: vloc(1)
 DOUBLE PRECISION, INTENT(OUT)            :: v(1)
 DOUBLE PRECISION, INTENT(OUT)            :: f(1)
 DOUBLE PRECISION, INTENT(OUT)            :: p(1)
 DOUBLE PRECISION, INTENT(OUT)            :: q(1)
 DOUBLE PRECISION, INTENT(OUT)            :: r(1)
 DOUBLE PRECISION, INTENT(OUT)            :: vec(nx,1)
 INTEGER, INTENT(IN OUT)                  :: nx
 INTEGER, INTENT(IN OUT)                  :: svec(1)
 INTEGER :: ENTRY,v2,xentry,pv,vector,vv,v1,mcb(7), sysbuf,mcb1(7),phia, path
 DOUBLE PRECISION :: value,w,x,y,z,dlmdas,  &
     rmult,rrmult,sft,sftinv,deps,vmult,zero,one
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm,uwm,uim
 COMMON /givn  /  title(1),mo,t3,mr,mt1,t6,mv1,t8(3),ENTRY,t12(5),  &
     rstrt,v2,t19,xentry,t21(80),vcom(30),t131(20)
 COMMON /packx /  it1,it2,ii,jj,incr
 COMMON /unpakx/  it3,iii,jjj,incr1
 COMMON /system/  sysbuf,iotpe,ksys(52),iprec
 COMMON /reigkr/  ioptn
 COMMON /mgivxx/  dlmdas
 EQUIVALENCE      (n     ,vcom( 1)), (pv    ,vcom( 5)),  &
     (nv    ,vcom( 7)), (nrigid,vcom(10)),  &
     (phia  ,vcom(12)), (nver  ,vcom(13)), (maxitr,vcom(15)), (iterm ,vcom(16))
 DATA    mul3  ,  mcb1,mcb / 0, 0,0,0,2,2,0,0, 7*0    /
 DATA    zero  ,  one      / 0.0D+0,  1.0D+0          /
 DATA    mgiv  /  4HMGIV   /
 
!     D        = DIAGONAL     TERMS OF THE TRIDIAGONAL MATRIX (N)
!     O        = OFF-DIAGONAL TERMS OF THE TRIDIAGONAL MATRIX (N)
!     VAL      = EIGENVALUES (NV)
!     VLOC     = ORIGINAL ORDERING OF THE EIGENVALUES (NV)
!     V,F,P,Q,R= N DIMENSIONAL ARRAYS
!     VEC      = THE REST OF OPEN CORE
 
!     MT       = TRANSFORMATION TAPE
!     N        = ORDER  OF PROBLEM
!     NV       = NUMBER OF EIGENVECTORS
!     RSTRT
!     V2       = NUMBER   OF EIGENVECTORS ALREADY CLLCULATED
!     VV       = POINTER  TO CURRENT VECTOR IN CORE VEC(1,VV)
!     NM2X     = MIDPOINT OF PROBLEM (SWITCH SINE SAVE TAPES)
 
 
!     INITALIZE VARIABLES
 
 deps  = 1.0D-35
 sft   = 1.0D+20
 sftinv= 1.0D+0/sft
 vmult = 1.0D-02
 nz    = korsz(svec)
 ibuf1 = nz - sysbuf + 1
 ibuf2 = ibuf1 - sysbuf
 im    = 1
 CALL makmcb (mcb1,phia,n,2,iprec)
 im1   = 2
 nz    = ibuf2 - 1
 path  = 0
 nv1   = (nz-1)/(n+n)
 im2   = 2
 nm1   = n - 1
 nm2   = n - 2
 nver  = 0
 v2    = nrigid
 
!     REARRANGE EIGENVALUES AND EXTRACTION ORDER FOR MULTIPLE ROOTS
!     TO GUARANTEE THAT THEY ARE IN SUBMATRIX ORDER FOR PURPOSES
!     OF TRIAL VECTOR AND ORTHOGONOLIZATION COMPUTATIONS
 
 rmult = vmult
 rrmult= vmult/100.0D0
 iclos = 0
 i = nrigid + 1
 10 IF (DABS(val(i))+DABS(val(i+1)) < rmult) GO TO 20
 IF (val(i) == zero) GO TO 90
 IF (DABS(one-val(i)/val(i+1)) > rrmult) GO TO 80
 20 IF (iclos /= 0) GO TO 90
 iclos = i
 GO TO 90
 30 CONTINUE
 DO  i1 = iclos,i
   MIN   = vloc(i1)
   value = val(i1)
   k = i1
   DO  j = i1,i
     IF (vloc(j) >= MIN) CYCLE
     k    = j
     MIN  = vloc(j)
     value= val(j)
   END DO
   vloc(k) = vloc(i1)
   val(k)  = val(i1)
   vloc(i1)= MIN
   val(i1) = value
 END DO
 iclos = 0
 80 IF (iclos /= 0) GO TO 30
 90 i = i + 1
 IF (i     < nv) GO TO 10
 IF (iclos /=  0) GO TO 30
 
!     START LOOP FOR CORE LOADS OF VECTORS
 
 100 CALL klock (ist)
 v1   = v2 + 1
 v2   = v2 + nv1
 mul2 = mul3
 mulp2= 0
 mul3 = 0
 IF (nv-v2 < 0) THEN
   GO TO   101
 ELSE IF (nv-v2 == 0) THEN
   GO TO   110
 ELSE
   GO TO   102
 END IF
 101 v2   = nv
 GO TO 110
 
!     SEARCH FOR MULTIPLICITIES OF EIGENVALUES V2 AND V2+1.
 
 102 vv = v2
 103 IF (DABS(val(v2))+DABS(val(v2+1)) < rmult) GO TO 1041
 IF (DABS(one-val(v2)/val(v2+1))   > rmult) GO TO 110
 1041 CONTINUE
 l1 = vloc(v2  )
 l2 = vloc(v2+1)
 n1 = MIN0(l1,l2)
 n2 = MAX0(l1,l2) - 1
 DO  k = n1,n2
   IF (o(k) == zero) GO TO 110
 END DO
 v2 = v2 - 1
 IF (v2+6 > n1 .AND. v2 > v1) GO TO 103
 v2 = vv
 mul3  = 1
 
!     FIND EIGENVECTORS V1 - V2.
 
 110 n1 = 0
 n2 = 0
 nv2= v2 - v1 + 1
 DO  vv = 1,nv2
   vector = v1 + vv - 1
   value  = val (vector)
   
!     FOR MGIV METHOD, USE ORIGINAL LAMBDA COMPUTED BY QRITER
!     IN EIGENVECTOR COMPUTATIONS
   
   IF (ioptn == mgiv) value = 1.0D0/(value + dlmdas)
   loc = vloc(vector)
   IF (loc >= n1 .AND. loc <= n2) GO TO 120
   
!     SEARCH FOR A DECOUPLED SUBMATRIX.
   
   mul1 = 0
   IF (loc == 1) GO TO 112
   DO  k = 2,loc
     n1 = loc - k + 2
     IF (o(n1-1) == zero) GO TO 113
   END DO
   112 n1 = 1
   113 IF (loc == n) GO TO 115
   DO  k = loc,nm1
     IF (o(k) == zero) GO TO 116
   END DO
   115 n2 = n
   GO TO 120
   116 n2 = k
   120 IF (mul1 /= 0 .OR. mul2 /= 0) GO TO 122
   DO  i = 1,n
     v(i) = zero
   END DO
   IF (n1 /= n2) GO TO 122
   v(loc) = one
   GO TO 152
   122 n2m1 = n2 - 1
   n2m2 = n2 - 2
   
!     SET UP SIMULTANEOUS EQUATIONS
   
   x = d(n1) - value
   y = o(n1)
   DO  k = n1,n2m1
     IF (x == zero) GO TO 125
     f(k) = -o(k)/x
     GO TO 126
     125 f(k) = -sft*o(k)
     126 IF (DABS(x)-DABS(o(k)) < 0.0) THEN
       GO TO   127
     ELSE IF (DABS(x)-DABS(o(k)) == 0.0) THEN
       GO TO   128
     ELSE
       GO TO   129
     END IF
     
!     PIVOT.
     
     127 p(k) = o(k)
     q(k) = d(k+1) - value
     r(k) = o(k+1)
     z    =-x/p(k)
     x    = z*q(k) + y
     y    = z*r(k)
     GO TO 130
     
!     DO NOT PIVOT.
     
     128 IF (x == zero) x = sftinv
     129 p(k) = x
     q(k) = y
     r(k) = zero
     x    = d(k+1) - (value+o(k)*(y/x))
     y    = o(k+1)
     130 CONTINUE
   END DO
   IF (mul1 /= 0 .OR. mul2 /= 0) GO TO 135
   DO  k = n1,n2m1
     v(k) = one
   END DO
   w    = one/DSQRT(DBLE(FLOAT(n2-n1+1)))
   v(n2)= one
   
!     SOLVE FOR AN EIGENVECTOR OF THE TRIDIAGONAL MATRIX.
   
   135 mul2   = 0
   maxitr = 3
   DO  iter = 1,maxitr
     
!     BACK SUBSTITUTION
     
     IF (x == zero) GO TO 136
     v(n2) = v(n2)/x
     GO TO 137
     136 v(n2  ) = v(n2)*sft
     137 v(n2-1) = (v(n2-1) - q(n2-1)*v(n2))/p(n2-1)
     MAX = n2
     IF (DABS(v(n2)) < DABS(v(n2-1))) MAX = n2m1
     IF (n2m2 < n1) GO TO 140
     DO  k = n1,n2m2
       l    = n2m2 - (k-n1)
       v(l) = (v(l)-q(l)*v(l+1) - r(l)*v(l+2))/p(l)
       IF (DABS(v(l)) > DABS (v(MAX))) MAX = l
     END DO
     
!     NORMALIZE THE VECTOR.
     
     140 y = DABS(v(MAX))
     z = zero
     DO  i = n1,n2
       v(i) = v(i)/y
       IF (DABS(v(i)) < deps) CYCLE
       z = z + v(i)*v(i)
     END DO
     z = DSQRT(z)
     DO  i = n1,n2
       v(i) = v(i)/z
     END DO
     
!     CHECK CONVERGENCE OF THE LARGEST COMPONENT OF THE VECTOR.
     
     y = DABS(v(MAX))
     IF (SNGL(w) == SNGL(y)) GO TO 152
     IF (iter    ==  maxitr) CYCLE
     w = y
     
!     PIVOT V.
     
     DO  i = n1,n2m1
       IF (p(i) == o(i)) GO TO 144
       v(i+1) = v(i+1) + v(i)*f(i)
       CYCLE
       144 z =  v(i+1)
       v(i+1) = v(i) + z/f(i)
       v(i  ) = z
     END DO
   END DO
   
!     TOO MANY ITERATIONS.
   
!     THE ACCURACY OF EIGENVECTOR XXXX CORRESPONDING TO THE EIGENVALUE
!     XXXXXXX  IS IN DOUBT.
   
   152 DO  i = 1,n
     vec(i,vv) = v(i)
   END DO
   
!     CHECK MULTIPLICITY OF THE NEXT EIGENVALUE IF IT IS IN THE SAME
!     SUBMATRIX AS THIS ONE.
   
   IF (vector == v2) GO TO 160
   
!     FOR MGIV METHOD, USE ADJUSTED LAMBDA COMING OUT OF QRITER
!     IN THE FOLLOWING CHECKS
   
   IF (DABS(val(vector+1))+DABS(val(vector)) < rmult) GO TO 154
   IF (DABS(val(vector+1)-val(vector)) > rmult*DABS(val(vector+1))) GO TO 160
   154 CONTINUE
   l1 = vloc(vector+1)
   IF (l1 < n1 .OR. l1 > n2) GO TO 160
   
!     A MULTIPLICITY DOES EXIT...THE INITIAL APPROXIMATION OF THE NEXT
!     EIGENVECTOR SHOULD BE ORTHOGONAL TO THE ONE JUST CALCULATED.
   
   IF (mul1 == 0) mul1 = vv
   mulp2 = mulp2 + 1
   mulp3 = mulp2 + mul1 - 1
   DO  kkk = n1,n2
     v(kkk) = one
   END DO
   DO  jjj = mul1,mulp3
     z = zero
     DO  kk = n1,n2
       DO  ii = n1,n2
         z = z + vec(ii,jjj)*v(ii)
       END DO
       v(kk) = v(kk) - z*vec(kk,jjj)
     END DO
   END DO
   CYCLE
   
!     DOES THIS EIGENVALUE = PREVIOUS ONE(S) IN THIS SUBMATRIX
   
   160 IF (mul1 == 0) CYCLE
   
!     A MULTIPLICITY OF EIGENVALUES OCCURRED...IMPROVE THE ORTHOGONALITY
!     OF THE CORRESPONDING EIGENVECTORS.
   
   mulp1 = mul1 + 1
   DO  l = mulp1,vv
     DO  i = n1,n2
       p(i) = vec(i,l)
       q(i) = zero
     END DO
     lm1 = l - 1
     DO  k = mul1,lm1
       z = zero
       DO  i = n1,n2
         z = z + p(i)*vec(i,k)
       END DO
       DO   i = n1,n2
         q(i) = q(i) + z*vec(i,k)
       END DO
     END DO
     z = zero
     DO   k = n1,n2
       q(k) = p(k) - q(k)
       IF (DABS(q(k)) < deps) CYCLE
       z = z + q(k)*q(k)
     END DO
     z = DSQRT(z)
     DO  k = n1,n2
       vec(k,l) = q(k)/z
     END DO
   END DO
   mul1  = 0
   mulp2 = 0
 END DO
 
!     CORE IS NOW FULL OF EIGENVECTORS OF THE TRIDIAGONAL MATRIX.
!     CONVERT THEM TO EIGENVECTORS OF THE ORIGINAL MATRIX.
 
 it1 = 2
 it2 = 2
 jj  = n
 incr= 1
 
!     IS THE ORIGINAL MATRIX A 2X2
 
 IF (nm2  == 0) GO TO 186
 mt = mt1
 IF (path /= 0) GO TO 176
 mt = mo
 176 CALL gopen (mt,svec(ibuf1),im2)
 IF (path == 0 .AND. v2 /= nv) CALL gopen (mt1,svec(ibuf2),1)
 it3  = 2
 jjj  = n
 incr1= 1
 DO  m = 1,nm2
   l1  = n - m
   iii = l1+ 1
   IF (path == 0) CALL bckrec (mt)
   CALL unpack (*167,mt,p)
   GO TO 180
   167 DO  i = 1,m
     p(i) = zero
   END DO
   180 IF (path /= 0 .OR. v2 == nv) GO TO 177
   ii = l1+1
   CALL pack (p,mt1,mcb)
   177 IF (path == 0) CALL bckrec (mt)
   DO  k = 1,m
     l2 = n - k + 1
     i  = m - k + 1
     y  = p(i)
     IF (y == zero) CYCLE
     x = zero
     IF (DABS(y) < one) x = DSQRT(one-y**2)
     DO  vv = 1,nv2
       z = x*vec(l1,vv) -y*vec(l2,vv)
       vec(l2,vv) = x*vec(l2,vv) + y*vec(l1,vv)
       vec(l1,vv) = z
     END DO
   END DO
 END DO
 CALL CLOSE (mt,1)
 IF (path /= 0) GO TO 186
 IF (v2 /= nv) WRITE (iotpe,1001) uim,n,nv,nv1
 1001 FORMAT (a29,' 2016A, WILVEC EIGENVECTOR COMPUTATIONS.', /37X,  &
     'PROBLEM SIZE IS',i6,', NUMBER OF EIGENVECTORS TO BE ',  &
     'RECOVERED IS',i6 , /37X,'SPILL WILL OCCUR FOR THIS ',  &
     'CORE AT RECOVERY OF',i6,' EIGENVECTORS.')
 path = 1
 CALL CLOSE (mt1,1)
 im2  = 0
 
!     WRITE THE EIGENVECTORS ONTO PHIA
 
 186 CALL gopen (phia,svec(ibuf1),im)
 ii  = 1
 it2 = iprec
 IF (im /= 1 .OR. nrigid <= 0) GO TO 205
 
!     PUT OUT ZERO VECTORS FOR RIGID BODY MODES
 
 jj = 1
 DO  vv = 1,nrigid
   CALL pack (zero,phia,mcb1)
 END DO
 jj = n
 205 CONTINUE
 im = 3
 IF (n == 1) GO TO 250
 DO  vv = 1,nv2
   CALL pack (vec(1,vv),phia,mcb1)
 END DO
 250 IF (v2 == nv) im1 = 1
 CALL CLOSE (phia,im1)
 xentry = -ENTRY
 
!     ANY TIME LEFT TO FIND MORE
 
 CALL tmtogo (itime)
 CALL klock  (ifin)
 IF (2*(ifin-ist) >= itime) GO TO 200
 IF (v2 /= nv) GO TO 100
 201 CALL wrttrl (mcb1)
 RETURN
 
!     MAX TIME
 
 200 iterm = 3
 GO TO 201
END SUBROUTINE wilvec
