SUBROUTINE qriter (val,o,loc,qr)
     
!     ORTEGA-KAISER QR ITERATION FOR A LARGE TRIDIAGONAL MATRIX
 
 
 DOUBLE PRECISION, INTENT(IN OUT)         :: val(1)
 DOUBLE PRECISION, INTENT(IN OUT)         :: o(1)
 INTEGER, INTENT(OUT)                     :: loc(1)
 INTEGER, INTENT(IN OUT)                  :: qr
 INTEGER :: sysbuf,msg(10)
 REAL :: lfreq
 DOUBLE PRECISION :: shift,zero,one,ones,epsi,g,r,s,t,u, dlmdas
 CHARACTER (LEN=5) :: below,above,belabv
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm,uwm,uim
 COMMON /system/  sysbuf,nout
 COMMON /givn  /  idum0(100),n,lfreq,idum3,idum4,hfreq,lama,nv,  &
     NE,idum9,nfound,idum11,idum12,idum13,never,MAX
 COMMON /reigkr/  ioptn
 COMMON /mgivxx/  dlmdas
 DATA    epsi  ,  zero,one,msg/ 1.0D-10, 0.0D+0, 1.0D+0, 53,9*0 /
 DATA    mgiv  ,  below,above / 4HMGIV, 'BELOW', 'ABOVE'        /
 
!     VAL    = DIAGONAL TERMS OF THE TRIDIAGONAL.
!              REORDERED EIGENVALUES UPON RETURN.
!     O      = SQUARE OF THE OFF-DIAGIONAL TERMS OF THE TRIDIAGONAL.
!     LOC    = ORIGINAL LOCATIONS OF THE REORDERED EIGENVALUES.
!     QR     = 1 MEANS  VAL= EIGENVALUES--JUST REORDER THEM
!     N      = ORDER OF THE PROBLEM = ALSO NO. OF FREQ. EXTRACTED
!     MAX    = MAXIMUM NUMBER OF ITERATIONS
!     SHIFT  = SHIFT FACTOR (SMALLEST DIAGONAL TERM
!     LFREQ  , HFREQ = FREQ. RANGE OF INTEREST IF NV IS ZERO
!     NV     = NUMBER OF EIGENVECTORS TO BE COMPUTED, SAVED AND OUTPUT.
!              IF NV IS ZERO (INPUT), AND LFREQ-HFREQ ARE PRESENT, NV IS
!              SET TO BE NO. OF MODES WITHIN THE FREQ. RANGE (OUTPUT)
!     NE     = NO. OF EIGENVALUES (INCLUDING RIGID MODES) TO BE PRINTED.
!              ALL, IF NE IS NOT SPECIFIED.
!              IF NE .LT. NV, NE IS SET EQUAL TO NV
 
 MAX = 100*n
 IF (nv >  n) nv = n
 IF (NE ==  0) NE = n
 IF (NE < nv) NE = nv
 
!     IS THIS AN ORDERING ONLY CALL
 
 never = 0
 IF (qr /= 0) GO TO 150
 
!     SEARCH FOR A DECOUPLED SUBMATRIX.
 
 m2 = n
 100 m2m1 = m2 - 1
 DO  k = 1,m2m1
   m1 = m2 - k
   IF (o(m1) /= zero) GO TO 102
 END DO
 
!     ALL OFF-DIAGONAL TERMS ARE ZEROS, JOB DONE. GO TO 150
!     THE DIAGONALS CONTAIN THE EIGENVALUES.
 
 GO TO 150
 
!     DECOUPLED SUBMATRIX
 
 102 m2m1 = m1
 m2   = m1 + 1
 IF (m2m1 == 1) GO TO 105
 DO  k = 2,m2m1
   m1 = m2 - k
   IF (o(m1) == zero) GO TO 104
 END DO
 GO TO 105
 104 m1 = m1 + 1
 105 mm = m1
 
!     Q-R ITERATION FOR THE DECOUPLED SUBMATRIX
 
 110 DO  iter = 1,MAX
   IF (DABS(val(m2))+o(m2m1) == DABS(val(m2))) GO TO 140
   DO  k = m1,m2m1
     IF (val(k) /= val(k+1)) GO TO 115
   END DO
   shift = zero
   GO TO 120
   
!     FIND THE SMALLEST DIAGONAL TERM = SHIFT
   
   115 shift = val(m2)
   DO  i = m1,m2m1
     IF (DABS(val(i)) < DABS(shift)) shift = val(i)
   END DO
   
!     REDUCE ALL TERMS BY SHIFT
   
   DO   i = m1,m2
     val(i) = val(i) - shift
   END DO
   
!     Q-R ITERATION
   
   120 r = val(m1)**2
   s = o(m1)/(r+o(m1))
   t = zero
   u = s*(val(m1) + val(m1+1))
   val(m1) = val(m1) + u
   IF (m1 == m2m1) GO TO 125
   m1p1 = m1 + 1
   DO  i = m1p1,m2m1
     g = val(i) - u
     r = (one-t)*o(i-1)
     ones = one - s
     IF (DABS(ones) > epsi) r = g*g/ones
     r = r + o(i)
     o(i-1) = s*r
     IF (o(i-1) == zero) mm = i
     t = s
     
!     IBM MAY FLAG AN EXPONENT UNDERFLOW ON NEXT LINE.
!     IT IS PERFECTLY OK SINCE O(I) SHOULD BE APPROACHING ZERO.
     
     s = o(i)/r
     u = s*(g + val(i+1))
     val(i) = u + g
   END DO
   
   125 val(m2) = val(m2) - u
   r = (one-t)*o(m2m1)
   ones = one - s
   IF (DABS(ones) > epsi) r = val(m2)**2/ones
   o(m2m1) = s*r
   
!     SHIFT BACK
   
   IF (shift == zero) GO TO 133
   DO  i = m1,m2
     val(i) = val(i) + shift
   END DO
   133 m1 = mm
 END DO
 
!     TOO MANY ITERATIONS
 
 
!     THE ACCURACY OF EIGENVALUE  XXXXX  IS IN DOUBT--QRITER FAILED TO
!     CONVERGE IN  XX  ITERATIONS
 
 never = never + 1
 CALL mesage (msg(1),val(m2),MAX)
 
!     CONVERGENCE ACHIEVED
 
 140 IF (m1 == m2m1) GO TO 145
 m2   = m2m1
 m2m1 = m2 -1
 GO TO 110
 145 IF (m1 <= 2) GO TO 150
 m2 = m1 - 1
 GO TO 100
 150 IF (n == 1) GO TO 205
 
!     REORDER EIGENVALUES ALGEBRAICALLY IN ASCENDING ORDER
 
 IF (ioptn /= mgiv) GO TO 155
 
!     FOR MGIV METHOD, RECOMPUTE LAMBDA
 
 DO  k = 1,n
   val(k) = (1.0D0/val(k)) - dlmdas
 END DO
 155 CONTINUE
 DO  k = 1,n
   DO  m = 1,n
     IF (val(m) /= -10000.0D0) EXIT
   END DO
   170 IF (m == n) GO TO 185
   mp1 = m + 1
   DO  i = mp1,n
     IF (val(i) == -10000.0D0) CYCLE
     IF (val(m) > val(i)) m = i
   END DO
   185 o(k)   = val(m)
   val(m) =-10000.0D0
   loc(k) = m
 END DO
 DO  i = 1,n
   val(i) = o(i)
 END DO
 
!     IF RIGID MODES WERE FOUND BEFORE, REPLACE RIGID FREQ. BY ZERO
 
 IF (nfound == 0) GO TO 205
 DO  i = 1,nfound
   val(i) = zero
 END DO
 
!     OUTPUT OPTION CHECK - BY FREQ. RANGE OR BY NO. OF FREQ.
!     REQUESTED
 
 205 ib    = 1
 IF (nv /= 0) GO TO 225
 IF (lfreq <= 0.0) GO TO 225
 
!     LOCATE PONTER THAT POINTS TO EIGENVALUE ABOVE OR EQUAL THE
!     LOWEST LFREQ. AS REQUESTED.
 
 DO  i = 1,n
   IF (val(i) >= lfreq) GO TO 220
 END DO
 i  = 0
 220 ib = i
 
!     OPEN LAMA FOR OUTPUT
!     PUT EIGENVALUES ON LAMA FOLLOWED BY ORDER FOUND
 
 225 ibuf1 = (korsz(o)-sysbuf+1)/2
 CALL gopen (lama,o(ibuf1),1)
 nn = 0
 IF (ib == 0) GO TO 240
 DO  i = ib,n
   valx = val(i)
   IF (nv /= 0 .AND.    i >   NE) GO TO 240
   IF (nv == 0 .AND. valx > hfreq) GO TO 240
   CALL WRITE (lama,valx,1,0)
   nn = nn + 1
 END DO
 
 240 CONTINUE
 
!     IF FREQ. RANGE IS REQUESTED, AND ALL FREQ. FOUND ARE OUTSIDE THE
!     RANGE, OUTPUT AT LEAST ONE FREQ.
 
 IF (nn > 0) GO TO 260
 IF (ib == 0) belabv = below
 IF (ib /= 0) belabv = above
 WRITE (nout,250) uim,belabv
 250 FORMAT (a29,', ALL ROOTS FOUND WERE ',a5,' FREQ. RANGE SPECIFIED',  &
     /5X,'HOWEVER, ONE EIGENVALUE OUTSIDE THIS FREQ. RANGE WAS',  &
     ' SAVED AND PRINTED')
 nn = 1
 IF (ib /= 0) ib = n
 IF (ib == 0) ib = 1
 CALL WRITE (lama,val(ib),1,0)
 260 CALL WRITE (lama,0,0,1)
 CALL WRITE (lama,loc(ib),nn,1)
 CALL CLOSE (lama,1)
 msg(2) = lama
 msg(3) = nn
 CALL wrttrl (msg(2))
 
!     IF FREQ. DOES NOT START FROM FIRST FUNDAMENTAL MODE, ADJUST VAL
!     AND LOC TABLES SO THAT WILVEC WILL PICK UP FREQUENCIES CORRECTLY
 
 IF (ib <= 1) GO TO 280
 j = 1
 DO  i = ib,n
   val(j) = val(i)
   loc(j) = loc(i)
   j = j + 1
 END DO
 
 280 IF (nv == 0 .AND. ib > 1 .AND. nn < nfound .AND. val(1) <= zero)  &
     nfound = 0
 IF (nv == 0) nv = nn
 RETURN
END SUBROUTINE qriter
