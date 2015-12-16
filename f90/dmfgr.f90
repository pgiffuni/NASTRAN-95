SUBROUTINE dmfgr (a,m,n,eps,irank,irow,icol)
     
!     DMFGR CALCULATES THE RANK AND LINEARLY INDEPENDENT ROWS AND
!     COLUMNS OF A M BY N MATRIX.  IT EXPRESSES A SUBMATRIX OF
!     MAXIMAL RANK AS A PRODUCT OF TRIANGULAR FACTORS, NONBASIC ROWS
!     IN TERMS OF BASIC ONES AND BASIC VARIABLES IN TERMS OF FREE ONES
 
!     DIMENSIONED DUMMY VARIABLES
 
 
 DOUBLE PRECISION, INTENT(IN OUT)         :: a(1)
 INTEGER, INTENT(IN)                      :: m
 INTEGER, INTENT(IN)                      :: n
 REAL, INTENT(IN OUT)                     :: eps
 INTEGER, INTENT(OUT)                     :: irank
 INTEGER, INTENT(OUT)                     :: irow(1)
 INTEGER, INTENT(OUT)                     :: icol(1)
 
 DOUBLE PRECISION :: piv,hold,SAVE
 
!     TEST OF SPECIFIED DIMENSIONS
 
 IF (m > 0) THEN
   GO TO    10
 ELSE
   GO TO    20
 END IF
 10 IF (n > 0) THEN
   GO TO    40
 END IF
 20 irank = -1
 30 RETURN
 
!     RETURN IN CASE OF FORMAL ERRORS
 
!     INITIALIZE COLUMN INDEX VECTOR
!     SEARCH FIRST PIVOT ELEMENT
 
 40 irank = 0
 piv = 0.d0
 jj  = 0
 DO  j = 1,n
   icol(j) = j
   DO  i = 1,m
     jj  = jj + 1
     hold = a(jj)
     IF (DABS(piv)-DABS(hold) < 0.0) THEN
       GO TO    50
     ELSE
       GO TO    60
     END IF
     50 piv = hold
     ir  = i
     ic  = j
   60 CONTINUE
   END DO
 END DO
 
!     INITIALIZE ROW INDEX VECTOR
 
 DO  i = 1,m
   irow(i) = i
 END DO
 
!     SET UP INTERNAL TOLERANCE
 
 tol = ABS(eps*SNGL(piv))
 
!     INITIALIZE ELIMINATION LOOP
 
 nm = n*m
 DO  ncol = m,nm,m
   
!     TEST FOR FEASIBILITY OF PIVOT ELEMENT
   
   IF (ABS(SNGL(piv))-tol > 0.0) THEN
     GO TO    90
   ELSE
     GO TO   220
   END IF
   
!     UPDATE RANK
   
   90 irank = irank + 1
   
!     INTERCHANGE ROWS IF NECESSARY
   
   jj = ir - irank
   IF (jj > 0) THEN
     GO TO   100
   ELSE
     GO TO   120
   END IF
   100 DO  j = irank,nm,m
     i  = j + jj
     SAVE = a(j)
     a(j) = a(i)
     a(i) = SAVE
   END DO
   
!     UPDATE ROW INDEX VECTOR
   
   jj = irow(ir)
   irow(ir) = irow(irank)
   irow(irank) = jj
   
!     INTERCHANGE COLUMNS IF NECESSARY
   
   120 jj = (ic-irank)*m
   IF (jj > 0) THEN
     GO TO   130
   ELSE
     GO TO   150
   END IF
   130 kk = ncol
   DO  j = 1,m
     i  = kk + jj
     SAVE  = a(kk)
     a(kk) = a(i)
     kk = kk - 1
     a(i) = SAVE
   END DO
   
!     UPDATE COLUMN INDEX VECTOR
   
   jj = icol(ic)
   icol(ic) = icol(irank)
   icol(irank) = jj
   150 kk = irank + 1
   mm = irank - m
   ll = ncol  + mm
   
!     TEST FOR LAST ROW
   
   IF (mm < 0) THEN
     GO TO   160
   ELSE
     GO TO   270
   END IF
   
!     TRANSFORM CURRENT SUBMATRIX AND SEARCH NEXT PIVOT
   
   160 jj   = ll
   SAVE = piv
   piv  = 0.d0
   DO  j = kk,m
     jj   = jj + 1
     hold = a(jj)/SAVE
     a(jj)= hold
     l    = j - irank
     
!     TEST FOR LAST COLUMN
     
     IF (irank-n < 0) THEN
       GO TO   170
     ELSE
       GO TO   200
     END IF
     170 ii = jj
     DO  i = kk,n
       ii = ii + m
       mm = ii - l
       a(ii) = a(ii) - hold*a(mm)
       IF (DABS(a(ii))-DABS(piv) > 0.0) THEN
         GO TO   180
       ELSE
         GO TO   190
       END IF
       180 piv = a(ii)
       ir  = j
       ic  = i
       190 CONTINUE
     END DO
     200 CONTINUE
   END DO
 END DO
 
!     SET UP MATRIX EXPRESSING ROW DEPENDENCIES
 
 220 IF (irank-1 < 0) THEN
   GO TO    30
 ELSE IF (irank-1 == 0) THEN
   GO TO   270
 END IF
 230 ir = ll
 DO  j = 2,irank
   ii = j - 1
   ir = ir - m
   jj = ll
   DO  i = kk,m
     hold = 0.d0
     jj = jj + 1
     mm = jj
     ic = ir
     DO  l = 1,ii
       hold = hold + a(mm)*a(ic)
       ic = ic - 1
       mm = mm - m
     END DO
     a(mm) = a(mm) - hold
   END DO
 END DO
 
!     TEST FOR COLUMN REGULARITY
 
 270 IF (n-irank > 0) THEN
   GO TO   280
 ELSE
   GO TO    30
 END IF
 
!     SET UP MATRIX EXPRESSING BASIC VARIABLES IN TERMS OF FREE
!     PARAMETERS (HOMOGENEOUS SOLUTION).
 
 280 ir = ll
 kk = ll + m
 DO  j = 1,irank
   DO  i = kk,nm,m
     jj = ir
     ll = i
     hold = 0.d0
     ii = j
     290 ii = ii - 1
     IF (ii > 0) THEN
       GO TO   300
     ELSE
       GO TO   310
     END IF
     300 hold = hold - a(jj)*a(ll)
     jj = jj - m
     ll = ll - 1
     GO TO 290
     310 a(ll) = (hold-a(ll))/a(jj)
   END DO
   ir = ir - 1
 END DO

 RETURN
END SUBROUTINE dmfgr
