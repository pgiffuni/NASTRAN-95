SUBROUTINE hsbg(n,a,ia,b)
     
!     ..................................................................
 
!        SUBROUTINE HSBG
 
!        PURPOSE
!           TO REDUCE A REAL MATRIX INTO UPPER ALMOST TRIANGULAR FORM
 
!        USAGE
!           CALL HSBG(N,A,IA)
 
!        DESCRIPTION OF THE PARAMETERS
!           N      ORDER OF THE MATRIX
!           A      THE INPUT MATRIX, N BY N
!           IA     SIZE OF THE FIRST DIMENSION ASSIGNED TO THE ARRAY
!                  A IN THE CALLING PROGRAM WHEN THE MATRIX IS IN
!                  DOUBLE SUBSCRIPTED DATA STORAGE MODE.  IA=N WHEN
!                  THE MATRIX IS IN SSP VECTOR STORAGE MODE.
 
!        REMARKS
!           THE HESSENBERG FORM REPLACES THE ORIGINAL MATRIX IN THE
!           ARRAY A.
 
!        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
!           NONE
 
!        METHOD
!           SIMILARITY TRANSFORMATIONS USING ELEMENTARY ELIMINATION
!           MATRICES, WITH PARTIAL PIVOTING.
 
!        REFERENCES
!           J.H. WILKINSON - THE ALGEBRAIC EIGENVALUE PROBLEM -
!           CLARENDON PRESS, OXFORD, 1965.
 
!     ..................................................................
 
 
 INTEGER, INTENT(IN)                      :: n
 DOUBLE PRECISION, INTENT(OUT)            :: a(1)
 INTEGER, INTENT(IN)                      :: ia
 REAL, INTENT(IN)                         :: b(1)
 
 DOUBLE PRECISION :: piv,s,t
 
!     MAKE THIS ROUTINE DOUBLE A AND B ARE SAME SPACE
 
 n2=n*n
 k=n2
 DO  i=1,n2
   a(k)=b(k)
   k=k-1
 END DO
 l=n
 nia=l*ia
 lia=nia-ia
 
!        L IS THE ROW INDEX OF THE ELIMINATION
 
 20 IF(l-3 < 0) THEN
   GO TO   360
 END IF
 40 lia=lia-ia
 l1=l-1
 l2=l1-1
 
!        SEARCH FOR THE PIVOTAL ELEMENT IN THE LTH ROW
 
 isub=lia+l
 ipiv=isub-ia
 piv=DABS(a(ipiv))
 IF(l-3 > 0) THEN
   GO TO    50
 ELSE
   GO TO    90
 END IF
 50 m=ipiv-ia
 DO  i=l,m,ia
   t=DABS(a(i))
   IF(t-piv > 0.0) THEN
     GO TO    60
   ELSE
     GO TO    80
   END IF
   60 ipiv=i
   piv=t
 END DO
 90 IF(piv == 0.0) THEN
   GO TO   320
 END IF
 100   IF(piv-DABS(a(isub)) > 0.0) THEN
   GO TO   120
 ELSE
   GO TO   180
 END IF
 
!        INTERCHANGE THE COLUMNS
 
 120 m=ipiv-l
 DO  i=1,l
   j=m+i
   t=a(j)
   k=lia+i
   a(j)=a(k)
   a(k)=t
 END DO
 
!        INTERCHANGE THE ROWS
 
 m=l2-m/ia
 DO  i=l1,nia,ia
   t=a(i)
   j=i-m
   a(i)=a(j)
   a(j)=t
 END DO
 
!        TERMS OF THE ELEMENTARY TRANSFORMATION
 
 180 DO  i=l,lia,ia
   a(i)=a(i)/a(isub)
 END DO
 
!        RIGHT TRANSFORMATION
 
 j=-ia
 DO  i=1,l2
   j=j+ia
   lj=l+j
   DO  k=1,l1
     kj=k+j
     kl=k+lia
     a(kj)=a(kj)-a(lj)*a(kl)
   END DO
 END DO
 
!        LEFT TRANSFORMATION
 
 k=-ia
 DO  i=1,n
   k=k+ia
   lk=k+l1
   s=a(lk)
   lj=l-ia
   DO  j=1,l2
     jk=k+j
     lj=lj+ia
     s=s+a(lj)*a(jk)
   END DO
   a(lk)=s
 END DO
 
!        SET THE LOWER PART OF THE MATRIX TO ZERO
 
 DO  i=l,lia,ia
   a(i)=0.0
 END DO
 320 l=l1
 GO TO 20
 360 RETURN
END SUBROUTINE hsbg
