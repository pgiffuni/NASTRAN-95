SUBROUTINE pla4b (ke,j)
!*****
! THIS ROUTINE IS THE INSERTION ROUTINE FOR THE PLA4 MODULE.  IT ADDS
! THE 6 X 6 DOUBLE PRECISION MATRIX KE TO THE SUBMATRIX OF ORDER
! 6 X JMAX
!*****
 
 DOUBLE PRECISION, INTENT(IN)             :: ke(36)
 INTEGER, INTENT(IN OUT)                  :: j
 DOUBLE PRECISION :: dz(1)
 
 
 
 INTEGER :: frowic             ,iz(1)
 
! VARIABLE CORE
 
 COMMON   /zzzzzz/ z(1)
 
! PLA42 COMMUNICATIONS BLOCK
 
 COMMON   /pla42c/ idum5(6)  &
     ,                  igpct              ,ngpct  &
     ,                  ipoint             ,npoint  &
     ,                  i6x6k              ,n6x6k ,                  idum11(12)  &
     ,                  jmax               ,frowic  &
     ,                  lrowic             ,nrowsc ,                  idum(121)
 
 
 
 EQUIVALENCE (dz(1),z(1),iz(1))
 
! SEARCH THE GPCT AND FIND AN INDEX M SUCH THAT
! IABS(GPCT(M)) .LE. J .LT. IABS(GPCT(M+1))
 
 low = igpct + 1
 lim = ngpct + low - 2
 IF (low > lim) GO TO 15
 DO  i = low,lim
   isave = i
   IF (j >= IABS(iz(i+1)) ) CYCLE
   IF (j >= IABS(iz(i)) ) GO TO 20
 END DO
 IF (j >= IABS(iz(isave+1)) ) isave = isave + 1
 GO TO 20
 15 isave = low
 
! ADD KE TO THE SUBMATRIX
 
 20 l1 = frowic - 1
 jj = ipoint + isave - igpct
 j2 = iz(jj) - 1
 i1 = 0
 lim = nrowsc - 1
 30 IF (i1 > lim) RETURN
 k1 = i6x6k + i1*jmax + j2
 j1 = 0
 l  = 6*l1
 k  = k1
 40 j1 = j1 + 1
 IF (j1 > 6) GO TO 50
 k  = k + 1
 l  = l + 1
 dz(k) = dz(k) + ke(l)
 GO TO 40
 50 i1 = i1 + 1
 l1 = l1 + 1
 GO TO 30
END SUBROUTINE pla4b
