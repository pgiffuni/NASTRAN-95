SUBROUTINE ds1b (ke,j)
     
!     THIS ROUTINE ADDS THE 6 X 6 DOUBLE PRECISION MATRIX KE TO THE
!     SUBMATRIX OF ORDER NROWSC X JMAX.
 
 
 DOUBLE PRECISION, INTENT(IN)             :: ke(36)
 INTEGER, INTENT(IN OUT)                  :: j
 INTEGER :: cstm  ,mpt   ,dit   ,ecptds,outrw ,eor   , clsrw ,frowic,iz(1)
 DOUBLE PRECISION :: dz(1)
 COMMON /zzzzzz/  z(1)
 COMMON /ds1aaa/  npvt  ,icstm ,ncstm ,igpct ,ngpct ,ipoint,  &
     npoint,i6x6k ,n6x6k ,cstm  ,mpt   ,dit   ,  &
     ecptds,gpct  ,kggd  ,inrw  ,outrw ,eor   ,  &
     neor  ,clsrw ,jmax  ,frowic,lrowic,nrowsc, nlinks,link(10)     ,nogo
 EQUIVALENCE      (dz(1),z(1),iz(1))
 
!     SEARCH THE GPCT AND FIND AN INDEX M SUCH THAT
!     IABS(GPCT(M)) .LE. J .LT. IABS(GPCT(M+1))
 
 low = igpct + 1
 lim = ngpct + low - 2
 IF (low > lim) GO TO 15
 DO  i = low,lim
   isave = i
   IF (j >= IABS(iz(i+1))) CYCLE
   IF (j >= IABS(iz(i  ))) GO TO 20
 END DO
 IF (j >= IABS(iz(isave+1))) isave = isave + 1
 GO TO 20
 15 isave = low
 
!     ADD KE TO THE SUBMATRIX
 
 20 l1  = frowic - 1
 jj  = ipoint + isave - igpct
 j2  = iz(jj) - 1
 i1  = 0
 lim = nrowsc - 1
 30 IF (i1 > lim) RETURN
 k1  = i6x6k + i1*jmax + j2
 j1  = 0
 l   = 6*l1
 k   = k1
 40 j1  = j1 + 1
 IF (j1 > 6) GO TO 50
 k   = k + 1
 l   = l + 1
 dz(k) = dz(k) + ke(l)
 GO TO 40
 50 i1  = i1 + 1
 l1  = l1 + 1
 GO TO 30
END SUBROUTINE ds1b
