SUBROUTINE sma2b (ke,j,ii,ifile,dumdp)
! ******
! SUBROUTINE SMA2B  ADDS A N X N DOUBLE PRECISION MATRIX, KE, TO THE
! SUBMATRIX OF ORDER NROWSC X JMAX, WHICH IS IN CORE.  N IS 1 IF EITHER
! NPVT, THE PIVOT POINT, IS A  SCALAR POINT, OR J, THE SECOND SUBSCRIPT
! OF KE CORRESPONDS TO A SCALAR POINT, OR J .NE. TO ANY ENTRY IN THE
! GPCT.  OTHERWISE N IS 6.
! ******
 
 DOUBLE PRECISION, INTENT(IN)      :: ke(36)
 INTEGER, INTENT(IN)               :: j
 INTEGER, INTENT(IN)               :: ii
 INTEGER, INTENT(IN OUT)           :: ifile
 DOUBLE PRECISION, INTENT(IN OUT)  :: dumdp
 DOUBLE PRECISION                  :: dz(1)
 INTEGER                           :: iz(1)  ,eor , &
                          clsrw  ,clsnrw , frowic , &
                          tnrows ,outrw  ,  ecpt
 COMMON   /BLANK/ nobgg
 COMMON   /system/ isys(21), linkno
 COMMON   /sem   /  mask(3) , lnknos(15)
 
! SMA2 I/O PARAMETERS
 
 COMMON   /sma2io/ ifcstm             ,ifmpt  &
     ,                  ifdit              ,idum1  &
     ,                  ifecpt             ,igecpt  &
     ,                  ifgpct             ,iggpct  &
     ,                  idum2              ,idum3  &
     ,                  ifmgg              ,igmgg  &
     ,                  ifbgg              ,igbgg  &
     ,                  idum4              ,idum5  &
     ,                  inrw               ,outrw  &
     ,                  clsnrw             ,clsrw  &
     ,                  neor               ,eor  &
     ,                  mcbmgg(7)          ,mcbbgg(7)
 
! SMA2 VARIABLE CORE
 
 COMMON   /zzzzzz /  z(1)
 
! SMA2 VARIABLE CORE BOOKKEEPING PARAMETERS
 
 COMMON   /sma2bk/ icstm              ,ncstm  &
     ,                  igpct              ,ngpct  &
     ,                  ipoint             ,npoint  &
     ,                  i6x6m              ,n6x6m  &
     ,                  i6x6b              ,n6x6b
 
! SMA2 PROGRAM CONTROL PARAMETERS
 
 COMMON   /sma2cl/ ioptb              ,bggind  &
     ,                  npvt               ,left  &
     ,                  frowic             ,lrowic  &
     ,                  nrowsc             ,tnrows  &
     ,                  jmax               ,nlinks  &
     ,                  link(10)           ,nogo
 
! ECPT COMMON BLOCK
 
 COMMON   /sma2et/ ecpt(100)
  
 EQUIVALENCE (z(1),iz(1),dz(1))

!     CALL EMG1B AND THEN RETURN IF THIS IS LINK 8.
!     PROCEED NORMALLY FOR OTHER LINKS.
 
 IF (linkno /= lnknos(8)) GO TO 100
 CALL emg1b (ke, j, ii, ifile, dumdp)
 RETURN
 
! DETERMINE WHICH MATRIX IS BEING COMPUTED.
 
 100 ibase = i6x6m
 IF (ifile == ifmgg) GO TO 5
 IF (ioptb < 0) RETURN
 ibase = i6x6b
 
! SEARCH THE GPCT AND FIND AN INDEX M SUCH THAT
! IABS(GPCT(M)) .LE. J .LT. IABS(GPCT(M+1))
 
 5 low = igpct + 1
 lim = ngpct + low - 2
 IF (low > lim) GO TO 15
 DO  i = low,lim
   isave = i
   IF (j >= IABS(iz(i+1)) ) CYCLE
   IF (j >= IABS(iz(i)) ) GO TO 20
 END DO
 IF ( j >= IABS(iz(isave+1)) )  isave = isave + 1
 GO TO 20
 
! IF II .GT. 0, WE ARE DEALING WITH A SCALAR POINT.
 
 15 isave = low
 20 IF (ii > 0) GO TO 60
 
! AT THIS POINT IT HAS BEEN DETERMINED THAT J IS A SCALAR INDEX NUMBER
! WHICH CORRESPONDS TO A GRID POINT.  HENCE THE DOUBLE PRECISION 6 X 6
! MATRIX, KE, WILL BE ADDED TO THE MATRIX.
 
 l1 = frowic - 1
 jj = ipoint + isave - igpct
 j2 = iz(jj) - 1
 i1 = 0
 lim = nrowsc - 1
 30 IF (i1 > lim) RETURN
 k1 = ibase + i1*jmax + j2
 j1 = 0
 l = 6*l1
 k = k1
 40 j1 = j1 + 1
 IF (j1 > 6) GO TO 50
 l = l + 1
 k = k + 1
 dz(k) = dz(k) + ke(l)
 GO TO 40
 50 i1 = i1 + 1
 l1 = l1 + 1
 GO TO 30
 
! AT THIS POINT WE ARE DEALING WITH A 1 X 1.
! FIRST COMPUTE THE ROW NUMBER, NROW
 
 60 nrow = ii - npvt + 1
 
! THE FOLLOWING 2 FORTRAN STATEMENTS ARE MERELY TO CHECK THE PROGRAM
! LOGIC.  EVENTUALLY THEY CAN BE DELETED.
 
 IF (nrow >= 1  .AND.  nrow <= tnrows)  GO TO 70
 CALL mesage (-30,22,ecpt(1))
 70 lrowic = frowic + nrowsc - 1
 
! IF NROW, THE ROW INTO WHICH THE NUMBER KE(1) IS TO BE ADDED IS NOT
! IN CORE IT CANNOT BE ADDED AT THIS TIME.
 
 IF (nrow < frowic  .OR.  nrow > lrowic) RETURN
 j2 = isave
 j3 = ipoint + isave - igpct
 INDEX = ibase + (nrow-1)*jmax + iz(j3) + j - IABS(iz(j2))
 dz(INDEX) = dz(INDEX) + ke(1)
 RETURN
 
END SUBROUTINE sma2b
