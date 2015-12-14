SUBROUTINE fbsi3 (BLOCK, y, mem, dmem, ibuff)
     
!     FBSI3 EXECUTES THE FORWARD/BACKWARD PASS FOR FBSI IN CSP
 
 
 INTEGER, INTENT(OUT)                     :: BLOCK(8)
 COMPLEX, INTENT(IN OUT)                  :: y(1)
 INTEGER, INTENT(IN)                      :: mem(2)
 REAL, INTENT(IN)                         :: dmem(2)
 INTEGER, INTENT(IN OUT)                  :: ibuff(2)
INTEGER :: dbl, buf(2), subnam, begn, END
INTEGER :: dbu, dbb, dbc
INTEGER :: rd, rdrew, wrt, wrtrew, rew
REAL :: l
COMPLEX :: yjk, sum, zero, ljj
CHARACTER (LEN=25) :: sfm
CHARACTER (LEN=29) :: uim
CHARACTER (LEN=25) :: uwm
CHARACTER (LEN=23) :: ufm
COMMON /names /  rd, rdrew, wrt, wrtrew, rew
COMMON /xmssg /  ufm, uwm, uim, sfm
COMMON /zzzzzz/  l(2)
COMMON /system/  sysbuf, nout
COMMON /fbsx  /  dbl(7), dbu(7), dbb(7), dbc(7)
COMMON /fbsm  /  nvec  , nvecsz, nwds  , lasind, ipos(7)
DATA             zero / (0.0, 0.0 )    /
DATA    subnam, begn, END / 4HFBS4, 4HBEGN, 4HEND /

ncol   = dbl(2)
buf(1) = subnam
buf(2) = begn
iopen  = 0
CALL conmsg (buf,2,0)
last   = nvec * nvecsz
nidlt  = 1
lcol   = ipos( 1 )
DO  j = 1,lcol
!      PRINT *,' FORWARD, PROCESSING COLUMN J=',J
  j1 = j - 1
  
! CHECK IF THIS ROW VALUE IS ZERO FOR ALL RIGHT HAND VECTORS
  
  DO  k = j,last,nvecsz
    IF (y(k) /= zero) GO TO 100
  END DO
  
! ALL VALUES FOR THIS ROW ARE ZERO, SKIP TO NEXT ROW OF RIGHT HAND VECTORS
  
  IF ( nidlt >= lasind ) EXIT
  kcol   = mem( nidlt )
  IF ( kcol /= j ) GO TO 7001
  40    nrows  = mem( nidlt+1 )
  nidlt  = nidlt + nrows*nwds + 4
  IF ( nidlt >= lasind ) EXIT
  kcol = mem( nidlt )
  IF ( kcol /= j ) CYCLE
  GO TO 40
  
!     GET 1ST STRING FOR COLUMN AND SAVE DIAGONAL ELEMENT
  
  100   CONTINUE
  kcol   = mem( nidlt )
  IF ( kcol /= j ) GO TO 7001
  nrows  = mem( nidlt + 1 )
  irow   = mem( nidlt + nrows*nwds + 2 )
  indxi  = nidlt + 2
  indxl  = indxi + nrows*2 - 1
  ljj   = 1.0 / CMPLX( dmem( indxi ), dmem( indxi+1 ) )
  IF (nrows == 1) GO TO 600
  indxi = indxi + 2
  irow  = irow + 1
  
!     PROCESS CURRENT STRING IN TRIANGULAR FACTOR AGAINST EACH
!     LOAD VECTOR IN CORE -- Y(I,K) = Y(I,K) + L(I,J)*Y(J,K)
  
  300   DO  k   = 1, last, nvecsz
    yjk        = y( j1+k )
    IF ( yjk == zero ) CYCLE
    iyrow      = irow + k - 1
    DO  ij  = indxi, indxl, 2
      y( iyrow ) = y( iyrow ) + CMPLX( dmem(ij), dmem(ij+1) ) * yjk
      iyrow      = iyrow + 1
    END DO
  END DO
  
!     GET NEXT STRING IN TRIANGULAR FACTOR
  
  600   CONTINUE
  nidlt  = nidlt + 4 + nrows*nwds
  IF ( nidlt >= lasind ) GO TO 800
  kcol   = mem( nidlt )
  IF ( kcol /= j ) GO TO 800
  nrows  = mem( nidlt + 1 )
  irow   = mem( nidlt + nrows*nwds + 2 )
  indxi  = nidlt + 2
  indxl  = indxi + nrows*2 - 1
  GO TO 300
  
!     END-OF-COLUMN ON TRIANGULAR FACTOR -- DIVIDE BY DIAGONAL
  
  800   DO  k = j,last,nvecsz
    y(k) = y(k)*ljj
  END DO
  
END DO
1005  CONTINUE
IF ( lcol == ncol ) GO TO 2005
ifcol = lcol + 1
CALL gopen  ( dbl, ibuff, rdrew )

! POSITION FILE TO APPROPRIATE COLUMN TO BE READ

CALL dsspos ( dbl, ipos(2), ipos(3), ipos(4) )
DO  j = ifcol, ncol
  j1 = j - 1
  
! CHECK IF THIS ROW VALUE IS ZERO FOR ALL RIGHT HAND VECTORS
  
  DO  k = j,last,nvecsz
    IF (y(k) /= zero) GO TO 1100
  END DO
  
! ALL VALUES FOR THIS ROW ARE ZERO, SKIP TO NEXT ROW OF RIGHT HAND VECTORS
  
  CALL skprec ( dbl, 1 )
  CYCLE
  
!     GET 1ST STRING FOR COLUMN AND SAVE DIAGONAL ELEMENT
  
  1100  CONTINUE
  BLOCK(8) = -1
  CALL getstr ( *7002, BLOCK )
  IF (BLOCK(4) /= j) GO TO 7002
  irow  = BLOCK(4)
  indxi = BLOCK(5)
  nrows = BLOCK(6)
  indxl = indxi + nrows*2 - 1
  1200  CONTINUE
  ljj   = 1.0 / CMPLX( l( indxi ), l( indxi+1 ) )
  IF (nrows == 1) GO TO 1600
  indxi = indxi + 2
  irow  = irow + 1
  
!     PROCESS CURRENT STRING IN TRIANGULAR FACTOR AGAINST EACH
!     LOAD VECTOR IN CORE -- Y(I,K) = Y(I,K) + L(I,J)*Y(J,K)
  
  1300  DO  k   = 1, last, nvecsz
    yjk        = y( j1+k )
    IF ( yjk == zero ) CYCLE
    iyrow      = irow + k - 1
    DO  ij  = indxi, indxl, 2
      y( iyrow ) = y( iyrow ) + CMPLX( l(ij), l(ij+1) ) * yjk
      iyrow      = iyrow + 1
    END DO
  END DO
  
!     GET NEXT STRING IN TRIANGULAR FACTOR
  
  1600  CONTINUE
  CALL endget ( BLOCK )
  CALL getstr ( *1800, BLOCK )
  irow  = BLOCK(4)
  indxi = BLOCK(5)
  nrows = BLOCK(6)
  indxl = indxi + nrows*2 - 1
  GO TO 1300
  
!     END-OF-COLUMN ON TRIANGULAR FACTOR -- DIVIDE BY DIAGONAL
  
  1800  DO  k = j,last,nvecsz
    y(k) = y(k)*ljj
  END DO
END DO
2005  CONTINUE
IF ( ncol  == 1 ) GO TO 7000
j = ncol - 1
IF ( lcol == ncol ) GO TO 3000

!     INITIALIZE FOR BACKWARD PASS BY SKIPPING THE NTH COLUMN

CALL bckrec (BLOCK)

!     GET A STRING IN CURRENT COLUMN. IF THIS STRING INCLUDES DIAGONAL,
!     ADJUST STRING TO SKIP IT.

2200  j1 = j - 1
BLOCK(8) = -1
2300  CALL getstb (*2900,BLOCK)
irow  = BLOCK( 4 )
nrows = BLOCK( 6 )
IF (irow-nrows == j1) nrows = nrows - 1
IF (nrows == 0) GO TO 2800
indxi = BLOCK( 5 )

!     PROCESS CURRENT STRING IN TRIANGULAR FACTOR AGAINST EACH
!     LOAD VECTOR IN CORE -- Y(J,K) = Y(J,K) + L(J,I)*Y(I,K)

DO  k = 1,last,nvecsz
  ji  = indxi + 2
  ik  = irow  + k
  sum = (0.0, 0.0)
  DO  ii = 1,nrows
    ji = ji - 2
    ik = ik - 1
    sum = sum + CMPLX( l(ji),l(ji+1) ) * y(ik)
  END DO
  y(j1+k) = y(j1+k) + sum
END DO

!     TERMINATE CURRENT STRING AND GET NEXT STRING

2800  CONTINUE
CALL endgtb (BLOCK)
GO TO 2300

!     END-OF-COLUMN -- TEST FOR COMPLETION

2900  IF (j == 1) GO TO 7000
j = j - 1
IF ( j == lcol ) GO TO 3010
GO TO 2200

3000  CONTINUE

!     INITIALIZE FOR BACKWARD PASS BY SKIPPING THE NTH COLUMN

3005  CONTINUE
nidlt = nidlt - 1
nrows = mem( nidlt )
nidlt = nidlt - nrows*nwds - 3
kcol  = mem( nidlt )
IF ( kcol == ncol ) GO TO 3005
nidlt = nidlt + nrows*nwds + 4
3010  CONTINUE

!     GET A STRING IN CURRENT COLUMN. IF THIS STRING INCLUDES DIAGONAL,
!     ADJUST STRING TO SKIP IT.

3200  j1 = j - 1
!      print *,' processing column in backward step, j=',j
3250  nidlt = nidlt - 1
IF ( nidlt <= 1 ) GO TO 3900
nrows = mem( nidlt )
irow  = mem( nidlt-1 )
nidlt = nidlt - nrows*nwds - 3
kcol  = mem( nidlt )
3260  CONTINUE
IF ( kcol /= j ) GO TO 3900
indxi = nidlt + nrows*2
irow  = irow + nrows - 1
IF ( (irow-nrows) == j1 ) nrows = nrows - 1

!     PROCESS CURRENT STRING IN TRIANGULAR FACTOR AGAINST EACH
!     LOAD VECTOR IN CORE -- Y(J,K) = Y(J,K) + L(J,I)*Y(I,K)

DO  k = 1,last,nvecsz
  ji  = indxi + 2
  ik  = irow  + k
  sum = 0.0
  DO  ii = 1,nrows
    ji = ji - 2
    ik = ik - 1
    sum = sum + CMPLX( dmem(ji),dmem(ji+1) ) * y(ik)
  END DO
  y(j1+k) = y(j1+k) + sum
END DO

!     TERMINATE CURRENT STRING AND GET NEXT STRING

GO TO 3250

!     END-OF-COLUMN -- TEST FOR COMPLETION

3900  IF (j == 1) GO TO 7000
j  = j - 1
j1 = j - 1
GO TO 3260

7000  buf(2) = END
CALL conmsg (buf,2,0)
CALL CLOSE ( dbl, rew )
RETURN

!     FATAL ERROR MESSAGE

7001  CONTINUE
7002  CONTINUE
WRITE  (nout,9001) sfm,subnam
9001  FORMAT (a25,' 2149, SUBROUTINE ',a4,/5X,'FIRST ELEMENT OF A COLU',  &
    'MN OF LOWER TRIANGULAR MATRIX IS NOT THE DIAGONAL ELEMENT')
CALL mesage (-61,0,0)
RETURN
END SUBROUTINE fbsi3
