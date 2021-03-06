SUBROUTINE fbs3 (BLOCK,y,yn,nwds)
     
!     FBS3 EXECUTES THE FORWARD/BACKWARD PASS FOR FBSF IN CSP
 
 
 INTEGER, INTENT(OUT)                     :: BLOCK(8)
 REAL, INTENT(IN OUT)                     :: y(1)
 REAL, INTENT(IN OUT)                     :: yn(1)
 INTEGER, INTENT(IN)                      :: nwds
INTEGER :: dbl, buf(2), subnam, begn, END
REAL :: ljjr, ljji, l
CHARACTER (LEN=25) :: sfm
CHARACTER (LEN=29) :: uim
CHARACTER (LEN=25) :: uwm
CHARACTER (LEN=23) :: ufm
COMMON /xmssg / ufm, uwm, uim, sfm
COMMON /system/ sysbuf, nout
COMMON /zzzzzz/ l(1)
COMMON /fbsx  / dbl   , n
EQUIVALENCE     (sumr,yjkr),  (sumi,yjki)
DATA    subnam, begn  , END / 4HFBS3, 4HBEGN, 4HEND /

buf(1) = subnam
buf(2) = begn
CALL conmsg (buf,2,0)
nbritm = nwds
j    = (locfx(yn)-locfx(y)+1)/nwds
last = MAX0(j,1)*nbritm
DO  j = 1,n
  j1 = j - 1
  DO  k = 1,last,nbritm
    yjkr = y(2*j+k-2)
    yjki = y(2*j+k-1)
    IF (yjkr /= 0.0 .OR. yjki /= 0.0) GO TO 7
  END DO
  CALL skprec (BLOCK(1),1)
  CYCLE
  
!     MAKE 1ST CALL FOR COLUMN AND SAVE DIAGONAL ELEMENT
  
  7 BLOCK(8) = -1
  CALL getstr (*80,BLOCK)
  IF (BLOCK(4) /= j) GO TO 80
  jstr = BLOCK(5)
  ljjr = l(jstr  )
  ljji = l(jstr+1)
  IF (BLOCK(6) == 1) GO TO 20
  nstr = jstr + 2*BLOCK(6) - 2
  jstr = jstr + 2
  BLOCK(4) = BLOCK(4) + 1
  
!     PROCESS CURRENT STRING IN TRIANGULAR FACTOR AGAINST EACH
!     LOAD VECTOR IN CORE -- Y(I,K) = Y(I,K) + L(I,J)*Y(J,K)
  
  10 DO  k = 1,last,nbritm
    yjkr = y(2*j+k-2)
    yjki = y(2*j+k-1)
    IF (yjkr == 0.0 .AND. yjki == 0.0) CYCLE
    ik   = 2*BLOCK(4) + k - 2
    DO  ij = jstr,nstr,2
      y(ik  ) =  y(ik  ) + l(ij)*yjkr - l(ij+1)*yjki
      y(ik+1) =  y(ik+1) + l(ij)*yjki + l(ij+1)*yjkr
      ik = ik + 2
    END DO
  END DO
  
!     GET NEXT STRING IN TRIANGULAR FACTOR
  
  20 CALL endget (BLOCK)
  CALL getstr (*30,BLOCK)
  jstr = BLOCK(5)
  nstr = jstr + 2*BLOCK(6) - 2
  GO TO 10
  
!     END-OF-COLUMN ON TRIANGULAR FACTOR -- DIVIDE BY DIAGONAL
  
  30 ssqr = 1.0/(ljjr**2 + ljji**2)
  DO  k = 1,last,nbritm
    yjkr       = (y(2*j+k-2)*ljjr + y(2*j+k-1)*ljji)*ssqr
    y(2*j+k-1) =-(y(2*j+k-2)*ljji - y(2*j+k-1)*ljjr)*ssqr
    y(2*j+k-2) = yjkr
  END DO
  
END DO

!     INITIALIZE FOR BACKWARD PASS BY SKIPPING THE NTH COLUMN

IF (n == 1) GO TO 65
CALL bckrec (BLOCK)
j = n - 1

!     GET A STRING IN CURRENT COLUMN.  IF STRING INCLUDES DIAGONAL,
!     ADJUST STRING TO SKIP IT.

40 j1 = j - 1
BLOCK(8) = -1
42 CALL getstb (*60,BLOCK)
IF (BLOCK(4)-BLOCK(6) == j1) BLOCK(6) = BLOCK(6) - 1
IF (BLOCK(6) == 0) GO TO 58
nterms = BLOCK(6)

!     PROCESS CURRENT STRING IN TRIANGULAR FACTOR AGAINST EACH
!     LOAD VECTOR IN CORE -- Y(J,K) = Y(J,K) + L(J,I)*Y(I,K)

DO  k = 1,last,nbritm
  ji = BLOCK(5) + 2
  ik = BLOCK(4)*2 + k
  jk = j1*2 + k
  sumr = 0.0
  sumi = 0.0
  DO  ii = 1,nterms
    ji = ji - 2
    ik = ik - 2
    sumr = sumr + l(ji)*y(ik  ) - l(ji+1)*y(ik+1)
    sumi = sumi + l(ji)*y(ik+1) + l(ji+1)*y(ik  )
  END DO
  y(jk  ) = y(jk  ) + sumr
  y(jk+1) = y(jk+1) + sumi
END DO

!     TERMINATE CURRENT STRING AND GET NEXT STRING

58 CALL endgtb (BLOCK)
GO TO 42

!     END-OF-COLUMN -- TEST FOR COMPLETION

60 IF (j /= 1) GO TO 70
65 buf(2) = END
CALL conmsg (buf,2,0)
RETURN

70 j = j - 1
GO TO 40

!     FATAL ERROR MESSAGE

80 WRITE  (nout,82) sfm,subnam
82 FORMAT (a25,' 2149, SUBROUTINE ',a4,/5X,'FIRST ELEMENT OF A COLU',  &
    'MN OF LOWER TRIANGULAR MATRIX IS NOT THE DIAGONAL ELEMENT')
CALL mesage (-61,0,0)
RETURN
END SUBROUTINE fbs3
