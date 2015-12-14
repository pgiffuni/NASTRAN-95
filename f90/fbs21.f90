SUBROUTINE fbs21 (BLOCK,y,yn,nwds)
     
!     FBS2 EXECUTES THE FORWARD/BACKWARD PASS FOR FBS IN RSP
!                                                        ===
 
 
 INTEGER, INTENT(OUT)                     :: BLOCK(20)
 REAL, INTENT(IN OUT)                     :: y(1)
 REAL, INTENT(IN OUT)                     :: yn(1)
 INTEGER, INTENT(IN)                      :: nwds
INTEGER :: dbl,buf(3),subnam(2),begn,END

DOUBLE PRECISION :: ljj,l
CHARACTER (LEN=25) :: sfm
CHARACTER (LEN=29) :: uim
CHARACTER (LEN=25) :: uwm
CHARACTER (LEN=23) :: ufm
COMMON /xmssg /  ufm,uwm,uim,sfm
COMMON /system/  sysbuf,nout
COMMON /zzzzzz/  l(1)
COMMON /fbsx  /  dbl, n
DATA    subnam,  begn, END /4HFBS2, 4H1   , 4HBEGN, 4HEND /

buf(1) = subnam(1)
buf(2) = subnam(2)
buf(3) = begn
CALL conmsg (buf,3,0)
nbritm = nwds/2
nbrvec = (locfx(yn) - locfx(y))/nwds + 1
last = 1 + (nbrvec-1)*nbritm
DO  j=1,n
  
!     MAKE 1ST STRING CALL FOR COLUMN AND SAVE DIAGONAL ELEMENT
  
  BLOCK(8) = -1
  CALL getstr (*81,BLOCK)
  IF (BLOCK(4) /= j) GO TO 81
  jstr = BLOCK(5)
  ljj  = l(jstr)
!WKBI
  xljj = ljj
  IF (BLOCK(6) == 1) GO TO 20
  nstr = jstr + BLOCK(6) - 1
  jstr = jstr + 1
  BLOCK(4) = BLOCK(4) + 1
  
!     PROCESS CURRENT STRING IN TRIANGULAR FACTOR AGAINST EACH
!     LOAD VECTOR IN CORE -- Y(I,K) = Y(I,K) + L(I,J)*Y(J,K)
  
  10 DO  k = 1,last,nbritm
    yjk = y(j+k-1)
    ik  = BLOCK(4) + k - 1
    DO  ij = jstr,nstr
!WKBI
      xlij = l(ij)
!WKBR Y(IK) = Y(IK) + L(IJ)*YJK
      y(ik) = y(ik) + xlij*yjk
      ik = ik + 1
    END DO
  END DO
  
!     GET NEXT STRING IN TRIANGULAR FACTOR
  
  20 CALL endget (BLOCK)
  CALL getstr (*30,BLOCK)
  jstr = BLOCK(5)
  nstr = jstr + BLOCK(6) - 1
  GO TO 10
  
!     END-OF-COLUMN ON TRIANGULAR FACTOR -- DIVIDE BY DIAGONAL
  
  30 DO  k = 1,last,nbritm
!WKBR Y(J+K-1) = Y(J+K-1)/LJJ
    y(j+k-1) = y(j+k-1)/xljj
  END DO
END DO

!     INITIALIZE FOR BACKWARD PASS BY SKIPPING THE NTH COLUMN

IF (n == 1) GO TO 65
CALL bckrec (BLOCK)
j = n - 1

!     GET A STRING IN CURRENT COLUMN. IF STRING INCLUDES DIAGONAL,
!     ADJUST STRING TO SKIP IT.

40 BLOCK(8) = -1
42 CALL getstb (*60,BLOCK)
IF (BLOCK(4)-BLOCK(6)+1 == j) BLOCK(6) = BLOCK(6) - 1
IF (BLOCK(6) == 0) GO TO 59
nterms = BLOCK(6)

!     PROCESS CURRENT STRING IN TRIANGULAR FACTOR AGAINST EACH
!     LOAD VECTOR IN CORE -- Y(J,K) = Y(J,K) + L(J,I)*Y(I,K)

DO  k = 1,last,nbritm
  ji = BLOCK(5)
  ik = BLOCK(4) + k - 1
  jk = j + k - 1
  DO  ii = 1,nterms
    y(jk) = y(jk) + l(ji)*y(ik)
    ji = ji - 1
    ik = ik - 1
  END DO
END DO

!     TERMINATE CURRENT STRING AND GET NEXT STRING

59 CALL endgtb (BLOCK)
GO TO 42

!     END-OF-COLUMN -- TEST FOR COMPLETION

60 IF (j /= 1) GO TO 70
65 buf(3) = END
CALL conmsg (buf,3,0)
RETURN

70 j = j - 1
GO TO 40

!     FATAL ERROR MESSAGE

81 WRITE  (nout,82) sfm,subnam
82 FORMAT (a25,' 2149, SUBROUTINE ',2A4,/5X,'FIRST ELEMENT OF A COL',  &
    'UMN OF LOWER TRIANGULAR MATRIX IS NOT THE DIAGONAL ELEMENT')
CALL mesage (-61,0,0)
RETURN
END SUBROUTINE fbs21
