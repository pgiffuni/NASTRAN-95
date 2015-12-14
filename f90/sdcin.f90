SUBROUTINE sdcin (BLOCK,ac,n,vecs,vecd)
     
!     SDCIN USES GETSTR/ENDGET TO READ A ROW OF A MATRIX AND ADD THE
!     TERMS OF THE ROW INTO A VECTOR
 
!     BLOCK = A 15-WORD ARRAY IN WHICH BLOCK (1) = GINO NAME
!     AC    = A VECTOR OF N COLUMN POSITIONS (COL NBRS MAY BE .LT. 0)
!     N     = NUMBER OF WORDS IN AC AND NUMBER OF TERMS IN VECS
!     VECS  = A VECTOR OF N TERMS. THE POS OF EACH TERM IS DEFINED BY
!     THE NUMBER STORED IN THE CORRESPONDING POSITION IN AC
!     VECD  = SAME VECTOR AS VECS
 
 
 INTEGER, INTENT(IN)                      :: BLOCK(15)
 INTEGER, INTENT(IN OUT)                  :: ac(1)
 INTEGER, INTENT(IN)                      :: n
 REAL, INTENT(OUT)                        :: vecs(1)
 DOUBLE PRECISION, INTENT(OUT)            :: vecd(1)
 INTEGER :: prc     ,words    ,rlcmpx   ,TYPE   , rc       ,prec
 REAL :: xns(1)
 DOUBLE PRECISION :: xnd
 COMMON /TYPE  /  prc(2)   ,words(4) ,rlcmpx(4)
 COMMON /system/  sysbuf   ,nout
 COMMON /zzzzzz/  xnd(1)
 EQUIVALENCE      (xnd(1),xns(1))
 
!     PERFORM GENERAL INITIALIZATION
 
 TYPE = BLOCK(2)
 prec = prc(TYPE)
 rc   = rlcmpx(TYPE)
 i    = 1
 
!     LOCATE POSITION IN VECTOR CORRESPONDING TO STRING
 
 10 IF (i > n) GO TO 92
 DO  j = i,n
   IF (IABS(ac(j)) == BLOCK(4)) GO TO 12
 END DO
 GO TO 90
 12 i = j + BLOCK(6)
 nn = BLOCK(4) + BLOCK(6) - 1
 IF (IABS(ac(i-1)) /= nn) GO TO 91
 
!     ADD TERMS FROM STRING INTO VECTOR
 
 ii = rc*(j-1)
 jstr = BLOCK(5)
 nstr = jstr + rc*BLOCK(6) - 1
 IF (prec == 2) GO TO 24
 
 DO  jj = jstr,nstr
   ii = ii + 1
   vecs(ii) = vecs(ii) + xns(jj)
 END DO
 GO TO 30
 
 24 DO  jj = jstr,nstr
   ii = ii + 1
   vecd(ii) = vecd(ii) + xnd(jj)
 END DO
 
!     CLOSE CURRENT STRING AND GET NEXT STRING
 
 30 CALL endget (BLOCK)
 CALL getstr (*99,BLOCK)
 GO TO 10
 
!     LOGIC ERRORS
 
 90 kerr = 1
 GO TO 97
 91 kerr = 2
 GO TO 97
 92 kerr = 3
 GO TO 97
 97 WRITE  (nout,98) kerr
 98 FORMAT (22H0*** sdcin fatal error ,i2)
 CALL mesage (-61,0,0)
 99 RETURN
END SUBROUTINE sdcin
