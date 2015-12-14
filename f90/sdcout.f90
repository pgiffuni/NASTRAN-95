SUBROUTINE sdcout (BLOCK,irw,ac,n,vecs,vecd)
     
!     SDCOUT WRITES A ROW OF A MATRIX IN STRING FORMAT USING
!     PUTSTR/ENDPUT.
 
!     BLOCK = A 15-WORD ARRAY IN WHICH BLOCK(1),(2),(3) HAVE ALREADY
!             BEEN COMPLETED WITH GINO NAME, TYPE AND FORMAT
!     IRW   = ZERO -- ROW NBR OF VECTOR = AC(1)
!           = N.Z. -- ROW NBR OF VECTOR IS IRW
!     AC    = A VECTOR OF N COLUMN POSITIONS (COL NBRS MAY BE .LT. 0)
!     N     = NUMBER OF WORDS IN AC AND NUMBER OF TERMS IN VECS
!     VECS  = A VECTOR OF N TERMS. THE POS OF EACH TERM IS DEFINED
!             BY THE NUMBER STORED IN THE CORRESPONDING POSITION IN AC
!     VECD  = SAME VECTOR AS VECS
 
 
 INTEGER, INTENT(OUT)                     :: BLOCK(15)
 INTEGER, INTENT(IN)                      :: irw
 INTEGER, INTENT(IN)                      :: ac(1)
 INTEGER, INTENT(IN OUT)                  :: n
 REAL, INTENT(IN)                         :: vecs(1)
 DOUBLE PRECISION, INTENT(IN)             :: vecd(1)
 INTEGER :: prc      ,words    ,rlcmpx   ,TYPE   , rc       ,prec
 REAL :: xns(1)
 DOUBLE PRECISION :: xnd
 COMMON /TYPE  /  prc(2)   ,words(4) ,rlcmpx(4)
 COMMON /zzzzzz/  xnd(1)
 EQUIVALENCE      (xnd(1),xns(1))
 
 BLOCK(8)  = -1
 BLOCK(12) = irw
 IF (irw == 0) BLOCK(12) = IABS(ac(1))
 ii   = 0
 TYPE = BLOCK(2)
 rc   = rlcmpx(TYPE)
 prec = prc(TYPE)
 i    = 1
 
!     DETERMINE LENGTH OF A STRING BY SCANNING AC
 
 10 BLOCK(4) = IABS(ac(i))
 j = BLOCK(4) - i
 k = i + 1
 12 IF (k > n) GO TO 14
 IF (IABS(ac(k)) /= j+k) GO TO 14
 k = k + 1
 GO TO 12
 14 nbrstr = k - i
 
!     WRITE STRING WITH PUTSTR/ENDPUT
 
 15 CALL putstr (BLOCK)
 BLOCK(7) = MIN0(BLOCK(6),nbrstr)
 jstr = BLOCK(5)
 nstr = jstr + rc*BLOCK(7) - 1
 IF (prec == 2) GO TO 18
 
 DO  jj = jstr,nstr
   ii = ii + 1
   xns(jj) = vecs(ii)
 END DO
 GO TO 22
 
 18 DO  jj = jstr,nstr
   ii = ii + 1
   xnd(jj) = vecd(ii)
 END DO
 
!     TEST FOR COMPLETION
 
 22 i = i + BLOCK(7)
 IF (i > n) GO TO 30
 CALL endput (BLOCK)
 IF (nbrstr == BLOCK(7)) GO TO 10
 nbrstr   = nbrstr - BLOCK(7)
 BLOCK(4) = IABS( ac(i) )
 GO TO 15
 
!     END LAST STRING
 
 30 BLOCK(8) = 1
 CALL endput (BLOCK)
 RETURN
END SUBROUTINE sdcout
