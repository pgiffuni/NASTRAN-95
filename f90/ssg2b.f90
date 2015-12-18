SUBROUTINE ssg2b (kfs,cdt,pabar,sr1,t1,iprec1,ia1,sr2)
     
 
 INTEGER, INTENT(IN)                      :: kfs
 INTEGER, INTENT(IN)                      :: cdt
 INTEGER, INTENT(IN)                      :: pabar
 INTEGER, INTENT(IN)                      :: sr1
 INTEGER, INTENT(IN)                      :: t1
 INTEGER, INTENT(IN OUT)                  :: iprec1
 INTEGER, INTENT(IN)                      :: ia1
 INTEGER, INTENT(IN)                      :: sr2
 CHARACTER (LEN=27) :: swm
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm,swm
 COMMON /zzzzzz/ core(1)
 COMMON /system/ ksystm(55)
 COMMON /mpyadx/ filea(7),fileb(7),filec(7),filed(7),nz,t,i1,i2, prec,scr2
 EQUIVALENCE     (ksystm(55),kprec1), (ksystm(1),sysbuf), (ksystm( 2),ioutpt)
 DATA    square, rect,diag,symm,ident / 1,2,3,6,8 /
 
 prec1 = MIN0(kprec1,iprec1)
 IF (prec1 <= 0) prec1 = kprec1
 nz = korsz(core)
 DO  i = 1,21
   filea(i) = 0
 END DO
 filea(1) = kfs
 scr2 = sr2
 IF (IABS(ia1)-1 < 0) THEN
   GO TO    40
 ELSE IF (IABS(ia1)-1 == 0) THEN
   GO TO    20
 ELSE
   GO TO    30
 END IF
 20 i2 = ia1
 i1 = ia1
 GO TO 50
 30 i2 =-1
 i1 = 1
 GO TO 50
 40 i1 =-1
 i2 = 1
 50 CALL rdtrl (filea)
 fileb(1) = cdt
 CALL rdtrl (fileb)
 IF (fileb(1) <= 0) fileb(4) = symm
 filec(1) = pabar
 CALL rdtrl (filec)
 IF (filec(1) <= 0) GO TO 70
 IF (filec(2) == fileb(2) .OR. fileb(1) <= 0) GO TO 80
 WRITE (ioutpt,60) swm,fileb(1),fileb(3),fileb(2),fileb(3),filec(2)
 60 FORMAT (a27,' 2363, SSG2B FORCED MPYAD COMPATIBILITY OF MATRIX ON'  &
     ,       i5,8H, from (,i5,1H,,i5,7H), TO (,i5,1H,,i5,1H))
 fileb(2) = filec(2)
 GO TO 80
 70 filec(1) = 0
 filec(4) = diag
 80 filed(4) = rect
 filed(1) = sr1
 
!     COMPUTE TYPE OF OUTPUT
 
 irc = 0
 IF (filea(5) > 2 .OR. fileb(5) > 2 .OR. (filec(5) > 2 .AND.  &
     filec(1) /= 0)) irc = 2
 filed(5) = prec1 + irc
 t = t1
 prec = prec1
 filed(3) = filea(3)
 IF (t /= 0) filed(3) = filea(2)
 IF (filea(1) <= 0 .OR. fileb(1) <= 0) filed(3) = filec(3)
 CALL mpyad (core,core,core)
 IF (filed(2) == filed(3) .AND. filed(4) /= symm) filed(4) = square
 IF (filed(4) == symm .OR. filed(4) /= square) GO TO 100
 
!     IF END RESULT IS A SYMMETRIC MATRIX, MAKE SURE THE FORM IS SET TO
!     6 (SYMM). IT COULD SAVE CPU TIME LATER AND WORTH ONE FINAL CHECK.
 
 k = 0
 DO  i = 1,21,7
   IF (filea(i) <= 0) CYCLE
   j = filea(i+3)
   IF (j == diag .AND. i == 15   ) CYCLE
   IF (j /= symm .AND. j /= ident) GO TO 100
   IF (j ==  symm) k = k + 10
   IF (j == ident) k = k + 1
 END DO
 IF (k >  0) filed(4) = ident
 IF (k >= 10) filed(4) = symm
 100 CALL wrttrl (filed)
 RETURN
END SUBROUTINE ssg2b
