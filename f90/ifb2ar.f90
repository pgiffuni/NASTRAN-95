SUBROUTINE ifb2ar (TYPE,ifb,ar,l)
     
!     THIS ROUTINE STORES IN ARRAY AR(L+1) THE BCD VALUE OF IFB, AND
!     UPDATE THE L COUNTER
 
!     IF TYPE=1, IFB IS AN INTEGER, AND 8 DIGITS ARE USED IN AR, AND
!                L IS INCREASED BY 2 (INTEGER IS RIGHT ADJUSTED)
!     IF TYPE=2, IFB IS A REAL NUMBER, 12 DIGITS ARE USED IN AR, AND
!                L IS INCREASED BY 3
!     IF TYPE=3, IFB IS A BCD WORD, 4 LETTERS ARE USE IN AR, AND
!                L IS INCREASED BY 1
 
 
 INTEGER, INTENT(IN)                      :: TYPE
 INTEGER, INTENT(IN)                      :: ifb
 INTEGER, INTENT(OUT)                     :: ar(1)
 INTEGER, INTENT(OUT)                     :: l
 INTEGER :: ia, sub(2),zero(2)
 REAL :: ra,x,xl
 CHARACTER (LEN=7) :: fmtx,FMT(10)
 CHARACTER (LEN=8) :: c8
 CHARACTER (LEN=10) :: fmty,fnt(9)
 CHARACTER (LEN=12) :: c12
 EQUIVALENCE  (ia,ra)
 DATA  FMT / '(F12.9)', '(F12.8)', '(F12.7)', '(F12.6)', '(F12.5)',  &
     '(F12.4)', '(F12.3)', '(F12.2)', '(F12.1)', '(F12.0)'/
 DATA  fnt /'(1X,F11.8)', '(1X,F11.7)', '(1X,F11.6)', '(1X,F11.5)',  &
     '(1X,F11.4)', '(1X,F11.3)', '(1X,F11.2)', '(1X,F11.1)', '(1X,F11.0)'/
 DATA  zero/  4H    ,    4H 0.0  /
 DATA  sub /  4HIFB2,    4HAR    /
 
 k = -1
 j = TYPE + 1
 SELECT CASE ( j )
   CASE (    1)
     GO TO 300
   CASE (    2)
     GO TO 200
   CASE (    3)
     GO TO 300
   CASE (    4)
     GO TO 250
 END SELECT
 100  k = k + 1
 IF (k < 0) THEN
   GO TO   150
 ELSE IF (k == 0) THEN
   GO TO   200
 ELSE
   GO TO   250
 END IF
 150  CALL mesage (-37,0,sub)
 
!     INTEGER, RIGHT ADJUSTED
 
 200  WRITE (c8,210,ERR=300) ifb
 READ  (c8,220) ar(l+1),ar(l+2)
 210  FORMAT (i8)
 220  FORMAT (2A4)
 l = l + 2
 RETURN
 
!     BCD WORD
 
 250  ar(l+1) = ifb
 l = l + 1
 RETURN
 
!     REAL NUMBER
 
 300  ia = ifb
 x  = ABS(ra)
 IF (x < 1.0E-36) GO TO 390
 xl = ALOG10(x)
 IF (xl > -4.0 .AND. xl < 10.0) IF (xl-1.0) 350,350,330
 310  WRITE  (c12,320,ERR=100) ra
 320  FORMAT (1P,e12.5)
 GO TO 370
 330  i = xl
 IF (ra < 0.) i = i + 1
 IF (i <= 0   .OR.   i > 9 ) GO TO 310
 IF (ra > 0. .AND. xl > 0.) GO TO 340
 fmtx = FMT(i)
 GO TO 360
 340  fmty = fnt(i)
 WRITE (c12,fmty) ra
 GO TO 370
 350  fmtx = FMT(1)
 360  WRITE  (c12,fmtx) ra
 370  READ   (c12,380) (ar(l+j),j=1,3)
 380  FORMAT (3A4)
 GO TO 400
 390  ar(l+1) = zero(1)
 ar(l+2) = zero(1)
 ar(l+3) = zero(2)
 400  l = l + 3
 RETURN
END SUBROUTINE ifb2ar
