SUBROUTINE ifppvc (*,ipvs,jr)
     
!     IFPPVC TAKES 1PARM AND 1VARY CARDS AND MAKES A SCRATCH FILE
!     TO USE IN MODIFYING OTHER BULK DATA CARDS
 
 
 INTEGER, INTENT(OUT)                     :: ipvs
 INTEGER, INTENT(IN OUT)                  :: jr(1)
 LOGICAL :: abort
 INTEGER :: dum,t1,BLANK, NAME(2)
 DIMENSION       z(1)
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /system/ ibuf,nout,abort,dum(79),jrun
 COMMON /ifpdta/ id(2),kn,d1(52),m(50),mf(50),m1(35),m1f(35),  &
     d2(3),nopen,d3(6),knt,d4(18)
 COMMON /zzzzzz/ kor(1)
 COMMON /ifpx1 / ncds,t1(2,1)
 EQUIVALENCE     (kor(1), z(1))
 DATA    ncdsmx, ifil,  iplus, istar, nptp,  ithr,  BLANK   /  &
     343,    213,   1H+,   1H*,   4HNPTP,4HTHRU,4H      /
 DATA    ivar,   ipar,  ivar1, ipar1, NAME                  /  &
     4HAVAR, 4HAPAR,4H1VAR,4H1PAR,4HIFPP,4HVC           /
 
 istop = 0
 ics   = 0
 ltj   = 0
 isort = 0
 ipp   = 2*ibuf + 2
 ln    = ipp
 ii    = ipp - 1
 lst   = 0
 nv    = 0
 idon  = 0
 if0   = 0
 ioldn = 0
 isv   = 0
 iplus = khrfn1(BLANK,4,iplus,1)
 istar = khrfn1(BLANK,4,istar,1)
 CALL sswtch (42,l42)
 GO TO 20
 
!     READ NEW CARD
 
 10 CALL READ (*410,*410,nptp,jr,20,1,kdum)
 knt = knt + 1
 20 it  = khrfn1(BLANK,4,jr(1),1)
 IF (it == iplus .OR. it == istar) GO TO 420
 IF (jr(1) == ivar1) GO TO 90
 IF (jr(1) /= ipar1) GO TO 300
 
!     1PARM CARDS
 
 jr(1) = ipar
 IF (l42 == 0) CALL rcard2 (m1,m1f,nw,jr)
 IF (l42 /= 0) CALL rcard  (m1,m1f,nw,jr)
 IF (nw /= 10) GO TO 430
 
!     CHECK FORMAT
 
 IF (m1f(2) /= 1 .OR. m1(3) < 0) GO TO 430
 IF (m1(3) < jrun) GO TO 10
 IF (m1f(3) /= 0 .AND. m1f(3) /= 1) GO TO 430
 IF (m1f(5) /= 0 .AND. m1f(5) /= 1) GO TO 430
 IF (m1f(7) /= 0 .AND. m1f(7) /= 1) GO TO 430
 IF (m1( 4) < 0  .OR. m1( 6) < 0) GO TO 430
 IF (m1( 8) < 0  .OR. m1f(9) /= 0) GO TO 430
 IF (m1f(4) /= 2 .AND. m1f(4) /= 0) GO TO 430
 IF (m1f(6) /= 2 .AND. m1f(6) /= 0) GO TO 430
 IF (m1f(8) /= 2 .AND. m1f(8) /= 0) GO TO 430
 IF (jrun == 0) GO TO 10
 IF (m1(3) /= jrun .AND. isort == 0) GO TO 80
 IF (m1(3) /= jrun) GO TO 40
 
!     FORM LIST OF K  SK PAIRS FOR THIS J
 
 isort = 1
 25 IF (ipp >= nopen) GO TO 440
 DO  i = 3,7,2
   IF (m1f(i) == 0 .AND. m1f(i+1) /= 0) GO TO 430
   IF (m1f(i) == 0) CYCLE
   kor(ipp  ) = m1(i+1)
   kor(ipp+1) = m1(i+2)
   ipp = ipp + 2
 END DO
 GO TO 10
 
!     SORT LIST ERROR IF DUPLICATE
 
 40 i  = ln
 it = 0
 n  = ipp - i
 IF (n < 3) GO TO 50
 CALL sort (0,0,2,-1,kor(i),n)
 it = kor(i)
 j  = n - 1
 DO  k = 2,j,2
   IF (kor(i+k) /= it) GO TO 42
   abort = .true.
   WRITE (nout,450) ufm,it
   42 it = kor(i+k)
 END DO
 50 isort = 0
 IF (ics /= 0) GO TO 55
 lst= n
 ln = ipp
 55 IF (jr(1) == ivar1) GO TO 100
 IF (idon == 1) GO TO 310
 
!     CHECK FOR DUPLICATE K ON 1PARM ON JRUN = 1
 
 80 IF (jrun /= 1) GO TO 10
 IF (ics  /= 0) GO TO 82
 ics = 1
 ltj = m1(3)
 82 IF (ltj == m1(3)) GO TO 25
 ltj = m1(3)
 GO TO 40
 
!     1VARY CARDS START BUILDING SCRATCH FILE
 
 90 IF (isort == 1 .OR. ltj /= 0) GO TO 40
 
!     IF LST = 0 USE ALL DEFAULT VALUES FOR SK
 
 100 ltj = 0
 nv  = nv + 1
 jr(1) = ivar
 IF (l42 == 0) CALL rcard2 (m1,m1f,nw,jr)
 IF (l42 /= 0) CALL rcard  (m1,m1f,nw,jr)
 IF (nw < 10 .OR. nw > 12) GO TO 430
 
!     CHECK FORMAT
 
 IF (m1f(2) /= 3) GO TO 430
 IF (m1f(3) /= 1  .OR. m1f(4) /= 1) GO TO 430
 IF (m1( 5) <= 0  .OR. m1( 6) <= 0) GO TO 430
 IF (m1f(5) /= 0 .AND. m1f(5) /= 2) GO TO 430
 IF (m1f(6) /= 0 .AND. m1f(6) /= 2) GO TO 430
 IF (m1f(7) /= 0 .AND. m1f(7) /= 1) GO TO 430
 IF (m1f(8) /= 0 .AND. m1f(8) /= 1 .AND. m1f(8) /= 3) GO TO 430
 IF (m1f(9) /= 0 .AND. m1f(9) /= 1) GO TO 430
 IF (m1f(7) == 0 .AND. m1f(8) == 0 .AND. m1f(9) == 0) GO TO 430
 IF (m1f(7) == 1 .AND. m1( 9) == 0) GO TO 430
 IF (m1f(8) == 1 .AND. m1(10) == 0) GO TO 430
 i = 0
 IF (m1f(8) == 3) i = 1
 IF (m1f(9) == 1 .AND. m1(i+11) == 0) GO TO 430
 IF (m1f(8) == 3 .AND. m1(10) /= ithr) GO TO 430
 IF (m1f(8) == 3 .AND. m1(9) > 0 .AND. m1(12) < 0) GO TO 430
 IF (m1f(8) == 3 .AND. m1(9) < 0 .AND. m1(12) > 0) GO TO 430
 IF (jrun == 0) GO TO 10
 DO  kn = 1,ncdsmx
   IF (m1(3) == t1(1,kn) .AND. m1(4) == t1(2,kn)) GO TO 110
 END DO
 GO TO 460
 110 IF (kn /= ioldn .AND. ioldn /= 0) GO TO 140
 112 ioldn = kn
 
!     START A LIST WITH THIS NUMONIC
 
 ifield = m1(5)
 k      = m1(6)
 ia     = m1(7)
 ib     = m1(8)
 IF (m1f(8) == 3) GO TO 120
 IF (lst+isv+18 > nopen) GO TO 440
 DO  i = 7,9
   IF (m1f(i) == 0) CYCLE
   kor(ln+isv  ) = kn
   kor(ln+isv+1) = m1(i+2)
   kor(ln+isv+2) = ifield
   kor(ln+isv+3) = k
   kor(ln+isv+4) = ia
   kor(ln+isv+5) = ib
   isv = isv + 6
 END DO
 GO TO 10
 
!     THRU OPTION
 
 120 n1 = m1(9)
 n2 = m1(12)
 IF (n2 >= n1) GO TO 125
 it = n1
 n1 = n2
 n2 = it
 125 IF (lst+isv+(IABS(n2-n1)*6) > nopen) GO TO 440
 130 kor(ln+isv  ) = kn
 kor(ln+isv+1) = n1
 kor(ln+isv+2) = ifield
 kor(ln+isv+3) = k
 kor(ln+isv+4) = ia
 kor(ln+isv+5) = ib
 isv = isv + 6
 n1  = n1  + 1
 IF (n1 <= n2) GO TO 130
 GO TO 10
 
!     THIS TYPE OF CARD IS DONE SORT LIST AND MAKE FILE
!     SORT ON ID THEN FIELD THEN K
 
 140 IF (isv == 6) GO TO 150
 CALL sort (0,0,6,-2,kor(ln),isv)
 CALL sort (0,0,6,-3,kor(ln),isv)
 CALL sort (0,0,6,-4,kor(ln),isv)
 
!     FIX UP CORE FOR THIS BUFFER AND OPEN FILE
 
 150 IF (if0 /= 0) GO TO 160
 ibuf1 = nopen + 2*ibuf
 nopen = nopen - ibuf
 if0   = 1
 IF (lst+isv > nopen) GO TO 440
 CALL OPEN (*470,ifil,kor(ibuf1+1),1)
 
!     TEST FOR DUPLICATE K FOR SAME FIELD AND ID PLUS SORT AND REG
 
 160 IF (isv == 6) GO TO 220
 it  = kor(ln+1)
 ics = kor(ln+2)
 ik  = kor(ln+3)
 DO  i = 7,isv,6
   IF (it == kor(ln+i).AND.ics == kor(ln+i+1).AND.ik == kor(ln+i+2))  &
       GO TO 170
   GO TO 180
   170 abort = .true.
   WRITE (nout,480) ufm,it,ics,ik
   180 IF (it < 0 .AND. kor(ln+i) > 0) GO TO 220
   IF (it > 0 .AND. kor(ln+i) < 0) GO TO 200
   190 it  = kor(ln+i  )
   ics = kor(ln+i+1)
   ik  = kor(ln+i+2)
   CYCLE
   200 j = kor(ln)
   WRITE (nout,485) ufm,t1(1,j),t1(2,j)
   GO TO 190
 END DO
 
!     PUT OUT CARDS SORT TYPE OF IDS (NEG) DO IN REVERSE
!     FIND VALUES OF SK FOR EACH K
 
 220 n = 6
 i = ln
 IF (kor(ln+1) > 0) GO TO 230
 n = -6
 i = ln + isv - 6
 230 a = 0.0
 IF (kor(i+3) == jrun .AND. lst == 0) a = 1.0
 IF (lst == 0) GO TO 250
 DO  k = 1,lst,2
   IF (kor(i+3) /= kor(ii+k)) CYCLE
   a = z(ii+k+1)
   EXIT
 END DO
 250 z(i+3) = a
 it = kor(i+1)
 ics= kor(i+2)
 IF (it > 0 .AND. ics == 2) GO TO 260
 j = ics/10
 j = j*10
 IF (j /= ics) GO TO 270
 260 abort = .true.
 j = kor(ln)
 WRITE (nout,500) ufm,t1(1,j),t1(2,j),it,ics
 GO TO 280
 270 j = (ics-1)/ 10
 j = j*10
 IF (j == ics-1) GO TO 260
 280 CONTINUE
 IF (abort .OR. a == 0.0) GO TO 290
 CALL WRITE (ifil,kor(i),6,0)
 290 i = i + n
 isv = isv - IABS(n)
 IF (isv > 0) GO TO 230
 isv = 0
 IF (idon == 1) GO TO 310
 GO TO 112
 
!     CARDS ARE DONE
 
 300 idon = 1
 IF (jrun == 0) GO TO 320
 IF (nv == 0 .AND. jrun > 0) GO TO 490
 IF (nv == 0) GO TO 310
 GO TO 140
 310 IF (jrun == 0) GO TO 320
 CALL WRITE (ifil,0,0,1)
 CALL CLOSE (ifil,1)
 ipvs = 1
 320 IF (istop == 0) RETURN
 IF (istop == 1) RETURN 1
 
!     ERROR MESSAGES
 
 400 abort = .true.
 IF (idon-1 == 0) THEN
   GO TO   310
 ELSE
   GO TO    10
 END IF
 410 WRITE  (nout,415) ufm
 415 FORMAT (a23,', NO BULK DATA CARDS TO MODIFY.  ERROR IN IFPPVC')
 istop = 1
 GO TO  310
 420 WRITE  (nout,425) ufm,jr
 425 FORMAT (a23,' 312, NO CONTINUATION CARD ALLOWED ON 1PARM OR ',  &
     '1VARY CARDS', /5X,'CARD- ',20A4)
 GO TO  400
 430 i = knt +1
 WRITE  (nout,435) ufm,m1(1),m1(2),i,jr
 435 FORMAT (a23,' 317, ILLEGAL DATA OR FORMAT ON CARD ',2A4,' SORTED',  &
     i8 ,/5X,'CARD- ' ,20A4)
 GO TO  400
 440 WRITE  (nout,445) ufm
 445 FORMAT (a23,' 3008, NOT ENOUGH CORE FOR 1PARM AND 1VARY CARDS')
 GO TO  400
 450 FORMAT (a23,' 314, DUPLICATE OR NO K  ON 1PARM CARDS FOR SOME J ',  &
     'K =',i9)
 460 WRITE  (nout,465) ufm,m1(3),m1(4)
 465 FORMAT (a23,'316, CARD TYPE ',2A4,' NOT LEGAL ON 1VARY')
 GO TO  400
 470 CALL mesage (-1,ifil,NAME)
 480 FORMAT (a23,'314, DUPLICATE K FOR  ID',i9,' FIELD',i9,' K',i9)
 485 FORMAT (a23,'316, ILLEGAL TO USE SORTED COUNT AND REGULAR ID ON ',  &
     'SAME TYPE OF CARD ',2A4)
 490 WRITE  (nout,495) ufm
 495 FORMAT (a23,', NO 1VARY CARDS TO GO WITH 1PARM CARDS.  ERROR IN ',  &
     'IFPPVC')
 IF (isort == 1 .OR. ltj /= 0) GO TO 40
 GO TO  400
 500 FORMAT (a23,'31, CARD TYPE ',2A4,' ID =',i9,' HAS ILLEGAL FIELD', i9)
END SUBROUTINE ifppvc
