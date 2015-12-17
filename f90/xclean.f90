SUBROUTINE xclean
     
 IMPLICIT INTEGER (a-z)
 EXTERNAL        lshift,rshift,andf,orf
 DIMENSION       nclean( 2),ddbn ( 1),dfnu ( 1),fcum ( 1),  &
     fcus  ( 1),fdbn ( 1),fequ ( 1),FILE ( 1),  &
     fknd  ( 1),fmat ( 1),fntu ( 1),fpun ( 1),  &
     fon   ( 1),ford ( 1),minp ( 1),mlsn ( 1),  &
     mout  ( 1),mscr ( 1),sal  ( 1),sdbn ( 1), sntu  ( 1),sord ( 1)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /system/ ibufsz,outtap
 COMMON /xfiat / fiat(7)
 COMMON /xfist / fist
 COMMON /xdpl  / dpd(6)
 COMMON /xsfa1 / md(401),sos(1501),comm(20),xf1at(5)
 EQUIVALENCE             (dpd  (1),dnaf    ),(dpd  (2),dmxlg   ),  &
     (dpd  (3),dculg   ),(dpd  (4),ddbn (1)),(dpd  (6),dfnu(1) ),  &
     (fiat (1),funlg   ),(fiat (2),fmxlg   ),(fiat (3),fculg   ),  &
     (fiat (4),fequ (1)),(fiat (4),FILE (1)),(fiat (4),ford (1)),  &
     (fiat (5),fdbn (1)),(fiat (7),fmat (1)),(md   (1),mlgn    ),  &
     (md   (2),mlsn (1)),(md   (3),minp (1)),(md   (4),mout (1)),  &
     (md   (5),mscr (1)),(sos  (1),slgn    ),(sos  (2),sdbn (1)),  &
     (sos  (4),sal  (1)),(sos  (4),sntu (1)),(sos  (4),sord (1)),  &
     (xf1at(1),fntu (1)),(xf1at(1),fon  (1)),(xf1at(2),fpun (1)),  &
     (xf1at(3),fcum (1)),(xf1at(4),fcus (1)),(xf1at(5),fknd (1))
 EQUIVALENCE             (comm (1),almsk   ),(comm (2),apndmk  ),  &
     (comm (3),cursno  ),(comm (4),entn1   ),(comm (5),entn2   ),  &
     (comm (6),entn3   ),(comm (7),entn4   ),(comm (8),flag    ),  &
     (comm (9),fnx     ),(comm(10),lmsk    ),(comm(11),lxmsk   ),  &
     (comm(12),macsft  ),(comm(13),rmsk    ),(comm(14),rxmsk   ),  &
     (comm(15),s       ),(comm(16),scornt  ),(comm(17),tapmsk  ),  &
     (comm(18),thcrmk  ),(comm(19),zap     )
 DATA  nclean / 4HXCLE,4HAN  /
 
!     ENTRY SIZE NUMBERS,  1=FIAT, 2=SOS, 3=MD
 
 ifail  = 0
 entn1x = entn1 - 1
 entn1y = entn1*2
 lmt1   = funlg*entn1
 lmt2   = lmt1 + 1
 lmt3   = fculg*entn1
 flag   = 0
 icursn = lshift(cursno,16)
 
!     MERGE FIAT BY REPLACING ANY UNIQUE FILES WITH MATCHED CURRENT
!     FILES ONLY CURRENT TAIL  AND  EXCEPT EQU FILES
 
 IF (funlg == fculg) GO TO 170
 ASSIGN 170 TO isw
 100 DO  i = lmt2,lmt3,entn1
   trial = andf(rmsk,FILE(i))
   IF (trial == zap) CYCLE
   
!     ERASE SCRATCH AND LTU EXPIRED FILES FROM CURRENT TAIL
   
   k = andf(lmsk,ford(i))
   IF (k == lmsk .OR. (icursn > k .AND. k /= 0 .AND. fequ(i) > 0)) GO TO 152
   DO  j = 1,lmt1,entn1
     IF (trial /= andf(rmsk,FILE(j))) CYCLE
     IF (fequ(j) < 0.0) THEN
       GO TO   160
     ELSE
       GO TO   140
     END IF
   END DO
   CYCLE
   140 k = andf(lmsk,ford(j))
   IF (k /= lmsk .AND. icursn <= k) CYCLE
   lmt4 = i + entn1x
   DO  k = i,lmt4
     FILE(j) = FILE(k)
     FILE(k) = 0
     fntu(j) = fntu(k)
     j = j + 1
   END DO
   j = j - entn1
   GO TO 156
   152 lmt4 = i + entn1x
   DO  k = i,lmt4
     FILE(k) = 0
   END DO
   j = i
   156 CALL xpolck (fdbn(j),fdbn(j+1),ik,l)
   IF (fequ(j) < 0) flag = -1
   IF (ik == 0) CYCLE
   ddbn(l  ) = 0
   ddbn(l+1) = 0
   160 CONTINUE
 END DO
 GO TO isw, (170,310)
 
!     REGENERATE ALL NTU VALUES (AND LTU IF EMPTY) IN FIAT BY SCANNING
!     SOS DELETE FIAT ENTRY IF NOT FOUND, OR A SCRATCH
 
 170 lmt4 = mlgn*entn3
 lmt2 = lmt1 + 1
 
!     FIAT LOOP
 
 loop250:  DO  i = 1,lmt3,entn1
   IF (andf(lmsk,ford(i)) == lmsk) GO TO 220
   trial = fdbn(i)
   IF (trial == 0) CYCLE loop250
   
!     SOS LOOP - BY MODULE
   
   lmt6 = 0
   DO  j = 1,lmt4,entn3
     lmt5 = lmt6 + 1
     lmti = lmt6 + minp(j)*entn2
     lmt6 = lmt6 +(minp(j) + mout(j) + mscr(j))*entn2
     
!     SOS LOOP - BY FILE WITHIN MODULE
     
     DO  k = lmt5,lmt6,entn2
       IF (trial /= sdbn(k) .OR. fdbn(i+1) /= sdbn(k+1)) CYCLE
       IF (andf(rmsk,FILE(i)) == zap) GO TO 190
       IF (k > lmti .OR. fmat(i) /= 0 .OR. fmat(i+1) /= 0 .OR.  &
           fmat(i+2) /= 0) GO TO 190
       IF (entn1 == 11 .AND. (fmat(i+5) /= 0 .OR. fmat(i+6) /= 0 .OR.  &
           fmat(i+7) /= 0)) GO TO 190
       
!     IF FIAT ENTRY IS INPUT WITH ZERO TRAILERS - PURGE IT
       
       IF (i <= lmt1) GO TO 185
       lmti = 0
       FILE(i) = orf(FILE(i),zap)
       GO TO 210
       
!     PURGE FILE --PUT ENTRY AT END OF FIAT
       
       185 IF (fculg == fmxlg) GO TO 186
       nfculg = fculg*entn1 + 1
       ifail  = 0
       fculg  = fculg + 1
       FILE(nfculg  ) = orf(FILE(i),zap)
       fdbn(nfculg  ) = fdbn(i  )
       fdbn(nfculg+1) = fdbn(i+1)
       GO TO 210
       
!     TRY TO PACK FIAT FOR MORE SPACE
       
       186 IF (ifail == 1) GO TO 900
       ifail = 1
       ASSIGN 170 TO ihop
       GO TO 310
       190 fntu(i) = andf(rmsk,mlsn(j))
       IF (andf(lmsk,ford(i)) == 0) ford(i) = orf(ford(i),andf(lmsk,sord(k)))
       CYCLE loop250
     END DO
   END DO
   
!     DELETE FIAT ENTRY (UNLESS LTU YET TO COME)
   
!     HOLD FILES UNTIL LARGEST LTU OF EQUIVALENCED GROUP EXPIRES
   
   IF (fequ(i) >= 0 .AND. icursn > andf(lmsk,ford(i))) GO TO 210
   fntu(i) = rshift(andf(lmsk,ford(i)),16)
   CYCLE loop250
   210 CALL xpolck (fdbn(i),fdbn(i+1),ik,l)
   IF (ik == 0) GO TO 215
   ddbn(l  ) = 0
   ddbn(l+1) = 0
   215 IF (lmti == 0) CYCLE loop250
   220 hold = andf(rxmsk,FILE(i))
   IF (fequ(i) < 0) flag = -1
   lmt6 = i + entn1x
   DO  k = i,lmt6
     FILE(k) = 0
   END DO
   IF (i > lmt1) CYCLE loop250
   FILE(i) = hold
   flag = -1
 END DO loop250
 lmt3 = fculg*entn1
 
!     CHECK EQU FILES FOR BREAKING OF EQU
 
 IF (funlg == fculg) RETURN
 loop300:  DO  i = 1,lmt3,entn1
   IF (fequ(i) >= 0 .OR. andf(lmsk,ford(i)) >= icursn) CYCLE loop300
   DO  j = 1,lmt3,entn1
     IF (fequ(j) >= 0) CYCLE
     IF (i == j) CYCLE
     IF (andf(rmsk,FILE(i)) == andf(rmsk,FILE(j)) .AND.  &
         icursn <= andf(lmsk,ford(j))) CYCLE loop300
     
   END DO
   fequ(i) = andf(almsk,fequ(i))
   flag = -1
 END DO loop300
 
!     IF BREAK HAS OCCURED, REPEAT FIAT MERGE
 
 ASSIGN 451 TO ihop
 IF (flag /= -1) GO TO 310
 ASSIGN 310 TO isw
 GO TO 100
 
!     CLOSE UP FILES(IF ANY) BELOW UNIQUE LENGTH - RESET FCULG
 
 310 lmt7 = lmt3 - entn1x
 lmt3 = lmt7 - 1
 330 IF (lmt7 < lmt2) GO TO 450
 IF (fdbn(lmt7) /= 0) GO TO 350
 lmt7 = lmt7 - entn1
 GO TO 420
 350 DO  i = lmt2,lmt3,entn1
   IF (fdbn(i) /= 0) CYCLE
   lmt4 = i + entn1x
   DO  k = i,lmt4
     FILE(k) = FILE(lmt7)
     FILE(lmt7) = 0
     fntu(k) = fntu(lmt7)
     lmt7 = lmt7 + 1
   END DO
   GO TO 410
 END DO
 GO TO 450
 410 lmt7 = lmt7  - entn1y
 lmt2 = i + entn1
 420 lmt3 = lmt3  - entn1
 fculg= fculg - 1
 GO TO 330
 
!     RESET ANY NECESSARY OFF SWITCHES
 
 450 GO TO ihop, (451,170)
 451 IF (funlg == fculg) RETURN
 lmt2 = lmt1 + 1
 lmt3 = fculg*entn1
 DO  i = lmt2,lmt3,entn1
   IF (fequ(i) < 0) CYCLE
   trial = andf(rmsk,FILE(i))
   IF (trial == rmsk) CYCLE
   ifordi = andf(lmsk,ford(i))
   DO  j = 1,lmt3,entn1
     IF (trial /= andf(rmsk,FILE(j))) CYCLE
     IF (i == j) CYCLE
     IF (andf(lmsk,ford(j))-ifordi < 0.0) THEN
       GO TO   452
     ELSE IF (andf(lmsk,ford(j))-ifordi == 0.0) THEN
       GO TO   460
     ELSE
       GO TO   454
     END IF
     452 fon(j) = orf(s,fon(j))
     CYCLE
     454 fon(i) = orf(s,fon(i))
     460 CONTINUE
   END DO
 END DO
 RETURN
 
 900 WRITE  (outtap,901) sfm
 901 FORMAT (a25,' 1021, FIAT OVERFLOW.')
 CALL mesage (-37,0,nclean)
 RETURN
END SUBROUTINE xclean
