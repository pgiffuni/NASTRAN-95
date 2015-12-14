SUBROUTINE sma3a (mcbcur)
!*****
! THIS ROUTINE BUILDS A GENERAL ELEMENT MATRIX (DOUBLE PRECISION AND
! SYMMETRIC) OF SIZE LUSET X LUSET.  MCBCUR IS THE MATRIX CONTROL BLOCK
! FOR THIS MATRIX.
!*****
 
 INTEGER, INTENT(IN OUT)                  :: mcbcur(7)
 DOUBLE PRECISION :: dq(1)              ,dpword  &
     ,                  det
 INTEGER :: iq(1)              ,eor  &
     ,                  outrw              ,clsrw  &
     ,                  clsnrw             ,small
 LOGICAL :: zonly
 
 DIMENSION  mcb(7)  &
     ,                  q(1)               ,ibuff3(3) ,                  NAME(2)
 
 COMMON /BLANK/ luset              ,ngenel  &
     ,                  noecpt
 COMMON   /system/ isys               ,dummy(53)  &
     ,                  iprec
 COMMON   /zzzzzz/ q
 COMMON   /genely/ ifgei              ,ifkggx  &
     ,                  ifout              ,ifa  &
     ,                  ifb                ,ifc  &
     ,                  ifd                ,ife ,                  iff  &
     ,                  inrw               ,outrw  &
     ,                  clsrw              ,clsnrw  &
     ,                  eor                ,neor  &
     ,                  mcba(7)            ,mcbb(7)  &
     ,                  mcbc(7)            ,mcbd(7)  &
     ,                  mcbe(7)            ,mcbf(7)  &
     ,                  mcbkgg(7) ,                  iui                ,iud  &
     ,                  izi                ,is  &
     ,                  izis               ,istzis  &
     ,                  ibuff3             ,left
 COMMON   /zblpkx/ dpword             ,dum2(2)  &
     ,                  INDEX
 
 EQUIVALENCE (iq(1),dq(1),q(1))  &
     ,                  (ibuff3(2),m)      ,(ibuff3(3),n)
 
 DATA               NAME(1)/4HSMA3/    ,NAME(2)/4HA   /
 
! MAKE THE ARGUMENT A LOCAL VARIABLE
 
 DO  i=1,7
   mcb(i) = mcbcur(i)
 END DO
 
! READ THE UI SET OF SCALAR INDEX NUMBERS INTO OPEN CORE.
 
 CALL fread(ifgei,iq(iui+1),m,0)
 
! IUD POINTS TO THE ZEROTH LOCATION OF THE UD ARRAY.
 
 iud = iui + m
 left = left - m
 
! SET UP ARITHMETIC CONSTANTS.
 
 mpn = m + n
 msq = m**2
 nsq = n**2
 zonly = .false.
 IF (n == 0) zonly = .true.
 IF (zonly) GO TO 20
 
! SINCE N .NE. 0, THE UD SET EXISTS.  READ IT INTO CORE.
 
 CALL fread(ifgei,iq(iud+1),n,0)
 left = left - n
 
! BUILD THE ARRAY IQ(IP+1),IQ(IP+2),...,IQ(IP+MPN) SUCH THAT
! IQ(IP+K) = L IMPLIES IQ(IUI+L) IS THE K TH SMALLEST NUMBER OF THE
! SET OF NUMBERS IQ(IUI+1),...,IQ(IUD+N)
 
 20 ip = iui + mpn
 k  = ip
 limk = ip + mpn
 low = iui + 2
 lim = iui + mpn
 30 small = iq(iui+1)
 ismall = iui + 1
 DO  j=low,lim
   IF (iq(j) >= small) CYCLE
   small  = iq(j)
   ismall = j
 END DO
 k = k + 1
 idiff = ismall - iui
 iq(k) = idiff
 iq(idiff) = iq(idiff) + luset
 IF (k < limk) GO TO 30
 low = iui + 1
 DO  i=low,lim
   IF (iq(i) <= luset) CALL mesage (-30,28,5)
   iq(i) = iq(i) - luset
 END DO
 
! READ INDICATOR OF Z OR K MATRIX
 
 CALL fread(ifgei,izk,1,0)
 
! SET UP POINTERS TO THE ZEROTH LOCATION OF THE DOUBLE PRECISION ARRAYS
!       -1
! K  ORZ  AND S
!  E    E      E
 
 izi = (iui + 2*mpn - 1) / 2  +  2
 is = izi + msq
 
! READ IN THE M**2 SINGLE PRECISION ELEMENTS OF THE SYMMETRIC Z OR K
! INTO A TEMPORARY BUFFER BEGINNING AT Q(IBUFF)
 
 ibuff = iui + 2 * (mpn + msq)
 
! IF ALL OF Z OR K CANNOT FIT INTO THIS BUFFER, READ BLOCKS OF M WORDS
 
 IF (ibuff + msq > left) GO TO 70
 ind = neor
 IF (zonly) ind = eor
 CALL fread(ifgei,iq(ibuff+1),msq,ind)
 
! STORE THE SINGLE PRECISION MATRIX IN ITS DOUBLE PRECISION LOCATION.
 
 lim = izi + msq
 i   = izi
 j   = ibuff
 60 i   = i + 1
 IF (i > lim) GO TO 100
 j   = j + 1
 dq(i) = q(j)
 GO TO 60
 
! READ Z OR K INTO THE BUFFER M WORDS AT A TIME AND STORE M WORDS
! AT A TIME
 
 70 ind = neor
 DO  k=1,m
   IF (k == m  .AND.  zonly) ind = eor
   CALL fread(ifgei,q(ibuff+1),m,ind)
   i = izi + (k - 1) * m
   j = ibuff
   lim = i + m
   80 i = i + 1
   IF (i > lim) CYCLE
   j = j + 1
   dq(i) = q(j)
   GO TO 80
 END DO
 
! IF K IS INPUT DO NOT COMPUTE INVERSE
 
 100 IF (izk == 2) GO TO 105
!*****
! COMPUTE THE INVERSE OF Z
!                        E
!*****
 
! THE 4TH ARGUMENT OF INVERD IS A DUMMY D.P. ARGUMENT WHILE 3 * M
! WORDS OF WORKING STORAGE ARE NEEDED FOR THE 8TH ARGUMENT OF SUBROUTINE
! INVERD.  SUBROUTINE INVERD WILL RETURN Z  INVERSE AT DQ(IZI+1)
!                                         E
 
 ibuff = iui + 2 * (mpn + msq) + 5
 ii = ibuff + 2 * m
!     NO NEED TO COMPUTE DETERMINANT SINCE IT IS NOT USED SUBSEQUENTLY.
 ising = -1
 CALL inverd (m,dq(izi+1),m,iq(ibuff+1),0,det,ising,iq(ii+1))
 
! ISING = 2 IMPLIES A SINGULAR Z
!                               E
 
 IF (ising == 2) CALL mesage (-5,idummy,NAME)
 105 CONTINUE
 
! READ IN THE M*N ELEMENTS OF THE M X N  S MATRIX IF N .GT. 0.
! THIS MATRIX IS SINGLE PRECISION AND ROW STORED.
 
 IF (zonly) GO TO 130
 ibuff = mpn + 2*msq + 2*m*n +5
 CALL fread(ifgei,q(ibuff+1),m*n,1)
 
! STORE THE S  MATRIX AT DQ(IS+1) MAKING S  DOUBLE PRECISION
!            E                            E
 
 low = is + 1
 lim = is + m*n
 j = ibuff
 DO  i=low,lim
   j = j + 1
   dq(i) = q(j)
 END DO
!                  -1
! COMPUTE K S  OR Z  S AND STORE AT DQ(IZIS+1)
!          E E     E  E
 
 izis = is + m*n
 CALL gmmatd (dq(izi+1),m,m,0, dq(is+1),m,n,0, dq(izis+1) )
 
!          T        T -1
! COMPUTE S K S OR S Z  S AND STORE AT DQ(ISTZIS+1)
!          E E E    E E  E
 
 istzis = izis + m*n
 CALL gmmatd (dq(is+1),m,n,1, dq(izis+1),m,n,0, dq(istzis+1) )
 
!              -1
! SET K S  OR Z  S  NEGATIVE
!      E E     E  E
 
 low = izis + 1
 lim = izis + m*n
 DO  i=low,lim
   dq(i) = -dq(i)
 END DO
!*****
! AT THIS POINT ALL MATRICES HAVE BEEN COMPUTED
!*****
 
! INITIALIZE FOR OUTPUT ONTO THE FILE
 
 130 izrow  = 1
 izscol = 1
 icol = 1
 limjui = iui + m
 limjud = iui + mpn
 jui  = iui + 1
 jud  = iud + 1
!******
! BEGIN OUTPUT LOOP
!******
 140 iloop = 0
 IF (zonly) GO TO 150
 IF (iq(jui) - iq(jud) < 0) THEN
   GO TO   150
 ELSE IF (iq(jui) - iq(jud) == 0) THEN
   GO TO   470
 ELSE
   GO TO   240
 END IF
 
! AT THIS POINT, WRITE OUT COLUMN(S) CORRESPONDING TO THE UI SET.
 
 150 IF (iq(jui) - icol < 0) THEN
   GO TO   160
 ELSE IF (iq(jui) - icol == 0) THEN
   GO TO   190
 ELSE
   GO TO   170
 END IF
 
! A TRANSFER TO STATEMENT NO. 1115 WILL BE MADE IF THE MAXIMUM OF THE
! UD SET IS LESS THAN THE MINIMUM OF THE UI SET AND THE COLUMNS
! CORRESPONDING TO THE UD SET HAVE BEEN OUTPUT.
 
 160 IF (iloop == 1  .OR.  zonly) GO TO 480
 iloop = 1
 GO TO 240
 
! SINCE IQ(JUI) .GT. ICOL, IQ(JUI) - ICOL COLUMNS OF ZERO VECTORS MUST
! BE OUTPUT.
 
 170 lim = iq(jui) - icol
 DO  i=1,lim
   CALL bldpk (2, iprec, mcb(1), 0, 0)
   CALL bldpkn(mcb(1),0,mcb)
 END DO
 
! INITIALIZE FOR THE OUTPUT OF THE CURRENT COLUMN BY CALLING BLDPK
 
 190 CALL bldpk (2, iprec, mcb(1), 0, 0)
 DO  i=1,mpn
   ippi = ip + i
   IF (iq(ippi) > m) GO TO 200
   
! SINCE IQ(IPPI).LE.M,OUTPUT AN ELEMENT OF K OR Z INVERSE
   
   jrow = izrow
   jcol = iq(ippi)
   k = (jrow - 1) * m  +  jcol  +  izi
   GO TO 210
   
! HERE WE ARE DEALING WITH A MEMBER OF THE UD SET.  HENCE AN ELEMENT OF
!                -1
! THE -K S  OR -Z  S  MATRIX MUST BE OUTPUT
!       E E      E  E
   
   200 jrow = izrow
   jcol = iq(ippi) - m
   k = (jrow - 1) * n  +  jcol  +  izis
   
! FILL ZBLPKI COMMON BLOCK
   
   210 kk = iq(ippi)
   INDEX = iq(kk)
   dpword = dq(k)
   IF (dpword /= 0.0D0) CALL zblpki
 END DO
 
! THE CURRENT COLUMN IS COMPLETE.  CALL BLDPKN TO WRAP UP.
 
 CALL bldpkn(mcb(1),0,mcb)
 izrow = izrow + 1
 icol  = iq(jui) + 1
 jui = jui + 1
 IF (jui > limjui) jui = limjui
 230 IF (izrow > m  .AND.  izscol > n) GO TO 320
 GO TO 140
 
! AT THIS POINT WRITE OUT A COLUMN(S) USING THE UD SET.
 
 240 IF (iq(jud) - icol < 0) THEN
   GO TO   250
 ELSE IF (iq(jud) - icol == 0) THEN
   GO TO   280
 ELSE
   GO TO   260
 END IF
 
! A TRANSFER TO STATEMENT NO. 1185 WILL BE MADE IF THE MAXIMUM OF THE
! UI SET IS LESS THAN THE MINIMUM OF THE UD SET AND THE COLUMNS
! CORRESPONDING TO THE UI SET HAVE BEEN OUTPUT.
 
 250 IF (iloop == 1) GO TO 490
 iloop = 1
 GO TO 150
 
! WRITE ZERO COLUMN(S).
 
 260 lim = iq(jud) - icol
 DO  i=1,lim
   CALL bldpk (2, iprec, mcb(1), 0, 0)
   CALL bldpkn(mcb(1),0,mcb)
 END DO
 280 CALL bldpk (2, iprec, mcb(1), 0, 0)
 
! OUTPUT A COLUMN WHOSE SIL NO. IS A MEMBER OF THE UD SET.
 
 DO  i=1,mpn
   ippi = ip + i
   IF (iq(ippi) > m) GO TO 290
   
!                                           -1
! SINCE IQ(IPPI).LE.M,AN ELEMENT OF -KS OR -Z  S MUST BE OUTPUT
   
   jrow = iq(ippi)
   jcol = izscol
   k = (jrow - 1) * n  +  jcol  +  izis
   GO TO 300
   
!                       T         T -1
! OUTPUT AN ELEMENT OF S K S  OR S Z  S
!                       E E E     E E  E
   
   290 jrow = iq(ippi) - m
   jcol = izscol
   k = (jrow - 1) * n  + jcol  +  istzis
   
! SET UP PARAMETERS IN ZBLPKI COMMON BLOCK
   
   300 kk = iq(ippi)
   INDEX = iq(kk)
   dpword = dq(k)
   IF (dpword /= 0.0D0) CALL zblpki
 END DO
 
! WRAP UP THIS COLUMN.
 
 CALL bldpkn(mcb(1),0,mcb)
 izscol = izscol + 1
 icol   = iq(jud)+ 1
 jud    = jud + 1
 IF (jud > limjud) jud = limjud
 GO TO 230
 
! DETERMINE IF ZERO COLUMNS ARE TO BE OUTPUT.
 
 320 k = iui + m
 l = iud + n
 MAX = iq(k)
 IF (iq(l) > MAX) MAX = iq(l)
 lim = MAX - luset
 IF (lim < 0) THEN
   GO TO   330
 ELSE IF (lim == 0) THEN
   GO TO   350
 ELSE
   GO TO   500
 END IF
 
! OUTPUT LIM ZERO COLUMNS
 
 330 lim = IABS(lim)
 DO  i = 1,lim
   CALL bldpk (2, iprec, mcb(1), 0, 0)
   CALL bldpkn(mcb(1),0,mcb)
 END DO
 350 DO  i=1,7
   mcbcur(i) = mcb(i)
 END DO
 RETURN
 
! FATAL ERROR MESSAGES
 
 470 CALL mesage (-30,28,1)
 480 CALL mesage (-30,28,2)
 490 CALL mesage (-30,28,3)
 500 CALL mesage (-30,28,4)
 RETURN
END SUBROUTINE sma3a
