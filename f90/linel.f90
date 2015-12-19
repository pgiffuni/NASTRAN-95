SUBROUTINE linel (iz,nwds,opcor,opt,x,pen,deform,gplst)
     
!     CALL TO LINEL IS AS FOLLOWS -
 
!     (1)
!     OPT = ZERO (INPUT) - TO CREATE COMPLETE LINE CONNECTION TABLE OF
!     **********           ELEMENTS OF ALL TYPES, TO BE USED BY SUPLT
!                          SUBROUTINE
!        INPUT-
!           OPCOR (INPUT) = NUMBER OF WORDS OF OPEN CORE FOR -IZ-
!        OUTPUT-
!           IZ   = LIST OF GRID POINT ELEMENET CONNECTIONS AND POINTERS
!                  TO EACH GRID POINT, FROM IZ(1) THRU IZ(NWDS). DATA
!                  COMPOSED OF   1. GPCT,  AND 2. NGP WORDS OF CONTROL
!                  POINTERS
!           NWDS = NO. OF WORDS IN IZ PRIOR TO POINTER ARRAY.
!                  I.E. 1 LESS THAN LOCATION OF POINTERS,
!                = 0 IF ARRAY NOT CREATED
!           OPT  = NWDS
 
!     (2)
!     OPT = NONZERO (INPUT) - LOAD INTO CORE THE GRID POINT CNNECTION
!     *************           LIST OF ALL ELEMENTS OF THE SAME TYPE
 
!        INPUT-
!           NWDS  = ETYP, 2 BCD WORDS (CALLING ROUTINE HAS ALREADY READ
!                   THIS WORD FROM DATA BLOCK ELSET)
!           OPT   = MO. OF GRID POINT CONNECTIONS PER ELEMENT, NGPEL
!                   (CALLING ROUTINE HAS ALREADY READ THIS WORD)
!           OPCOR = OPEN CORE AVAILABLE W.R.T. IZ(1)
!           GPLST = A SUBSET LIST OF GRID POINTS PERTAINING TO THOSE
!                   POINTS USED ONLY IN THIS PLOT
!        OUTPUT-
!           IZ    = GRID POINT CONNECTION LIST FOR ALL ELEMENTS OF THIS
!                   TYPE, OR AS MANY ELEMS OF THIS TYPE AS CORE ALLOWS.
!           NWDS  = TOTAL LENGTH OF TABLE IZ
!           OPT   = NUMBER OF CONNECTIONS PER ELEMENT
!           (IF INSUFF. CORE TO READ ALL THE ELEMENTS, BOTH NWDS AND OPT
!           ARE SET TO NEGATIVE UPON RETURN. FURTHER CALLS MUST BE MADE
!           TO COMPLETE THIS ELEMENT
!           IF ILLEGAL ELEMENT IS ENCOUNTERED, NWDS AND OPT ARE SET TO
!           ZERO, AND ELSET IS SPACED OVER THE ELEMENT)
 
!           (NOTE THAT  'DO 100 I=1,NWDS,OPT'  MAY THEN BE USED
!           BUT IT IS MORE EFFICIENT TO USE  'DO 100 I=1,NWDS' AND CHECK
!           ZERO AS THE COMMAND TO LIFT THE PEN)
 
!     EACH ELEMENT TYPE HAS THE FOLLOWING DATA IN ELSET FILE
!           ELTYP = BCD SYMBOL (1 WORD)
!           NGPEL = NUM. GRID POINTS.
!                   IF NEGATIVE OR .GT. 4 NOT A CLOSED LOOP
!           ELID  = ELEMENT ID
!           G     = NGPEL GRIDS.
!           LOOP THRU ELID AND G UNTIL ELID = 0 (I.E. NO MORE ELEMS OF
!                                              THIS TYPE)
!     (3)
!     ELEMENT OFFSET PLOT (UNDEFORMED PLOT ONLY, PEDGE=3),
!     *******************
!     IF ELEMENTS WITH OFFSET ARE PRESENT, CALL OFSPLT TO PLOT THEM OUT
!     AND DO NOT INCLUDE THEM IN THE IZ TABLE
!     IF OFFSET COMMAND IS REQUESTED BY USER VIA THE PLOT CARD
!     (PEDGE = 3), SKIP COMPLETELY THE GENERATION OF THE IZ TABLE
 
!     OFFSET n OPTION (ON PLOT CONNAND CARD IN CASE CONTROL SECTION) -
!       n .LT. 0, SKIP OFFSET VALUES ON GENERAL PLOTS. (PEDGE.NE.3)
!       n =    0, OFFSET VALUES INCLUDED IN ALL GENERAL PLOTS (PEDGE=3)
!       n .GT. 0, PLOT ONLY THOSE ELEMENTS HAVING OFFSET DATA, OFFSET
!                 DATA ARE MAGNIFIED n TIMES. (PEDGE=3)
!     SUBROUTINE PLOT SETS THE PEDGE FLAG, AND PLTSET SETS THE OFFSCL.
 
 
 INTEGER, INTENT(OUT)                     :: iz(1)
 INTEGER, INTENT(IN OUT)                  :: nwds
 INTEGER, INTENT(IN OUT)                  :: opcor
 INTEGER, INTENT(IN OUT)                  :: opt
 REAL, INTENT(IN OUT)                     :: x(3,1)
 INTEGER, INTENT(IN OUT)                  :: pen
 INTEGER, INTENT(IN OUT)                  :: deform
 INTEGER, INTENT(IN OUT)                  :: gplst(1)
 INTEGER :: elid,elset,etyp,g, m1(16),NAME(2),  &
     ng(121), TYPE,ngtyp(2,13),ldx(9),offscl, offset, pedge
 
 COMMON /BLANK / ngp,skp1(9),skp2(2),elset,skp3(7),merr
 COMMON /system/ skp4,iout
 COMMON /pltscr/ nnn,g(3)
 COMMON /drwdat/ skp5(15),pedge
 COMMON /xxparm/ skp6(235),offscl
 DATA    NAME  / 4HLINE, 1HL  /,  nm1,m1 / 16,  &
     4H(33X, 4H,13H, 4HELEM, 4HENT , 4HTYPE, 4H ,a5,  &
     4H,4HW, 4HITH,, 4HI8,2, 4H4H g, 4HRIDS, 4H ski,  &
     4HPPED, 4H in , 4HLINE, 4HL.)   /
 
!     SPECIAL ELEMENT CONNECTION PATTERNS
 
 DATA ldx  / 2HD1,2HD2,2HD3,2HD4,2HD5,2HD6,2HD7,2HD8,2HD9      /
 DATA ktet / 2HTE /, kweg / 2HWG /, khx1 / 2HH1 /, khx2 / 2HH2 /,  &
     kix1 / 2HXL /, kix2 / 2HXQ /, kix3 / 2HXC /, kae  / 2HAE /,  &
     ktm6 / 2HT6 /,ktrplt/ 2HP6 /,ktrshl/ 2HSL /, kfh1 / 2HFA /,  &
     kfh2 / 2HFB /, kfwd / 2HFW /, kfte / 2HFT /, k2d8 / 2HD8 /,  &
     khb  / 2HHB /, kbar / 2HBR /, kt3  / 2HT3 /, kq4  / 2HQ4 /
 
!     NGTYP(1,TYPE) = LOCATION WORD 1 IN -NG-, +N = POINTER TO G
!                                              -N = THRU POINTER TO G
!     BE SURE TO KEEP PEN DOWN                  0 = LIFT PEN
!     AS MUCH AS POSSIBLE.
!     NGTYP(2,TYPE) = NUMBER OF ENTRIES/ELEMENT MINUS 1 IN TABLE IZ
 
 DATA ngtyp/ 0,0,  3,9,  10,14,  22,19,  37,30,  56,43,  79,6,  &
     83,7, 86,9,  95,10, 102, 8, 108, 2, 110, 7/
 DATA  ng  / &
!    1 - LINE,TRIANGLE,QUAD  &
 1,-5, &
!    2 - TETRA (WORD 3)  &
 1,-4,1,3,0,2,4, &
!    3 - WEDGE (WORD 10)  &
 1,-3,1,4,-6,4,0,5,2,0,3,6, &
!    4 - HEXA  (WORD 22)  & 
 1,-4,1,5,-8,5,0,6,2,0,3,7,0,8,4, &
!    5 - IHEXA2 (WORD 37)  & 
 1,-8,1,9,13,-20,13,0,15,10,3,0,5,11,17,0,19,12,7, &
!    6 - IHEXA3 (WORD 56)  &
 1,-12,1,13,17,21,-32,21,0,24,18,14,4,0,7,15,19,27,0,30,20,16,10 &
!    7 - AREO (WORD 79)  &
 ,  1,-4,1,0, &
!    8 - TRIM6, TRPLT1, AND TRSHL (WORD 83)  &
 1,-6,1, &
!    9 - IS2D8 (WORD 86)  &
 1,5,2,6,3,7,4,8,1, &
!   1O - POINT (WORD 95)  &
 2,-6,7,2,0,1,8, &
!   11 - LINE (WORD 102)  &
 3,-6,3,0,7,8, &
!   12 - REV OR ELIP CYL. (WORD 108)  &
 1,2, &
!   13 - AREA3 (WORD 110)  &
 1,-3,1,0,4,5, &
!   14 - AREA4 (WORD 116)  &
 1,-4,1,0,5,6 /
 
 k    = 1
 IF (opt == 0) GO TO 20
 etyp = nwds
 i    = opt
 GO TO 30
 
 20 IF (opt /= 0) GO TO 170
 CALL READ (*420,*190,elset,etyp,1,0,i)
 CALL fread (elset,i,1,0)
 
 30 ngpel  = IABS(i)
 ngpelx = ngpel
 offset = 0
 IF (etyp == kbar) offset = 6
 IF (etyp == kt3 .OR. etyp == kq4) offset = 1
 
 TYPE  = 1
 IF (etyp == ktet .OR. etyp == kfte) TYPE = 2
 IF (etyp == kweg .OR. etyp == kfwd) TYPE = 3
 IF (etyp == khx1 .OR. etyp == khx2 .OR. etyp == kfh1 .OR.  &
     etyp == kfh2 .OR. etyp == kix1) TYPE = 4
 IF (etyp == kix2) TYPE = 5
 IF (etyp == kix3) TYPE = 6
 IF (etyp ==  kae) TYPE = 7
 IF (etyp == ktm6 .OR. etyp == ktrplt .OR. etyp == ktrshl) TYPE = 8
 IF (etyp == k2d8) TYPE = 9
 IF (etyp == khb ) TYPE = 10
!                   CHBDY TYPE = 10,11,12,13,14
 
 IF (TYPE /= 1) GO TO 40
 
!     SIMPLE ELEMENT
 
 IF (ngpel > 2 .AND. i > 0) ngpelx = ngpel + 1
 IF (ngpel > 4) GO TO 131
 l1 = 1
 m  = ngpelx
 GO TO 50
 
!     COMPLEX ELEMENT
 
 40 l1 = ngtyp(1,TYPE)
 m  = ngtyp(2,TYPE)
 50 IF (ngpelx > nnn) GO TO 140
 
!     READ THE ELEMENT DATA
 
 55 CALL fread (elset,elid,1,0)
 IF (elid <= 0) GO TO 20
 CALL fread (elset,lid,1,0)
 CALL fread (elset,g,ngpel,0)
 IF (ngpel /= ngpelx) g(ngpelx) = g(1)
 
!     CALL OFSPLT TO PROCESS OFFSET PLOT
 
 IF (offset /= 0) CALL ofsplt (*55,etyp,elid,g,offset,x,deform,gplst)
 IF (TYPE < 10 .OR. TYPE > 14) GO TO 57
 
!     SPECIAL HANDLING FOR CHBDY
 
 TYPE = 9 + g(ngpel)
 l1 = ngtyp(1,TYPE)
 m  = ngtyp(2,TYPE)
 
 57 l = l1
 
 IF (opt /= 0) GO TO 70
 
!     CREATING CONNECTION ARRAY FOR SUPLT
 
 ll = 0
 i1 = 0
 60 i2 = ng(l)
 IF (i1 == 0) GO TO 66
 IF (i2 < 0) THEN
   GO TO    62
 ELSE IF (i2 == 0) THEN
   GO TO    64
 ELSE
   GO TO    65
 END IF
 
!     THRU RANGE
 
 62 i2 =-i2
 i2 = MIN0(i2,m)
 j  = i1 + 1
 i1 = g(i1)
 IF (2*(i2-j+1)+k > opcor) GO TO 390
 DO  i = j,i2
   iz(k  ) = MIN0(g(i),i1)
   iz(k+1) = MAX0(g(i),i1)
   k  = k  + 2
   ll = ll + 1
   i1 = g(i)
 END DO
 IF (ll == m-1) ll = ll - 1
 GO TO 66
 
 64 i1 = 0
 l  = l + 1
 GO TO 60
 
 65 IF (k+1 > opcor) GO TO 180
 iz(k  ) = MIN0(g(i2),g(i1))
 iz(k+1) = MAX0(g(i2),g(i1))
 k  = k  + 2
 66 ll = ll + 1
 i1 = i2
 IF (ll >= m) GO TO 55
 l  = l + 1
 GO TO 60
 
!     ON CONVERSION REMOVE ABOVE CODE
 
!     LOAD ELEMENT INTO CORE
 
 70 n = k + m
 
!     THIS TEST PROTECTS THE CORE FOR THE FIRST ELEMENT READ
 
 IF (n+1 > opcor) GO TO 140
 i1 = 0
 i2 = ng(l)
 GO TO 125
 80 IF (i1 == 0) GO TO 90
 IF (i2 < 0) THEN
   GO TO   110
 ELSE IF (i2 == 0) THEN
   GO TO   100
 END IF
 
 90 iz(k) = g(i2)
 GO TO 120
 100 iz(k) = i2
 GO TO 120
 110 i2 =-i2
 
!     NEXT LINE FOR ELEMENTS WITH MORE THAN ONE THRU POINTER
 
 IF (n /= k+m) i1 = i1 + 1
 DO  i = i1,i2
   iz(k) = g(i)
   k = k + 1
 END DO
 k = k - 1
 120 k = k + 1
 IF (k >= n) GO TO 130
 125 i1 = i2
 l  = l + 1
 i2 = ng(l)
 GO TO 80
 
!     STORE ZERO AT THE END OF EACH ELEMENT
 
 130 iz(k) = 0
 k = k + 1
 IF (k+m+1 > opcor) GO TO 180
 GO TO 55
 
!     CHECK FOR PDUM ELEMENTS BEFORE REJECTING
 
 131 DO  ii = 1,9
   IF (etyp == ldx(ii)) CALL pdumi (*20,*180,*140,ii,m,opcor,ngpel,  &
       k,elset,opt)
 END DO
 
!     ILLEGAL ELEMENT, NO CORE FOR 1 ELEMENT
 
 140 g(1) = 2
 g(2) = etyp
 g(3) = ngpel
 CALL wrtprt (merr,g,m1,nm1)
 
!     READ TO THE END OF THIS ELEMENT
 
 150 CALL fread (elset,elid,1,0)
 IF (elid <= 0) GO TO 160
 j = 1 + ngpel + offset
 CALL fread (elset,0,-j,0)
 GO TO 150
 160 CONTINUE
 
!     NOTE THAT BOTH OPT AND NWDS=0 FOR ILLEGAL ELEMENTS
 
 IF (opt /= 0) GO TO 390
 GO TO 20
 
!     END OF OPT.NE.0
 
 170 nwds = k - 1
 opt  = m + 2
 GO TO 410
 
!     INSUFFICIENT CORE FOR ALL ELEMENTS
 
 180 IF (opt == 0) GO TO 390
 nwds = 1 - k
 opt  = -(m+2)
 GO TO 410
 
!     SORT
 
 190 IF (pedge == 3) GO TO 400
 IF (opt   /= 0) GO TO 170
 IF (k     <= 1) GO TO 400
 CALL sort (0,0,2,1,iz,k-1)
 
!     NWDS IS SET TO NO. OF WORDS PRIOR TO ELIMINATING DUPLICATES
 
 nwds = k - 1
 IF (nwds <= 2) GO TO 310
 ASSIGN 310 TO iret
 
!     ELIMINATE DUPLICATE ENTRIES FROM LIST SORTED ON FIRST ENTRY
 
 200 CONTINUE
 i = 1
 l = 1
 ll= iz(l)
 
 
 DO  j = 3,nwds,2
   IF (iz(j) == ll) GO TO 220
   
!     NEW PIVOT
   
   l  = i + 2
   ll = iz(j)
   GO TO 230
   220 IF (iz(j+1)-iz(i+1) < 0) THEN
     GO TO   240
   ELSE IF (iz(j+1)-iz(i+1) == 0) THEN
     GO TO   300
   END IF
   
!     UNIQUE ENTRY FOR PIVOT FOUND
   
   230 iz(i+2) = ll
   iz(i+3) = iz(j+1)
   GO TO 290
   
!     SECOND COLUMN OUT-OF-SORT
!     LOAD ENTRY SORTED.  CHECK PREVIOUS ENTRIES
!     L = LOWER LIMIT OF COLUMN 1 FOR MERGING
!     K SET TO FIRST ENTRY OF NEXT NEW ENTRY IN LIST INITIALLY
   
   240 k = i
   250 IF (k <= l) GO TO 270
   IF (iz(j+1)-iz(k-1) < 0) THEN
     GO TO   260
   ELSE IF (iz(j+1)-iz(k-1) == 0) THEN
     GO TO   300
   ELSE
     GO TO   270
   END IF
   260 k = k - 2
   GO TO 250
   
!     LOAD ENTRY INTO LOCATION
   
   270 n = iz(j+1)
   m = i + 2
   280 iz(m+1) = iz(m-1)
   m = m - 2
   IF (m > k) GO TO 280
   iz(k+1) = n
   iz(i+2) = ll
   
!     INCREMENT FOR ENTRY LOADED
   
   290 i = i + 2
   300 CONTINUE
 END DO
 
!     NWDS RESET TO NO. WORDS AFTER ELIMINATING DUPLICATE ENTRIES
 
 nwds = i + 1
 GO TO iret, (310,330)
 
 
!     K IS SET TO THE NEXT PART OF CORE WHICH WILL BE FILLED WITH THE
!     HIGHER ENTRY IN THE FIRST POSITION
 
 310 k = nwds + 1
 IF (2*nwds > opcor) GO TO 400
 DO  i = 1,nwds,2
   iz(k  ) = iz(i+1)
   iz(k+1) = iz(i  )
   k = k + 2
 END DO
 nwds = k - 1
 CALL sort (0,0,2,1,iz,nwds)
 ASSIGN 330 TO iret
 GO TO 200
 
 330 CONTINUE
 IF (nwds+ngp+1 > opcor) GO TO 400
 k = 1
 j = 1
 l = 1
 m = 1 + nwds
 i = 0
 iz(m) = 1
 
!     CREATE A GPCT --- M = POINTER FOR POINTER ARRAY
!                       L = SIL NUMBER
!                       J = POINTER TO NEXT GPCT ENTRY
 
 340 IF (iz(k) == l) GO TO 360
 
!     NEW PIVOT
 
 350 m = m + 1
 iz(m) = iz(m-1) + i
 l = l + 1
 i = 0
 IF (l >  ngp) GO TO 370
 IF (k > nwds) GO TO 350
 GO TO 340
 
!     CONNECTED POINT
 
 360 iz(j) = iz(k+1)
 k = k + 2
 j = j + 1
 i = i + 1
 GO TO 340
 
!     EFFICIENCY PLOT POSSIBLE
 
 370 CONTINUE
 opt = nwds
 GO TO 410
 
 390 opt  = 0
 400 nwds = 0
 410 RETURN
 
 420 CALL mesage (-2,elset,NAME)
 GO TO 410
 
END SUBROUTINE linel
