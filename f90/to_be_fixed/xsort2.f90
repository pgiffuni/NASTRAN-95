SUBROUTINE xsort2
     
!     XSORT2 REPLACES XSORT FOR SPEED AND EFFICIENCY
 
!     XSORT2 REQUIRES IFP MODULE TO USE RCARD2 ROUTINE INSTEAD OF
!     RCARD (DUE TO THE ASTERISK POSITION IN DOUBLE FIELD INPUT
!     CARD HAS NOT BEEN MOVED TO COLUMN 8)
 
!     XSORT2 READS BULKDATA CARDS FROM THE INPUT TAPE, ADJUSTS THE
!     FIELDS, PERFORMS AN ALPHA-NUMERIC SORT ON THE CARD IMAGES FROM
!     LEFT TO RIGHT, INSERTS CONTINUATION CARDS IN THEIR PROPER
!     POSITION, AND PLACES THE RESULTING SORTED IMAGES ON THE NEW
!     PROBLEM TAPE, NPTP.
 
!     THIS ROUTINE DOES NOT USE XRECPS, RPAGE, INITCO, XFADJ, XFADJ1,
!     XBCDBI, XPRETY, EXTINT, INTEXT, CRDFLG, ISFT, AND THE CHARACTER
!     FUNCTIONS KHRFNi.
!     IT CALLS ONLY SORT2K - TO SORT IN-CORE DATA USING TWO SORT KEYS
!              AND  BISLC2 - BINARY SEARCH USING TWO SORTED KEYS
 
!     XSORT2 NEW LOGIC -
 
!     1.  INPUT BULKDATA CARDS ARE READ INTO OPEN CORE, EXCEPT CONTINU-
!         ATION (* OR +), DELETE (/), COMMENT ($), AND BLANK CARDS.
!     2.  WHEN CORE IS FULL, OR LAST INPUT DATA READ, SORT DATA IN CORE
!         AND WRITE THE ENTIRE SORTED DATA TO SEQUENTIAL GINO FILE 303.
!     3.  REPEAT 1 AND 2, AND WRITE DATA TO GINO FILES 304,305,306 ETC.
!         IF NECESSARY. UP TO 30 FILES ARE ALLOWED.
!     4.  ALL CONTINUATION CARDS ARE WRITEN TO GINO FILE 302. ALL
!         DELETES TO 301. BLANK AND COMMENT CARDS ARE IGNORED.
!     5.  WHEN ALL INPUT DATA CARDS ARE READ AND SAVED IN GINO FILE(S),
!         RE-LOAD THE DELETE CARDS FROM 301 INTO OPEN CORE SPACE, AND
!         COPY OPTP TO 301 WITH DESIGNATED CARDS DELETED.
!     6.  COMPUTE BUFFER SPACE (AT THE END OF OPEN CORE) AND THE WORK
!         SPACE (AT THE BEGINNING OF OPEN CORE) NEEDED FOR FILE MERGE
!         OPERATION, AND READ INTO CORE ALL CONTINUATION CARDS USING
!         THE REMAINING CORE SPACE.
!     7.  IF CORE SPACE IS NOT BIG ENOUGH TO HOLD ALL CONTINUATION
!         CARDS, CREATE A CONTINUATION-INDEX TABLE IN CORE, AND MOVE THE
!         CONTINUATION CARDS TO A NEW GINO FILE, WITH LARGE BLOCKS OF
!         CONTINUATION CARDS
!     8.  PRE-MERGE BULKDATA GINO FILES TO SAVE BUFFER SPACE IF MORE
!         THAN 9 GINO FILES WERE USED IN STEP 3.
!         PERFORM A 2-TO-1 MERGE IF 10 TO 17 FILES WERE INVOLVED, OR
!         A 3-TO-1 MERGE IF MORE THAN 17 FILES WERE USED IN STEP 3.
!         THE MERGE FILES ARE SAVED IN 302,303,304,305 ETC.
!     9.  MERGE ALL FILES IN SORTED ORDER, AND INSERT CONTINUATION CARDS
!         WHEN NECESSARY. THE MERGED RESULTS ARE WRITTEN TO NPTP
!     10. ECHO ANY CONTINUATION CARD WHICH HAS NO PARENT AND THEREFORE
!         NOT USED. MAKE SURE NO REDUNDANT MESSAGE FOR THE 'REMAINING'
!         CONTINUATION CARDS OF ONE 'PARENT'
 
!     NOTES FOR XREAD AND FFREAD ROUTINES, WHICH HAVE DONE SOME
!     IMPORTANT PRELIMINARY TASK -
 
!      1. XSORT2 CALLS XREAD WHICH CALLS FFREAD TO READ ALL INPUT DATA,
!         IN BOTH FIXED-FIELD AND FREE-FIELD FORMATS. UNSORTED INPUT
!         DATA IS NOW PRINTED BY FFREAD IF 'ECHO=UNSORT' IS REQUESTED.
!      2. ALL 10 BULKDATA FIELDS ARE LEFT-ADJUSTED IF INPUT ARE IN
!         FREE-FIELD FORMAT. XREAD LEFT-ADJUSTED ALL FIELDS FOR THE
!         FIXED-FIELD INPUT CASE.
!      3. XREAD PRE-CHECK ANY CONTINUATION, COMMENT, DELETE, BLANK, AND
!         ENDDATA CARDS, AND SET APPROPRIATE FLAGS IN BUF4 CONTROL ARRAY
!      4. THE FIRST THREE BULKDATA FIELDS ARE CONVERTED TO INTERNAL
!         INTEGER CODES AND SAVED IN BUF4 CONTROL ARRAY. THESE INTERNAL
!         CODES ARE READY FOR SORTING.
!      5. XREAD HANDLES BOTH SINGLE-FIELD AND/OR DOUBLE-FIELD INPUT
!         AND PASS ON THE FIRST 3 BULKDATA FIELD INFORMATION INDENTI-
!         CALLY TO THE BUF4 CONTROL ARRAY.
!      6. XREAD/FFREAD COMPLETELY ELIMINATE THE REVERSE-STORAGE PROBLEM
!         OF THE VAX MACHINE.  I.E.
!         THE CONSTANT 'ABCD' IS STORED INTERNALLY AS 'DCBA' IN THE VAX
!      7. IN DOUBLE-FIELD INPUT, THE ASTERISK (*) IN FIELD 1 REMAINS
!         WHERE IT IS. (THE OLD XSORT MOVED IT TO COL. 8 THEN TO COL. 1.
!         SUBROUTINE RCARD MUST BE MODIFIED TO HANDLE THIS DOUBLE-FIELD
!         CASE)
!      8. NO LEADING BCD-ZEROS IN FIELD 2 IF THAT FIELD CONTAINS AN
!         INTEGER NUMBER, AND THE NUMBER IS NOT RIGHT ADJUSTED (I.E.
!         XSORT2 TREATS FIELD 2 INTEGER THE SAME WAY AS INTEGERS IN ALL
!         OTHER FILEDS, NAMELY LEFT ADJUSTED WITH TRAILING BLANKS
!      9. IF THE 1ST FIELD OF THE 2ND CARD IS BLANK, A UNIQUE CONTINUA-
!         TION SYMBOL IS INSERTED INTO THE 1ST FIELD, AND THE SAME
!         SYMBOL IS ADDED TO THE 10TH FIELD OF THE PREVIOUS CARD
 
!     SCRATCH FILE LIMITATION IN LINK1 -
!     SEMDBD ALLOCATES ONLY 15 SCRATCH FILES. SINCE XCSA AND XGPI USE
!     THE LAST SCRATCH FILE FOR RIGID FORMAT, XSORT2, PROGRAMMED UP TO
!     30 FILES, IS THEREFORE PHYSICALLY LIMITTED TO 14 SCRATCH FILES.
 
!     WRITTEN BY G.CHAN/UNISYS   10/1987
 
 IMPLICIT INTEGER (a-z)
 EXTERNAL        lshift,rshift,andf,orf
 LOGICAL :: only1,debug
 INTEGER :: y(25,1),buf(50),ibufx(10),itape(10),temp(2),  &
     NAME(2),bulkda(2),param(2),cdcnt(3),ksmb(3), fub(25)
 CHARACTER (LEN=1) :: ufm*23,uwm*25,uim*29,sfm*25,head4*28,head(3)*56
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /machin/ mach,ijhalf(2),lqro
 COMMON /xsortx/ buf4(4),table(255)
 COMMON /system/ bufsz,nout,nogo,in,dum3(10),  &
     date(4),echo,dum4,apprc,dum5(9),hicore, dum6(7),  &
     nbpc,nbpw,dum7(28),subs,dum8(12),cpflg,dum9(8), lpch
 COMMON /output/ dum10(96),head1(32),head2(32),head3(32)
!ZZ   COMMON /ZZXST2/ Z(1)
 COMMON /zzzzzz/ z(196605)
 COMMON /names / rd,rdrew,wrt,wrtrew,rew,norew,eofnrw
 COMMON /stapid/ dum11(12),kumf
 COMMON /xechox/ ffflag,echou,echos,echop,ixsort,wasff,ncard, f3long,dum12
 COMMON /ifpx0 / dum13(2),ibits(1)
 COMMON /ifpx1 / numx1,icards(2)
 COMMON /two   / itwo(32)
 EQUIVALENCE     (y(1,1),z(1)),      (buf41,buf4(1)),  &
     (ibufx(1),buf(26)), (itape(1),buf(38))
 DATA    head  , head4 /  &
     ' I N P U T   B U L K   D A T A   D E C K   E C H O      ',  &
     '     S O R T E D   B U L K    D A T A    E C H O        ',  &
     ' ---1--- +++2+++ ---3--- +++4+++ ---5--- +++6+++ ---7---',  &
     ' +++8+++ ---9--- +++10+++   '/    ,i25    /25            /
 DATA    NAME         ,cdcnt               ,optp   ,nptp   ,BLANK /  &
     4HXSOR,4HT2  ,4HCARD,4HCOUN,4HT   ,4HOPTP ,4HNPTP ,4H    /
 DATA    tape1 ,tape2 ,tape3 ,maxscr,bulkda        ,param         /  &
     301   ,302   ,303   ,314   ,4HBULK,4HDATA ,4HPARA ,4HM   /
 DATA    ksmb  /4H+c0n,4H+cqn,4H+con  /    ,debug  /.false.       /
 
!     DIAG 47 CAN BE RE-ACTIVATED FOR PROGRAM DEBUG CHECKING
 
!     CALL SSWTCH (47,J)
!     IF (J .EQ. 1) DEBUG = .TRUE.
 
!     TURN ON XSORT FLAG AND FREE-FIELD FLAG FOR XREAD AND FFREAD
 
 ixsort = 1
 ffflag = 1234
 
!     CHECK UMF REQUEST
 
 IF (kumf <=  0) GO TO 110
 WRITE  (nout,100) ufm
 100 FORMAT (a23,' - USER MASTER FILE, UMF, IS  NOT SUPPORTED BY NEW ',  &
     'XSORT ROUTINE', /5X, 'ADD A ''DIAG 42'' CARD AND RESUBMIT YOUR NASTRAN JOB')
! 100 FORMAT (A23,' - USER MASTER FILE, UMF, IS NO LONGER SUPPORTED BY',
!    1        ' NASTRAN',/5X,'(NOTE - RELEASE 87 WAS THE LAST VERSION ',
!    2        'THAT SUPPORTED UMF OPERATION)')
 CALL mesage (-37,0,NAME)
 
!     INITIALIZE XSORT2
 
 110 echou = 0
 echos = 0
 echop = 0
 ncard = 0
 cmmt  = 0
 ncont = 0
 ndele = 0
 full  = 0
 exh   = 0
 tapecc= 0
 bsize = 3
 restr = 0
 case  = 1
 kontn = 10010000
 ksmbi = ksmb(1)
 IF (apprc < 0) restr = 1
 IF (restr == 1) ksmbi = khrfn3(ksmb(1),date(2),-2,0)
 j     = complf(0)
 large = rshift(j,1   )
 les1b = rshift(j,nbpc)
 IF (MOD(lqro,10) == 1) les1b = lshift(j,nbpc)
 IF (echo < 0) GO TO 120
 echou = andf(echo,1)
 echos = andf(echo,2)
 echop = andf(echo,4)
 IF (cpflg /= 0) echos = 1
 
!     SET UP UNSORTED HEADING
 
!     (UNSORTED INPUT DATA IS NOW PRINTED BY FFREAD ROUTINE BECAUSE
!      XREAD HAS BEEN MODIFIED TO RETURN ALL 10 DATA FIELDS LEFT-
!      ADJUSTED)
 
 120 DO  j = 1,32
   head2(j) = BLANK
   head3(j) = BLANK
   head1(j) = BLANK
 END DO
 imhere = 130
 IF (debug) WRITE (nout,140) imhere,restr,apprc,subs
 140 FORMAT (//,' *** XSORT2/IMHERE =',6I5)
 READ (head(1),150) (head1(j),j=11,24)
!WKBR 9/93 READ (HEAD(3),150) (HEAD3(J),J= 7,20)
 READ (head(3),150) (head3(j),j= 8,21)
!WKBR 9/93 READ (HEAD 4 ,150) (HEAD3(J),J=21,27)
 READ (head 4 ,150) (head3(j),j=22,28)
 150 FORMAT (14A4)
 IF (echou /= 0) CALL page
 
!     GET AVAILABLE CORE
!     IF IBM MACHINE, LIMIT AVAILABLE CORE SIZE TO 1,000,000 WORDS, SUCH
!     THAT DATA WILL BE SAVED IN PRIMARY FILES ONLY, AND NO SPILL INTO
!     SECONDARY FILES.
 
 nzz   = korsz(z)
 ibuf1 = nzz   - bufsz
 ibuf2 = ibuf1 - bufsz
 ibuf3 = ibuf2 - bufsz
 nz    = ibuf3 - 1
 IF (mach  == 2) nz = MIN0(nz,1000000)
 IF (nz < 2500) CALL mesage (-8,2500,NAME)
 nz25  = nz/25
 
!     OPEN TAPE1, GINO FILE 301 FOR DELETE (SLASH) CARDS
!     AND  TAPE2, GINO FILE 302 FOR CONTINUATION CARDS
!     SET  TAPE TO TAPE3, GINO FILE 303, FOR BULKDATA CARDS
!     UP TO 30 FILES ARE ALLOWED FOR REGUALR BULKDATA CARDS
!     (CURRENTLY /XFIST/ IN SEMDBD IS SET UP ONLY TO SCRATCH FILE 315.
!     I.E. UP TO 13 (OR 12, IF DECK CONTAINS MANY CONTINUATION CARDS)
!     FILES CAN BE USED HERE)
 
 imhere = 170
 IF (debug) WRITE (nout,140) imhere,nz25
 CALL OPEN (*2900,tape1,z(ibuf1),wrtrew)
 CALL OPEN (*2910,tape2,z(ibuf2),wrtrew)
 tape = tape3 - 1
 170 tape = tape  + 1
 IF (tape <= 314) GO TO 180
 IF (debug) WRITE (nout,2955)
 CALL mesage (-8,-nzz,NAME)
 180 CALL OPEN (*2960,tape,z(ibuf3),wrtrew)
 wrttn = 0
 
 
!     START READING INPUT CARDS VIA XREAD/FFREAD.
 
 
!     ADDITIONAL INFORMATION FROM XREAD NOT MENTIONED PREVIOUSLY -
 
!      1. BUF4(1) = BUF4(2) =-1 INDICATE BULKDATA IS A COMMENT CARD
!         BUF4(1) = BUF4(2) =-2 INDICATE BULKDATA IS A CONTINUATION CARD
!         BUF4(1) = BUF4(2) =-3 INDICATE BULKDATA IS A DELETE CARD, WITH
!                   DELETE RANGE SAVED IN BUF4(3) AND BUF4(4)
!         BUF4(1) =-3 AND BUF4(4) =-4 IF TRASH WAS FOUND IN DELETE CARD.
!                   THAT IS, TRASH AFTER SLASH IN BULKDATA FIELD 1
!         BUF4(1) = BUF4(4) =-5 INDICATE A  BLANK   CARD WAS READ
!         BUF4(1) = BUF4(4) =-9 INDICATE AN ENDDATA CARD WAS READ
!      2. IF BULKDATA FIELD 2 IS AN INTEGER INPUT, THE CORRECT INTEGER
!                 VALUE IS SAVED IN BUF4(3)
!         IF BULKDATA FIELD 3 IS AN INTEGER INPUT, THE CORRECT INTEGER
!                 VALUE IS SAVED IN BUF4(4)
!      3. IF THE DATA IN FIELD 2 AND/OR 3 ARE F.P. NUMBER, THEIR INTEGER
!                 VALUES (NOT EXACT) ARE SAVED IN BUF4(3) AND/OR BUF4(4)
!                 THESE VALUES ARE USED ONLY FOR SORTING
!      4. IF BULKDATA FIELD 2 IS NOT NUMERIC, THE FIRST 6 CHARACTERS ARE
!                 CONVERTED TO INTERNAL INTEGER CODE AND SAVED IN BUF4(3
!         IF THE REMAINING 2 CHARACTERS ARE NOT BLANKS, THEY ARE SAVED
!                 IN BUF4(4)
!      5. IF BUF4(4) IS NOT USED BY 4, IT HOLDS THE INTERNAL CODE OR THE
!                 INTEGER VALUE FOR FIELD 3 OF THE ORIGINAL BULKDATA.
 
!     WORK SPACE -                                     NZ
!      1                                               /
!     ------------------------------------------------------------------
!     !                 OPEN CORE, Z                    !    !    !    !
!     ------------------------------------------------------------------
!     !<----------INPUT CARDS, 25 WORDS EACH----------->!<----GINO---->!
!                (20-WORD CARD IMAGE, 4 CONTRL               BUFFERS
!               CONTROL WORDS, 1 INDEX POINTER)
 
 
!     SUMMARY OF COUNTERS -
 
!     NCONT = TOTAL CONTINUATION CARDS COUNT, ON INPUT BULK DATA DECK
!             AND ON RESTART OPTP FILE
!     NDELE = TOTAL COUNT ON RESTART DELETE CARDS
!     CMMT  = TOTAL COUNT ON NON-ESSENTIAL CARDS (COMMENTS, BLANKS, AND
!             RESTART DELETE CARDS) OF INPUT BULK DATA DECK
!     KONTN = SYMBOL COUNTER FOR AUTO-CONTINUAION GENERATION
!     KOUNT = DELETE RANGE COUNTER, USED ONLY IN 800-820 AREA
!     NCARD = TOTAL INPUT BULK DATA CARDS COUNT, INCLUDING NON-ESSENTIAL
!             CARDS; CONTINUATION CARDS AND CARDS ON OPTP ARE EXCLUDED
!     COUNT = CURRENT CORE COUNT ON INPUT CARDS FROM BULK DATA DECK, ALL
!             NON-ESSENTIAL AND CONTINUATION CARDS ARE EXCLUDED
!     NBULK = NO. OF ACTIVE BULK DATA INPUT CARDS
!           = NCARD-CMMT = SUM OF ALL COUNT's
!     NOTE  - NO CARD COUNT ON THE OPTP FILE BEFORE ASSEMBLING NPTP FILE
 
 count = 0
 200 IF (count < nz25) SELECT CASE ( case )
   CASE (    1)
     GO TO 212
   CASE (    2)
     GO TO 214
   CASE (    3)
     GO TO 207
   CASE (    4)
     GO TO 210
   CASE (    5)
     GO TO 210
   CASE (    6)
     GO TO 210
 END SELECT
!                                   1,  2,  3,  4,  5,  6 = CASE
 case = 1
 IF (wasff <= 0) GO TO 340
 
!     (200 THRU 215) SPECIAL HANDLING OF CONTINUATION CARD(S) WITH FIRST
!     FIELD BLANK DURING FREE-FIELD INPUT.   REGULAR CONTINUATION CARD
!     (FIRST FIELD NOT BLANK) OR FIXED-FIELD INPUT CARDS (BOTH PARENT
!     AND CHILD) ARE NOT CONSIDERED HERE.
 
!        EXAMPLE -     CBAR,10,20, 1 2 3  9)2
!                      ,,, .5 .5 .5
 
!     WE NEED TO CREATE A UNIQUE CONTINUATION SYMBOL FOR THE 1ST FIELD,
!     AND ADD THE SAME SYMBOL TO THE 10TH FIELD OF THE PREVIOUS CARD.
!     SET BUF41 FLAG TO -2.
!                                                                WAITING
!     AT THIS POINT,                                             CARD IN
!        CASE 1, NO CARD IS WAITING FOR PROCESSING               -------
!        CASE 2, CORE WAS FULL AND WAS EMPTIED OUT. A NON-           BUF
!                CONTINUATION CARD WAS READ AND AWAITS PROCESSING
!        CASE 3, CORE WAS FULL AND EMPTIED. A CONTINUATION CARD      BUF
!                WAS READ AND AWAITS PROCESSING.
!        CASE 4, CORE NOT FULL, A CONT.CARD WAS READ. THE NEXT CARD  FUB
!                IS NOT A CONT.CARD. THE CONT.CARD WAS PROCESSED,
!                AND THE NON-CONT. CARD  AWAITS PROCESSING.
!        CASE 5, CORE NOT FULL, A CONT.CARD WAS READ AND THE NEXT    FUB
!                CARD IS ALSO A CONT.CARD. THE FIRST CONT.CARD
!                WAS PROCESSED, AND THE SECOND CONT.CARD AWAITS
!                PROCESSING.
!        CASE 6, CONTINUE FROM PROCESSING CASES=4,5                  FUB
 
! ... CASES 2 AND 3 -
!     CORE IS FULL, READ ONE MORE CARD AND SEE THE NEW CARD IS A SPECIAL
!     CONTINUATION CARD OR NOT
!     IF IT IS, UPDATE THE 10TH FIELD OF THE PARENT CARD BEFORE
!     SENDING THE ENTIRE CORE FOR SORTING
 
 imhere = 202
 202 CALL xread (*208,buf)
 IF (buf41 == -1 .OR. buf41 == -5) GO TO 202
 case = 2
 IF (buf(1) /= BLANK .OR. buf(2) /= BLANK) IF (buf41+2) 340,203,340
 203 buf41x = -2
 case = 3
 GO TO 205
 
! ... CASES 4 AND 5 -
!     CORE IS NOT FULL, A SPECIAL CONTINUATION CARD WAS JUST READ
 
 204 IF (wasff <= 0) GO TO 214
 case  = 4
 buf41 = -2
 205 kontn = kontn + 1
 IF (kontn == 10020000) ksmbi = ksmb(2)
 IF (kontn == 10030000) ksmbi = ksmb(3)
 imhere = 205
 IF (debug) WRITE (nout,140) imhere,kontn,count,nz25,case
 CALL int2a8 (*3140,kontn,buf(1))
 buf(1) = ksmbi
 IF (count <= 0) GO TO 208
 y(19,count) = buf(1)
 y(20,count) = buf(2)
 IF (case-3 > 0.0) THEN
   GO TO   207
 ELSE
   GO TO   340
 END IF
 
 206 case = 6
 IF (buf41 == -9) GO TO 350
 
 207 CALL xread (*207,fub)
 IF (buf41 == -1 .OR. buf41 == -5) GO TO 207
 fub41 = buf41
 buf41 = -2
 IF (fub(1) /= BLANK .OR. fub(2) /= BLANK) GO TO 215
 fub41 = -2
 case  = 5
 kontn = kontn + 1
 IF (kontn == 10020000) ksmbi = ksmb(2)
 IF (kontn == 10030000) ksmbi = ksmb(3)
 imhere  = 207
 IF (debug) WRITE (nout,140) imhere,kontn,count,nz25,case
 CALL int2a8 (*3140,kontn,fub(1))
 fub(1)  = ksmbi
 buf(19) = ksmbi
 buf(20) = fub(2)
 GO TO 217
 
 208 nogo  = 1
 WRITE  (nout,209) sfm,imhere
 209 FORMAT (a25,'.  IMHERE =',i6)
 GO TO 214
 
 210 DO  i = 1,25
   buf(i) = fub(i)
 END DO
 buf41  = fub41
 IF (case-5 == 0.0) THEN
   GO TO   206
 ELSE
   GO TO   214
 END IF
 212 CALL xread (*3120,buf)
 IF (buf(1) == BLANK .AND. buf(2) == BLANK) GO TO 204
 214 case = 1
 
!     IGNORE COMMENT CARD (-1) OR BLANK CARD (-5)
 
 215 IF (buf41 /= -1 .AND. buf41 /= -5) GO TO 216
 cmmt = cmmt + 1
 GO TO 212
 
!     TEST FOR ENDDATA CARD (-9)
 
 216 IF (buf41 == -9) GO TO 350
 
!     IF THIS IS A CONTINUATION CARD (-2), ADD ONE CONTROL WORD ABOUT
!     RESTART, AND WRITE IHE CARD OUT TO TAPE2
!     (THE CONTROL WORD WILL FLAG THE PARENT BIT TO BE SET FOR RESTART
!     WHEN THIS CONTINUATION CARD IS MERGED INTO NPTP)
 
 IF (buf41 /= -2) GO TO 230
 217 buf(21) = restr
 CALL WRITE (tape2,buf(1),21,0)
 IF (debug) WRITE (nout,220) buf(1),buf(2),buf(21)
 220 FORMAT (5X,'A CONTINUATION CARD - ',2A4,',  CONT.FLAG=',i9)
 ncont = ncont + 1
 GO TO 200
 
!     IF THIS IS A DELETE CARD (-3), REJECT IT IF EXTRANEOUS DATA IN
!     FIELD 1 OTHERWISE WRITE THE RANGE OF DELETION ON TAPE1
 
 230 IF (buf41 /= -3) GO TO 300
 cmmt = cmmt + 1
 IF (buf4(4) /= -4) GO TO 250
 CALL page2 (2)
 WRITE  (nout,240) ufm
 240 FORMAT (a23,' 221, EXTRANEOUS DATA IN FIELD 1 OF BULK DATA ',  &
     'DELETE CARD.')
 nogo = -2
 
 250 IF (buf4(3) == -3) GO TO 270
 IF (buf4(4) == -3) buf4(4) = buf4(3)
 buf4(3) = buf4(3) - 2000000000
 buf4(4) = buf4(4) - 2000000000
 CALL WRITE (tape1,buf4(3),2,0)
 IF (debug) WRITE (nout,260) buf4(3),buf4(4)
 260 FORMAT (5X,'A DELETE CARD -',i11,1H,,i11)
 ndele  = ndele + 1
 GO TO 200
 270 WRITE  (nout,280) ufm
 280 FORMAT (a23,' 221, NO DATA IN FIELD 2 OF BULK DATA DELETE CARD')
 nogo = -1
 GO TO 200
 
!     REGULAR BULKDATA CARDS.
!     SAVE 20 WORDS OF BUF, 4 WORDS FROM BUF4 AND CORE COUNTER IN OPEN
!     CORE SPACE Y (25 WORDS TOTAL)
!     SET RESTART BITS IF THIS IS A RESTART RUN
!     RETURN TO READ NEXT BULKDATA CARD
 
 300 count = count + 1
 wrttn = 1
 DO  i = 1,20
   y(i,count) = buf(i)
 END DO
 DO  i = 1,4
   y(i+20,count) = buf4(i)
 END DO
 y(25  ,count) = count
 IF (debug) WRITE (nout,330) count,y(1,count),y(2,count)
 330 FORMAT (5X,'SAVED IN CORE   COUNT=',i5,3X,2A4)
 IF (restr == 0) GO TO 200
 ASSIGN 200 TO crdflg
 from = 330
 GO TO 2800
 
!     OPEN CORE BUFFER FULL, ENDDATA CARD HAS NOT BEEN ENCOUNTERED
 
 340 full = 1
 GO TO 400
 
!     ENDDATA CARD FOUND, SET FLAG
 
 350 full  = -1
 imhere= 350
 ncard = ncard - 1
 nbulk = ncard - cmmt
 IF (debug) WRITE (nout,140) imhere,ncard,ncont,ndele
 CALL page2 (2)
 IF (echou /= 1) GO TO 370
 WRITE  (nout,360) ncard
 360 FORMAT (//24X,'TOTAL COUNT=',i7)
 GO TO 400
 370 WRITE  (nout,380) ncard,cmmt
 380 FORMAT (//24X,'(NO. OF UNSORTED BULK DATA CARDS READ =',i6,  &
     ', INCLUDING',i4,' COMMENT CARDS)')
 
!     SORT CARD IMAGES SAVED IN THE OPEN CORE SPACE BY MODIFIED SHELL
!     METHOD.
!     SORT BY 21ST, 22ND, 23RD, AND 24TH CONTROL WORDS ONLY
!     ONLY THE LAST 5 WORDS (21ST THRU 25TH) ARE MOVED INTO SORTED
!     ORDER, THE FIRST 20 WORDS REMAIN STATIONARY.
 
 400 IF (wrttn ==    0) GO TO 580
 IF (count > nz25) CALL mesage (-37,0,NAME)
 m  = count
 imhere = 400
 IF (debug) WRITE (nout,140) imhere,count
 410 m  = m/2
 IF (m == 0) GO TO 500
 j  = 1
 k  = count - m
 420 i  = j
 430 n  = i + m
 IF (y(21,i) - y(21,n) < 0.0) THEN
   GO TO   490
 ELSE IF (y(21,i) - y(21,n) == 0.0) THEN
   GO TO   440
 ELSE
   GO TO   470
 END IF
 440 IF (y(22,i) - y(22,n) < 0.0) THEN
   GO TO   490
 ELSE IF (y(22,i) - y(22,n) == 0.0) THEN
   GO TO   450
 ELSE
   GO TO   470
 END IF
 450 IF (y(23,i) - y(23,n) < 0.0) THEN
   GO TO   490
 ELSE IF (y(23,i) - y(23,n) == 0.0) THEN
   GO TO   460
 ELSE
   GO TO   470
 END IF
 460 IF (y(24,i) - y(24,n) > 0.0) THEN
   GO TO   470
 ELSE
   GO TO   490
 END IF
 470 DO  l = 21,25
   tempx  = y(l,i)
   y(l,i) = y(l,n)
   y(l,n) = tempx
 END DO
 i = i - m
 IF (i >= 1) GO TO 430
 490 j = j + 1
 IF (j-k > 0) THEN
   GO TO   410
 ELSE
   GO TO   420
 END IF
 
!     END OF CORE SORT.
!     WRITE THE SORTED BULKDATA CARDS TO FILE, 24 WORDS EACH RECORD
!     IN ORDER GIVEN BY THE 25TH WORD.
!     IF ONLY ONE SCRATCH FILE (TAPE3) IS USED IN RECEIVING BULKDATA,
!     CHECK ANY DUPLICATE CARD.
 
 500 imhere = 500
 IF (debug) WRITE (nout,140) imhere,count,maxc
 only1 = .false.
 IF (full == -1 .AND. tape == tape3) only1=.true.
 base = 25
 DO  i = 1,count
   IF (only1) base = MOD(i,2)*25
   j = y(25,i)
   DO  k = 1,20
     buf(k+base) = y(k,j)
   END DO
   DO  k = 21,24
     buf(k+base) = y(k,i)
   END DO
   IF (.NOT.only1) GO TO 550
   IF (i  ==  1) GO TO 540
   DO  k = 1,20
     IF (buf(k+base) /= buf(k+obase)) GO TO 540
   END DO
   buf(21+base) = -6
   buf(22+base) = -6
   540 obase = base
   550 CALL WRITE (tape,buf(base+1),24,0)
   IF (debug) WRITE (nout,560) tape,(buf(k+base),k=1,8)  &
       ,(buf(k+base),k=21,24)
   560 FORMAT (5X,'WRITE TO ',i3,4(2X,2A4), /9X,'INT.CODE=',4I12)
 END DO
 CALL WRITE (tape,0,0,1)
 580 CALL CLOSE (tape,rew)
 imhere = 580
 IF (debug) WRITE (nout,140) imhere
 
!     REPEAT READING BULKDATA CARDS INTO CORE IF NECESSARY
 
!     IF NO DATA WRITTEN TO CURRENT FILE (e.g. UN-MODIFIED RESTART),
!     REDUCE TAPE COUNT BY ONE
 
 IF (full /= -1) GO TO 170
 IF (wrttn == 0) tape = tape - 1
 
!     CLOSE DELETE CARD FILE, TAPE 1.
!     CONTINUATION CARD FILE, TAPE 2, IS STILL IN USE
 
 CALL WRITE (tape1,0,0,1)
 CALL CLOSE (tape1,rew  )
 
!     TEST FOR COLD-START WITH NO BULKDATA
 
!     APPRC = APPROACH FLAG (1 DMAP, 2 DISP, 3 HEAT, 4 AERO)
!     SUBS  = SUBSTRUCTURING FLAG
 
 imhere = 585
 IF (debug) WRITE (nout,140) imhere,count,apprc,wrttn,restr,subs
 IF (wrttn == 1 .OR. restr == 1 .OR. subs /= 0) GO TO 600
 CALL CLOSE (tape2,rew)
 echos = 1
 IF (apprc == 1) GO TO 1600
 CALL page2 (2)
 WRITE  (nout,590) ufm
 590 FORMAT (a23,' 204, COLD START NO BULK DATA.')
 nogo = -2
 GO TO 3200
 
!     IF MODIFIED RESTART - TURN ON SORT ECHO FLAG IF ECHO IS NOT 'NONO'
!     IF NOT A RESTART JOB - JUMP TO 1000
 
 600 IF (nbulk > 1 .AND. restr == 1) echos = 1
!     IF (APPRC.EQ.1 .OR.  SUBS .NE.0) ECHOS = 1
 IF (echo  == -2) echos = 0
 IF (restr ==  0) GO TO 1000
 
!     THIS IS A RESTART JOB, PROCESS OPTP FILE -
 
!     OPEN OPTP AND LOCATE WHERE BULK DATA BEGINS
 
 imhere = 610
 IF (debug) WRITE (nout,140) imhere
 CALL OPEN (*3080,optp,z(ibuf3),rdrew)
 610 CALL skpfil (optp,+1)
 CALL READ (*3040,*3040,optp,buf(1),2,1,j)
 IF (buf(1) /= bulkda(1) .OR. buf(2) /= bulkda(2)) GO TO 610
 IF (nbulk > 0 .OR. ndele /= 0) GO TO 640
 
!     UN-MODIFIED RESTART, WITH NO NEW BULKDATA CARD AND NO DELETE -
!     SETUP SORTED HEADER FOR OLD BULK DATA CARDS IF ECHO FLAG IS ON,
!     COPY THE REST OF OPTP DIRECTLY TO NPTP, AND JOB DONE
 
 
 imhere = 620
 IF (debug) WRITE (nout,140) imhere
 CALL OPEN  (*3100,nptp,z(ibuf1),wrt)
 CALL WRITE (nptp,bulkda,2,1)
 ncard = 0
 IF (echos == 0) GO TO 620
 READ (head(2),150) (head1(j),j=11,24)
!WKBR 9/93 HEAD2(4) = CDCNT(1)
 head2(5) = cdcnt(1)
!WKBR 9/93 HEAD3(4) = CDCNT(2)
 head3(5) = cdcnt(2)
!WKBR 9/93 HEAD3(5) = CDCNT(3)
 head3(6) = cdcnt(3)
 CALL page
 620 CALL READ  (*630,*630,optp,buf(1),20,1,j)
 CALL WRITE (nptp,buf(1),20,1)
 ncard = ncard + 1
 IF (echop /= 0) WRITE (lpch,1750) (buf(j),j=1,20)
 IF (echos == 0) GO TO 620
 CALL page2 (-1)
 WRITE (nout,1730) ncard,(buf(j),j=1,20)
 GO TO 620
 630 CALL eof   (nptp)
 CALL CLOSE (nptp,  rew)
 CALL CLOSE (optp,norew)
 CALL CLOSE (tape2, rew)
 IF (echop /= 0) WRITE (lpch,2320)
 CALL page2 (-1)
 IF (echos /= 0) WRITE (nout,2300)
 IF (echos == 0) WRITE (nout, 635) uim,ncard
 635 FORMAT (a29,1H,,i8,' SORTED BULKD DATA CARDS PROCESSED FROM OPTP',  &
     ' FILE TO NPTP, UN-MODIFIED')
 GO TO 2700
 
!     MODIFIED RESTART WITH NEW BULKDATA CARDS, WITH OR WITHOUT DELETE
 
 640 imhere = 640
 IF (debug) WRITE (nout,140) imhere
 ic   = 1
 left = nz
 IF (ndele == 0) GO TO 710
 IF (restr == 1) GO TO 660
 CALL page2 (-1)
 WRITE  (nout,650) uwm
 650 FORMAT (a25,' 205, COLD START, DELETE CARDS IGNORED.')
 GO TO 710
 
!     RESTART WITH DELETE CARD(S) -
!     MOVE THE DELETE CARDS  INTO CORE AND FREE TAPE1.
!     SORT THE DELETE CARDS, CHECK FOR AND ELIMINATE OVERLAPS AND
!     REDUNDANCIES
 
 660 CALL OPEN (*2900,tape1,z(ibuf1),rdrew)
 CALL READ (*2900,*670,tape1,z(1),left,1,LEN)
 CALL mesage (-8,tape1,NAME)
 670 CALL CLOSE  (tape1,rew    )
 
 CALL sort (0,0,2,1,z(1),LEN)
 z(LEN+1) = large
 DO  i = 2,LEN,2
   z(i) = z(i)+1
   IF (z(i) < z(i-1)) z(i) = z(i-1)
   IF (z(i) < z(i+1)) CYCLE
   z(i  ) = -1
   z(i+1) = -1
 END DO
 j = 0
 DO  i = 1,LEN
   IF (z(i) < 0) CYCLE
   j = j + 1
   z(j) = z(i)
 END DO
 IF (j > 0) LEN = j
 left = nz - LEN
 ic   = LEN + 1
 z(ic)  = large
 imhere = 700
 IF (debug) WRITE (nout,700) imhere,(z(i),i=1,LEN)
 700 FORMAT (/,' *** IMHERE =',i5,(/,3X,10(i7,i5)))
 IF (MOD(LEN,2) /= 0) GO TO 3140
 GO TO 800
 
!     IF MODIFIED RESTART WITH NO DELETE, SET DELETE RANGE BEGINNING AT
!     INFINITY
 
 710 z(1)   = large
 imhere = 710
 IF (debug) WRITE (nout,140) imhere
 
!     WE ARE STILL IN PROCESSING RESTART - COPY OPTP TO TAPE1, SKIP
!     APPROPRIATE RECORDS AS SPECIFIED BY THE DELETE CARDS NOW IN
!     OPEN CORE, Z(1) THRU Z(LEN)
 
!     SEND A CARD FROM OPTP TO YREAD (AN ENTRY POINT IN XREAD) FOR
!     RE-PROCESSING. UPON RETURN FROM YREAD, BUF4 ARRAY HOLDS THE
!     INTERNAL INTEGER CODE GOOD FOR SORTING AND OTHER FUNCTIONS.
 
!     IF IT IS A CONTINUATION CARD, COPY THE FULL CARD (20 WORDS)
!     AND ONE CONTROL WORD TO TAPE2.
!     OTHERWISE COPY 24 WORDS (20-BUF AND 4-BUF4) TO TAPE1.
 
!     IF A CONTINUATION CARD IS DELETED, THE RESTART BITS OF THE
!     PARENT CARD SHOULD BE FLAGGED
 
 800 imhere  = 800
 IF (debug) WRITE (nout,140) imhere,restr,tape1
 CALL OPEN (*2900,tape1,z(ibuf1),wrtrew)
 kount   = 0
 point   = 1
 onoff   = 1
 zpoint  = z(point)
 buf(19) = 0
 810 temp(1) = buf(19)
 temp(2) = buf(20)
 CALL READ (*900,*900,optp,buf(1),20,1,j)
 kount   = kount + 1
 IF (kount < zpoint) GO TO 820
 point   = point + 1
 zpoint  = z(point)
 onoff   = onoff*(-1)
 820 CALL yread (*3060,buf)
 imhere  = 830
 IF (debug .AND. onoff == -1) WRITE (nout,830) imhere,kount, (buf(j),j=1,6)
 830 FORMAT (' IMHERE=',i5,'.  DELETED FROM OPTP ==>',i5,2H- ,6A4)
 IF (buf41 == -2) GO TO 870
 IF (onoff == +1) GO TO 840
 
!     ANY DELETED CARD, EXCEPT CONTINUATION CARD, MUST RESET
!     RESTART CARD FLAG
 
 ASSIGN 810 TO crdflg
 from = 830
 GO TO 2800
 
!     REGULAR BULKDATA CARD FROM OPTP -
!     SAVE FIRST FIELD IN KARD1/2 JUST IN CASE THIS IS A PARENT OF
!     A CONTINUATION CARD WHICH FALLS INSIDE A DELETE RANGE.
 
!     NOTE- CARDS FROM OPTP ARE IN SORTED ORDER, AND NO CARD COUNT HERE
 
 840 DO  j = 1,4
   buf(j+20) = buf4(j)
 END DO
 CALL WRITE (tape1,buf(1),24,0)
 IF (debug) WRITE (nout,860) (buf(j),j=1,6),buf(21)
 860 FORMAT (' IMHERE=860, OPTP==>TAPE1  ',6A4,'==>',i9)
 kard1 = buf(1)
 kard2 = buf(2)
 IF (kard1 /= param(1) .OR. kard2 /= param(2)) GO TO 810
 kard1 = buf(3)
 kard2 = buf(4)
 GO TO 810
 
!     CONTINUATION CARD FROM OPTP -
 
!     IF BOTH PARENT AND THIS CONTINUATION CARD IN NOT IN DELETE RANGE
!     SEND THIS CONTINUATION CARD TO TAPE2 WITH RESTART CONTROL WORD
!     SET TO ZERO.
!     IF PARENT IS NOT DELETED, BUT THIS CONTINUATION CARD IS, WE NEED
!     TO FLAG PARENT
!     IF PARENT IS ALSO IN DELETE RANGE, SKIP THIS CONTINUATION CARD.
 
 870 IF (onoff == +1) GO TO 890
 IF (kard1 == -1) GO TO 810
 IF (buf(1) == temp(1) .AND. buf(2) == temp(2)) GO TO 810
 from  = 860
 ASSIGN 880 TO crdflg
 GO TO 2810
 880 kard1 = -1
 GO TO 810
 890 buf(21) = 0
 CALL WRITE (tape2,buf(1),21,0)
 ncont = ncont + 1
 GO TO 810
 
!     OPTP IS SUCCESSFULLY MOVED TO TAPT1 AND TAPE2. CLOSE FILES
 
 900 CALL CLOSE (optp ,norew)
 CALL WRITE (tape1,0,0,1)
 CALL WRITE (tape2,0,0,1)
 CALL CLOSE (tape1,rew  )
 
!     PREPARE FOR FILE MERGE -
 
!     SELECT METHOD USED TO BRING CONTINUATION CARDS INTO CORE AND
!     COMPUTE NUMBER OF BUFFERS NEEDED FOR FILE PRE-MERGE.
 
!     METHOD 1 - NO FILE PRE-MERGE IF THERE IS NO CONINUATION CARDS, OR
!                ENOUGH SPACE IN CORE TO HOLD ALL CONTINUATION CARDS,
!                BUFFERS AND SCRATCH ARRAYS FOR ALL SCRATCH DATA FILES
!     METHOD 2 - ALL CONTINUATION CARDS, IN 3-WORD TABLE AND 20-WORD
!                CARD IMAGES, AND ALL GINO BUFFERS, OR REDUCED GINO
!                BUFFERS, FIT INTO CORE
!     METHOD 3 - CONTINUATION 3-WORD TABLE AND ALL GINO BUFFERS, OR
!                REDUCED GINO BUFFERS, FIT INTO CORE
!     METHOD 4 - FATAL, INSUFFICIENT CORE
 
 1000 CALL CLOSE (tape2,rew)
 method = 1
 n23    = 1
 nfiles = tape - tape3 + 1
 reduce = 1
 IF (nfiles >= 10) reduce = 2
 IF (nfiles > 17) reduce = 3
 j      = 0
 IF (restr == 1) j = 1
 maxc   = (nzz-(bufsz+25)*(nfiles+j))/21
 IF (ncont <= maxc) reduce = 1
 nfiler = (nfiles+reduce-1)/reduce + j
 imhere = 1010
 IF (debug) WRITE (nout,140) imhere,reduce,nfiles,nfiler
 IF (ncont == 0) THEN
   GO TO  1100
 ELSE
   GO TO  1020
 END IF
 1010 size   = (nfiler+1)*bufsz + nfiler*25
 size   = size + bufsz
 left   = nzz - size
 maxc   = left/n23
 imhere = 1020
 IF (debug) WRITE (nout,140) imhere,method,nfiles,nfiler,n23,ncont
 IF (ncont <= maxc) GO TO 1100
 SELECT CASE ( method )
   CASE (    1)
     GO TO 1020
   CASE (    2)
     GO TO 1030
   CASE (    3)
     GO TO 1040
 END SELECT
 1020 method = 2
 n23    = 23
 GO TO 1010
 1030 method = 3
 n23    = 3
 GO TO 1010
 
!     INSUFFICIENT CORE, COMPUTE HOW MUCH MORE NEEDED
 
 1040 j = ncont*n23 - left
 CALL mesage (-8,j,NAME)
 
!     ALLOCATE BUFFER SPACE AND REDEFINE AVAILABLE CORE SPACE, NZ
!     ALLOCATE SPACES AT THE BEGINNING OF CORE SPACE FOR BULKDATA
!     TO BE BROUGHT BACK FROM VARIOUS FILES.
 
!     IC     = POINTER, WHERE CONTINUATION TABLE BEGINS
!     IB     = POINTER, WHERE CONTINUATION  DATA BEGINS
!     NFILES = TOTAL NUMBER OF FILES USED BEFORE FILE REDUCTION,
!              RESTART TAPE1 NOT INCLUDED
!     NFILER = REDUCED NUMBER OF FILES THAT HOLD BULKDATA INPUT CARDS,
!              RESTART TAPE1 INCLUDED
!     TAPECC = AN ADDITIONAL FILE USED ONLY IN METHOD 3 (NOT INCLUDED
!              IN NFILES AND NFILER)
 
 1100 imhere = 1100
 IF (debug .OR. nfiles > 10 .OR. ncont > 1000)  &
     WRITE (nout,1110) uim,method,nfiler,hicore,ncont
 1110 FORMAT (a29,' FROM XSORT -  METHOD',i3,' WAS SELECTED TO PROCESS',  &
     ' CONTINUATION CARDS', /5X,'NO. OF FILES USED =',i4,4X,  &
     'HICORE =',i7,' WORDS', 4X,'NO. OF CONT. CARDS =',i7)
 nz   = ibuf1
 DO  i = 1,nfiler
   nz   = nz - bufsz
   ibufx(i) = nz
 END DO
 IF (ncont > 0) nz = nz - bufsz
 ibufc= nz
 nz   = nz - 1
 ic   = nfiler*25 + 1
 ib   = ic + ncont*3
 nzib = nz - ib + 1
 left = nz - ic + 1
 
!     NEED A STORAGE SPACE FOR AT LEASE 100 CONTINUATION CARDS
 
 IF (nzib < 2100) CALL mesage (-8,-2100+nzib,NAME)
 
!     METHOD 1, NO CONTINUATION CARD IN BULKDATA, SKIP TO 1280
 
 IF (method == 1) GO TO 1280
 
!     WORKING SPACE FOR THE CONTINUATION TABLE AND CONTINUATION CARD
!     IMAGES -
 
!                  IC                 IB                   NZ
!                  /                  /                    /
!     ------------------------------------------------------------------
!     ! ! ! !..Y..!                  !                     !  !  !  !  !
!     ------------------------------------------------------------------
!     ! SPACE FOR !<--CONTINUATION-->!<--AVAILABLE SPACE-->!<--GINO--->!
!       DATA FROM     INDEX TABLE        FOR CONTINUATION     BUFFERS
!       FILES 303,   (3 WORDS EACH)      CARD IMAGES
!       304,...                          (21 WORDS EACH)
!       FOR FILE      (PART 1 AERA)
!       MERGE                            (PART 2 AREA)
 
 imhere = 1125
 IF (debug) WRITE (nout,140) imhere,method,n23
 CALL OPEN (*2910,tape2,z(ibuf2),rdrew)
 IF (method == 3) GO TO 1200
 
!     METHOD 2 -
 
!     OPEN CORE IS DIVIDED INTO 2 PARTS - A 3-WORD CONTINUATION TABLE
!     IN PART 1, AND 21-WORD CONTINUATION CARD IMAGES IN PART 2.
 
!     3-WORD TABLE IN PART 1 HOLDS THE 2-BCD CONTIUATION SYMBOLS, WITH
!     THE FIRST BYTE (A + OR *) ZERO OUT, AND AN INDEX POINTER. THIS
!     TABLE WILL BE SORTED, AND WILL BE USED BY BISLC2 TO LOCATE THE
!     CARD IMAGES SAVED EITHER IN PART 2, OR IN TAPECC FILE.
 
 imhere = 1130
 IF (debug) WRITE (nout,140) imhere,method,ncont,ic,ib
 CALL READ (*3000,*1130,tape2,z(ib),nzib,1,LEN)
 CALL mesage (-8,0,NAME)
 1130 k = LEN + ib - 1
 i = ic
 DO  j = ib,k,21
   z(i  ) = andf(z(j),les1b)
   z(i+1) = z(j+1)
   z(i+2) = j
   i = i + 3
 END DO
 GO TO 1270
 
!     METHOD 3 -
 
!     COMPUTE NCCI (NO. OF CONTINUATION CARD IMAGES) THAT PART 2 AREA
!     (FROM Z(IB) THRU Z(NZ)) CAN HOLD AT A GIVEN TIME.
!     CREATE IN CORE A CONTINUATION TABLE WITH INDEX POINTERS (SAME
!     AS METHOD 2) IN PART 1 AREA.
!     FILL THE REMAINING PART 2 AREA WITH NCCI CARDS, AND WRITE THIS
!     BLOCK OF CARDS OUT TO A NEW SCRATCH FILE, TAPECC. REPEAT THIS
!     PROCESS FOR THE REST OF THE CONTINUATION CARDS.
!     THE INDEX POINTERS IN PART 1 (METHOD 3 ONLY) ALSO INCLUDE THE
!     DATA BLOCK NUMBER INFORMATION
 
 1200 ncci = nzib/21
 IF (ncci >= 10000000) ncci = 10000000 - 1
 nzib   = ncci*21
 tapecc = nfiles + tape3
 imhere = 1200
 IF (debug) WRITE (nout,140) imhere,method,tapecc,ncci
 IF (tapecc > maxscr) GO TO 2951
 CALL OPEN (*2950,tapecc,z(ibufc),wrtrew)
 bk  = 0
 i   = ic
 IF (ncci >= 750 .OR. mach <= 2 .OR. nbpw == 64) GO TO 1220
 j   = ((ncont*23 - nz+ic +999)/1000)*1000
 WRITE  (nout,1210) uim,j,hicore
 1210 FORMAT (a29,', DUE TO UNUSUAL LARGE NUMBER OF CONTINUATION CARDS',  &
     ' PRESENT IN THE BULKDATA DECK', /5X,'AN ADDITION OF',i7,  &
     ' WORDS TO OPEN CORE SPACE COULD MAKE LINK1 MORE EFFICIENT'  &
     ,      /5X,'CURRENTLY NASTRAN HICORE IS',i7,' WORDS')
 IF (ncci < 100) nogo = -3
 1220 bk  = bk + 10000000
 j   = ib
 top = nzib
 CALL READ (*1260,*1230,tape2,z(ib),top,0,LEN)
 GO TO 1240
 1230 top = LEN
 1240 top = top + ib - 1
 1250 z(i  ) = andf(z(j),les1b)
 z(i+1) = z(j+1)
 z(i+2) = j + bk
 i   = i + 3
 j   = j + 21
 IF (j < top) GO TO 1250
 CALL WRITE (tapecc,z(ib),nzib,1)
 GO TO 1220
 1260 CALL CLOSE (tapecc,rew)
 1270 CALL CLOSE (tape2 ,rew)
 LEN = i - ic
 IF (LEN > 3) CALL sort2k (0,0,3,1,z(ic),LEN)
 
!     NO PRE-MERGING FILES IF REDUCE IS 1 (I.E. LESS THAN 10 SCRATCH
!     FILES WERE USED TO HOLD THE RAW BULKDATA, OR ENOUGH CORE TO HOLD
!     EVERYTHING)
 
 1280 IF (reduce == 1) GO TO 1600
 
!     PRE-MERGE
!     =========
 
!     AT THIS POINT, CONTINUATION CARD IMAGES ARE EITHER IN CORE OR IN
!     SCRATCH FILE TAPECC, AND TAPE2 IS FREE FOR RE-USE.
!     ALL GINO BUFFERS ARE FREE
 
!     IF TOO MANY FILES WERE USED TO SAVE BULKDATA, MERGE THEM TO REDUCE
!     THE TOTAL NUMBER OF FILES GOING TO BE USED (I.E. TO REDUCE BUFFER
!     SPACE IN THE MERGE PHASE COMING NEXT)
 
!     PERFORM A 2-TO-1 MERGE IF NUMBER OF FILES PRESENTLY IS 10-17.
 
!     FILEB + FILEC == FILEA      E.G.  303 + 304 == 302
!                                       305 + 306 == 303
!                                       307 + 308 == 304  ETC.
!     OR
!     PERFORM A 3-TO-1 MERGE IF NUMBER OF FILES PRESENTLY IS 18-30.
 
!     FILEB+FILEC+FILED == FILEA  E.G.  303+304+305==302
!                                       306+307+308==303
!                                       309+310+311==304  ETC.
 
!     NOTE - 301 IS EITHER NOT USED, OR USED BY THE 'MODIFIED' OPTP
 
 imhere = 1290
 IF (debug) WRITE (nout,140) imhere,nfiles,nfiler,reduce
 filea  = 301
 FILE   = 302 - reduce
 
 DO  iii = 1,nfiles,reduce
   FILE = FILE+reduce
   
! ... CHECK LAST DO-LOOP CONDITION
!     IF ONE   FILE  LEFT, QUIT MERGING
!     IF TWO   FILES LEFT, DO A 2-TO-1 MERGE
!     IF THREE FILES LEFT, CONTINUE
   
   IF (nfiles-iii <= 0) GO TO 1420
   
   filea = filea + 1
   CALL OPEN (*2930,filea,z(ibuf1),wrtrew)
   imhere= 1300
   exh   = 0
   DO  l = 1,reduce
     filex = FILE + l
     ibufl = ibufx(l)
     itape(l) = 1
     IF (debug) WRITE (nout,140) imhere,filex,j
     CALL OPEN (*2940,filex,z(ibufl),rdrew)
     CALL READ (*3000,*2980,filex,y(1,l),24,0,i)
   END DO
   
!     PICK THE SMALLEST CONTROL WORDS FROM Y(21,22,23,24 OF A,B,C)
   
   1310 ii = 1
   loop1380:  DO  l = 2,reduce
     IF (y(21,l) - y(21,ii) < 0.0) THEN
       GO TO  1370
     ELSE IF (y(21,l) - y(21,ii) == 0.0) THEN
       GO TO  1320
     ELSE
       GO TO  1380
     END IF
     1320 IF (y(21,l) == large) CYCLE loop1380
     IF (y(22,l) - y(22,ii) < 0.0) THEN
       GO TO  1370
     ELSE IF (y(22,l) - y(22,ii) == 0.0) THEN
       GO TO  1330
     ELSE
       GO TO  1380
     END IF
     1330 IF (y(23,l) - y(23,ii) < 0.0) THEN
       GO TO  1370
     ELSE IF (y(23,l) - y(23,ii) == 0.0) THEN
       GO TO  1340
     ELSE
       GO TO  1380
     END IF
     1340 IF (y(24,l) - y(24,ii) < 0.0) THEN
       GO TO  1370
     ELSE IF (y(24,l) - y(24,ii) == 0.0) THEN
       GO TO  1350
     ELSE
       GO TO  1380
     END IF
     
!     FIRST 3 BULKDATA FIELDS THE SAME, CHECK POSSIBLE DUPLICATE CARD
!     SET 21ST AND 22ND CONTROL WORDS TO -6 IF IT IS A DUPLICATE
     
     1350 DO  j = 7,20
       IF (y(j,l) /= y(j,ii)) CYCLE loop1380
     END DO
     y(21,ii) = -6
     y(22,ii) = -6
     nogo = -1
     GO TO 1370
     
     1370 ii = l
   END DO loop1380
   imhere = 1380
   IF (debug) WRITE (nout,140) imhere,ii
   
   IF (y(1,ii) == large) CALL mesage (-61,0,NAME)
   CALL WRITE (filea,y(1,ii),24,0)
   filex = ii + FILE
   CALL READ (*2980,*1400,filex,y(1,ii),24,0,j)
   IF (debug) WRITE (nout,1390) filex,y(1,ii),y(2,ii)
   1390 FORMAT (5X,'TO PRE-MERGE FILE',i5,3X,2A4)
   GO TO 1310
   
! ... ONE OF THE FILES IS EXHAUSTED
   
   1400 exh = exh + 1
   itape(ii) = 0
   IF (exh >= reduce-1) GO TO 1420
   DO  j = 1,24
     y(j,ii) = large
   END DO
   imhere  = 1410
   IF (debug) WRITE (nout,140) imhere,exh
   GO TO 1310
   
! ... ONLY ONE FILE LEFT WHICH HAS NOT BEEN EXHAUSTED
   
   1420 filex = FILE + 1
   IF (itape(2) == 1) filex = FILE + 2
   IF (itape(3) == 1) filex = FILE + 3
   imhere = 1420
   IF (debug) WRITE (nout,140) imhere,filex
   DO  j = 1,24
     z(j) = y(j,filex)
   END DO
   
!     THIS REMAINING FILE COULD BE VERY BIG. IT COULD BE OPTP
   
   left24 = ((left-24)/24)*24
   1440 full = 1
   CALL READ (*3000,*1450,filex,z(i25),left24,0,LEN)
   full = 0
   LEN  = left24
   1450 IF (LEN < 24) GO TO 1560
   
! ... CHECK ANY DUPLICATE IN THIS GROUP, SET THE 21ST AND 22ND CONTROL
!     WORDS TO -6 IF DUPLICATE
!     THEN WRITE THE REST TO FILEA
   
   DO  l = 1,LEN,24
     i = l - 1
     k = i + 24
     DO  j = 21,24
       IF (z(i+j) /= z(k+j)) GO TO 1520
     END DO
     DO  j = 7,20
       IF (z(i+j) /= z(k+j)) GO TO 1520
     END DO
     z(i+21) = -6
     z(i+22) = -6
     1520 CALL WRITE (filea,z(l),24,0)
     IF (debug) WRITE (nout,1530) filea,z(l),z(l+1)
     1530 FORMAT (5X,'TO FILEA',i5,3X,2A4)
   END DO
   
!     IF FILE HAS NOT BEEN EXHAUSTED, GO BACK FOR MORE
   
   IF (full == 1) GO TO 1560
   DO  j = 1,24
     z(j) = z(LEN+j)
   END DO
   GO TO 1440
   
   1560 CALL WRITE (filea,z(LEN+1),24,1)
   IF (debug) WRITE (nout,1530) filea,z(LEN+1),z(LEN+2)
   DO  l = 1,reduce
     filex = FILE + l
     CALL CLOSE (filex,rew)
   END DO
   
   FILE = FILE + reduce
 END DO
 
!     END OF PRE-MERGE
 
 
!     SET UP SORTED HEADING IF APPLICABLE
 
 1600 IF (nbulk <= 1) GO TO 1620
 CALL page2 (2)
 WRITE  (nout,1610) uim
 1610 FORMAT (a29,' 207, BULK DATA DECK IS NOT SORTED. NASTRAN WILL ',  &
     'RE-ORDER THE INPUT DECK.')
 1620 IF (f3long == 0 .OR. echos == 0) GO TO 1640
 CALL page2 (2)
 WRITE  (nout,1630) uim
 1630 FORMAT (a29,' 207A, SIX CHARACTERS OF NASTRAN BCD NAME IN THE ',  &
     'THIRD FIELD WERE USED DURING RE-ORDERING DECK')
 1640 IF (echos == 0) GO TO 1650
 READ (head(2),150) (head1(j),j=11,24)
!WKBR 9/93 HEAD2(4) = CDCNT(1)
 head2(5) = cdcnt(1)
!WKBR 9/93 HEAD3(4) = CDCNT(2)
 head3(5) = cdcnt(2)
!WKBR 9/93 HEAD3(5) = CDCNT(3)
 head3(6) = cdcnt(3)
 CALL page
 
!     FINAL FILE MERGE, ADD CONTINUATION CARD AS NEEDED. RESULTS IN NPTP
!           ==========
 
!     ASSIGN BUFFER SPACES FOR THE SCRATCH FILES, RESERVE IBUF1 FOR NPTP
 
!     OPEN SCRATCH DATA FILES (303,304,305... OR        ==METHODS 1,2==
!     PREVIOUSLY SAVED         303,304,305...301  OR
!                              302,303,304,305... OR    ==METHOD  3  ==
!                              302,303,304,305,...,301)
!     AND READ INTO Y SPACE THE FIRST RECORD OF EACH SCRATCH FILE
 
!     OPEN NPTP FOR MERGED RESULT
 
 
 1650 CALL OPEN  (*3100,nptp,z(ibuf1),wrt)
 CALL WRITE (nptp,bulkda,2,1)
 IF (nbulk+ndele == 0) GO TO 2290
 IF (tapecc /= 0) CALL OPEN (*2950,tapecc,z(ibufc),rd)
 recx   = large
 ncard  = 0
 exh    = 0
 imhere = 1700
 IF (debug) WRITE (nout,140) imhere,ncont,nfiler
 
!     IF NO CONTINUATION CARDS, AND ONLY ONE FILE IS USED TO STORE
!     BULKDATA INPUT CARDS, MOVE DATA FROM TAPE3 (COLD START JOB), OR
!     FROM TAPE1 (RESTART JOB WITH DELETE ONLY AND NO NEW BULK DATA)
!     INTO NPTP DIRECTLY. OTHERWISE, JUMP TO 1760
 
 IF (.NOT.(ncont == 0 .AND. nfiler == 1)) GO TO 1760
 tape = tape3
 IF (restr == 1) tape = tape1
 CALL OPEN (*2920,tape,z(ibuf2),rdrew)
 left24 = ((ibuf2-1)/24)*24
 1700 full = 1
 k    = 1
 CALL READ (*3000,*1710,tape,z(1),left24,0,j)
 full = 0
 j    = left24
 1710 CALL WRITE (nptp,z(k),20,1)
 IF (debug) WRITE (nout,1720) z(k),z(k+1)
 1720 FORMAT (5X,'WRITE TO NPTP',4X,2A4)
 ncard = ncard + 1
 l = k + 19
 IF (echos == 0) GO TO 1740
 CALL page2 (-1)
 WRITE  (nout,1730) ncard,(z(i),i=k,l)
 1730 FORMAT (13X,i8,1H-,8X,20A4)
 1740 IF (echop /= 0) WRITE (lpch,1750) (z(i),i=k,l)
 1750 FORMAT (20A4)
 k = k + 24
 IF (k < j) GO TO 1710
 imhere = 1750
 IF (debug) WRITE (nout,140) imhere,full,j
 IF (full == 0) GO TO 1700
 CALL eof (nptp)
 CALL CLOSE (nptp,rew)
 CALL CLOSE (tape,rew)
 IF (echop /= 0) WRITE (lpch,2320)
 IF (echos == 0) GO TO 2700
 CALL page2 (-1)
 WRITE (nout,2300)
 GO TO 2700
 
!     OPEN AND READ IN THE FIRST DATA RECORD FROM ALL FILES
 
 1760 imhere = 1760
 tape = tape2
 IF (reduce > 1) tape = tape2 - 1
 IF (debug) WRITE (nout,140) imhere,reduce,nfiler,tape
 empty = 0
 DO  ii = 1,nfiler
   tape = tape + 1
   IF (ii == nfiler .AND. restr == 1) tape = tape1
   itape(ii) = tape
   iibuf = ibufx(ii)
   CALL OPEN (*2960,tape,z(iibuf),rdrew)
   CALL READ (*3000,*1780,tape,y(1,ii),24,0,j)
   IF (debug) WRITE (nout,1770) tape,ii,y(1,ii),y(2,ii)
   1770 FORMAT (5X,'SETTING MERGE TABLE.  TAPE,II =',2I4,2X,2A4)
   CYCLE
   1780 empty = empty + 1
   CALL CLOSE (tape,rew)
   DO  i = 1,24
     y(i,ii) = large
   END DO
 END DO
 exh = -1
 DO  ii = 1,nfiler
   IF (y(21,ii) == -6) GO TO 1830
 END DO
 1820 exh = empty
 ii  = 1
 IF (nfiler-1 > 0) THEN
   GO TO  1900
 ELSE
   GO TO  1980
 END IF
 1830 l = ii
 GO TO 2220
 
!     START MERGING FILES
 
!     PICK THE SMALLEST CONTROL WORDS IN 21ST, 22ND, 23RD AND 24TH
!     WORDS OF EACH Y RECORD AND WRITE IT TO MERGE FILE NPTP, 20 WORDS
!     EACH. REPLACE THE CHOSEN RECORD BY NEXT RECORD OF THE SAME FILE
 
 1900 ii = 1
 loop1970:  DO  l = 2,nfiler
   IF (y(21,l) - y(21,ii) < 0.0) THEN
     GO TO  1960
   ELSE IF (y(21,l) - y(21,ii) == 0.0) THEN
     GO TO  1910
   ELSE
     GO TO  1970
   END IF
   1910 IF (y(1,l)  == large) CYCLE loop1970
   IF (y(22,l) - y(22,ii) < 0.0) THEN
     GO TO  1960
   ELSE IF (y(22,l) - y(22,ii) == 0.0) THEN
     GO TO  1920
   ELSE
     GO TO  1970
   END IF
   1920 IF (y(23,l) - y(23,ii) < 0.0) THEN
     GO TO  1960
   ELSE IF (y(23,l) - y(23,ii) == 0.0) THEN
     GO TO  1930
   ELSE
     GO TO  1970
   END IF
   1930 IF (y(24,l) - y(24,ii) < 0.0) THEN
     GO TO  1960
   ELSE IF (y(24,l) - y(24,ii) == 0.0) THEN
     GO TO  1940
   ELSE
     GO TO  1970
   END IF
   
! ... FIRST 3 BULKDATA FIELDS ARE THE SAME, CHECK POSSIBLE DUPLICATE
!     CARDS
   
   1940 DO  j = 7,20
     IF (y(j,ii) /= y(j,l)) CYCLE loop1970
   END DO
   GO TO 2220
   
   1960 ii = l
 END DO loop1970
 
 1980 CALL WRITE (nptp,y(1,ii),20,1)
 ncard = ncard + 1
 IF (echos == 0) GO TO 1990
 CALL page2 (-1)
 WRITE (nout,1730) ncard,(y(j,ii),j=1,20)
 1990 IF (echop /= 0) WRITE (lpch,1750) (y(j,ii),j=1,20)
 IF (ncont == 0) GO TO 2200
 IF (restr == 0) GO TO 2000
 
!     IF THIS IS A RESTART JOB, SAVE THE FIRST FIELD, IN CASE THIS IS
!     THE PARENT OF A CONTINUATION CARD THAT CAME FROM NEW BULK DATA
 
 kard1 = y(1,ii)
 kard2 = y(2,ii)
 IF (kard1 /= param(1) .OR. kard2 /= param(2)) GO TO 2000
 kard1 = y(3,ii)
 kard2 = y(4,ii)
 
!     INSERT CONTINUATION CARD IF NEEDED
 
 2000 IF (nogo ==  -3) GO TO 2200
 tempx   = y(19,ii)
 temp(1) = andf(tempx,les1b)
 temp(2) = y(20,ii)
 2010 IF (tempx == BLANK .AND. temp(2) == BLANK) GO TO 2200
 CALL bislc2 (*2140,temp(1),z(ic),ncont,bsize,loc)
 k = loc*bsize + ic - 1
 l = z(k)
 IF (l < 0) GO TO 2150
 z(k) = -l
 IF (l > 10000000) GO TO 2050
 2020 DO  i = 1,20
   buf(i) = z(l)
   l = l + 1
 END DO
 IF (restr == 0 .OR. kard1 == -1 .OR. z(l) == 0) GO TO 2120
!         ----------     -------------    -----------
!    I.E. NO RESTART     ALREADY DONE     BULKDATA CARD
!                                         NOT FLAGGED
 
!     SET THE PARENT'S RESTART BIT IF ABOVE CONDITIONS NOT MET
 
 ASSIGN 2040 TO crdflg
 from = 2040
 GO TO 2810
 2040 kard1 = -1
 GO TO 2120
 
!     READ IN CONTINUATION CARD IMAGE FROM TAPECC FILE
 
 2050 REC = l/10000000
 l   = l - REC*10000000
 IF (REC-recx < 0.0) THEN
   GO TO  2060
 ELSE IF (REC-recx == 0.0) THEN
   GO TO  2020
 ELSE
   GO TO  2110
 END IF
 2060 CALL REWIND (tapecc)
 IF (REC == 1) GO TO 2090
 skip = REC - 1
 2070 DO  j = 1,skip
   CALL fwdrec (*3020,tapecc)
 END DO
 2090 CALL READ (*3020,*2100,tapecc,z(ib),nzib,1,LEN)
 recx = REC
 GO TO 2020
 2100 CALL mesage (-37,0,NAME)
 2110 skip = REC - recx - 1
 IF (skip < 0.0) THEN
   GO TO  2100
 ELSE IF (skip == 0.0) THEN
   GO TO  2090
 ELSE
   GO TO  2070
 END IF
 
!     GOT THE CONTINUATION CARD, WRITE IT OUT TO NPTP
!     CHECK WHETHER IT ASKS FOR MORE CONTINUATION CARD
 
 2120 CALL WRITE (nptp,buf,20,1)
 ncard = ncard + 1
 IF (echos == 0) GO TO 2130
 CALL page2 (-1)
 WRITE (nout,1730) ncard,(buf(j),j=1,20)
 2130 IF (echop /= 0) WRITE (lpch,1750) (buf(j),j=1,20)
 tempx   = buf(19)
 temp(1) = andf(tempx,les1b)
 temp(2) = buf(20)
 GO TO 2010
 
!     CONTINUATION CARD NOT FOUND. ASSUME THE 10TH FIELD IS USER'S
!     COMMENT
 
 2140 GO TO 2200
 
!     DUPLICATE PARENT - ERROR
 
 2150 CALL page2 (-1)
 IF (echos /= 0) GO TO 2155
 WRITE  (nout,2152) ufm,z(-l),z(-l+1)
 2152 FORMAT (a23,' 208A, ',2A4,' IS DUPLECATE CONTINUATION MARK.')
 GO TO 2180
 2155 WRITE  (nout,2160) ufm
 2160 FORMAT (a23,' 208, PREVIOUS CARD IS A DUPLICATE PARENT.')
 IF (debug) WRITE (nout,2170) loc,bsize,ic,k,l,tempx,temp(2)
 2170 FORMAT ('  LOC,BSIZE,IC,K,L =',5I8,2(2H /,a4),1H/)
 2180 nogo = -1
 
!     REPLACE THE MERGED RECORD BY THE NEXT RECORD OF THE SAME FILE
 
 2200 tape   = itape(ii)
 imhere = 2200
 IF (debug) WRITE (nout,140) imhere,tape,ii
 CALL READ (*3000,*2270,tape,y(1,ii),24,0,j)
 IF (debug) WRITE (nout,2210) tape,ii,y(1,ii),y(2,ii), (y(j,ii),j=21,24)
 2210 FORMAT (5X,'REPLACING - TAPE,II=',2I4,3X,2A4,4I12)
 IF (y(21,ii) /= -6) IF (exh) 1820,1900,1900
 2220 CALL page2 (-2)
 ncard = ncard + 1
 CALL WRITE (nptp,y(1,ii),20,1)
 WRITE  (nout,1730) ncard,(y(j,ii),j=1,20)
 WRITE  (nout,2230) uwm
 2230 FORMAT (a25,' 208, PREVIOUS CARD IS A DUPLICATE')
!     NOGO = -1
 IF (.NOT.debug) GO TO 2200
 DO  k = 1,nfiler
   WRITE  (nout,2240) k,(y(j,k),j=1,24)
   2240 FORMAT (1X,i2,3H)  ,20A4,2H /,4I8)
 END DO
 WRITE  (nout,2260) ii,l
 2260 FORMAT (//5X,'DUPLICATE  II,L=',2I8)
 GO TO 2200
 
!     A SCRATCH FILE IS JUST EXHAUSTED, SET THE CORRESPONDING RECORD
!     A SET OF VERY LARGE NUMBERS
!     IF ALL FILES ARE EXHAUSTED, MERGING DONE
 
 2270 exh = exh + 1
 CALL CLOSE (tape,rew)
 imhere = 2270
 IF (debug) WRITE (nout,140) imhere,tape,exh,nfiler,ncard
 IF (exh >= nfiler) GO TO 2290
 DO  i = 1,24
   y(i,ii) = large
 END DO
 GO TO 1900
 
!     MERGING DONE. EVERY THING IN NPTP.
 
 2290 CALL eof (nptp)
 CALL CLOSE (nptp,rew)
 imhere = 2290
 IF (debug) WRITE (nout,140) imhere,exh,nfiler
 IF (echos == 0) GO TO 2310
 CALL page2 (-1)
 WRITE  (nout,2300)
 2300 FORMAT (30X,'ENDDATA')
 2310 IF (echop /= 0) WRITE (lpch,2320)
 2320 FORMAT ('ENDDATA')
 
!     CHECK AND IDENTIFY PARENTLESS CONTINUATION CARDS
!     MAKE SURE TO EXCLUDE ANY BROKEN CONTINUATION CARDS SUPPOSEDLY
!     CONNECTED TO ONE PARENT
 
 IF (ncont == 0 .OR. nogo == -3) GO TO 2700
 imhere = 2330
 IF (debug) WRITE (nout,140) imhere,ncont,ic
 recx = large
 j = ic + bsize - 1
 DO  i = 1,ncont
   l = z(j)
   2400 IF (l  < 0) GO TO 2490
   imhere = 2400
   IF (debug) WRITE (nout,2480) imhere,z(j-2),z(j-1),l
   IF (l <= 10000000) GO TO 2470
   REC = l/10000000
   l   = l - REC*10000000
   IF (REC-recx < 0.0) THEN
     GO TO  2410
   ELSE IF (REC-recx == 0.0) THEN
     GO TO  2470
   ELSE
     GO TO  2450
   END IF
   2410 CALL REWIND (tapecc)
   IF (REC == 1) GO TO 2440
   skip = REC - 1
   2420 DO  k = 1,skip
     CALL fwdrec (*3020,tapecc)
   END DO
   2440 CALL READ (*3020,*2620,tapecc,z(ib),nzib,1,LEN)
   recx = REC
   GO TO 2470
   2450 skip = REC - recx - 1
   IF (skip < 0.0) THEN
     GO TO  2460
   ELSE IF (skip == 0.0) THEN
     GO TO  2440
   ELSE
     GO TO  2420
   END IF
   2460 CALL mesage (-37,0,NAME)
   2470 temp(1) = andf(z(l+18),les1b)
   temp(2) = z(l+19)
   imhere  = 2470
   IF (debug) WRITE (nout,2480) imhere,temp,l
   2480 FORMAT ('  IMHERE=',i5,'  LOOKING FOR - ',2A4,i14)
   IF (temp(1) == BLANK .AND. temp(2) == BLANK) GO TO 2490
   loc = loc + 1
   IF (temp(1) /= z(loc+ic) .OR. temp(2) /= z(loc*ncont+ic))  &
       CALL bislc2 (*2490,temp(1),z(ic),ncont,bsize,loc)
   k = loc*bsize + ic - 1
   l = z(k)
   z(k) = -IABS(z(k))
   GO TO 2400
   2490 j = j + bsize
 END DO
 
 j  = ic + bsize - 1
 ii = 0
 recx   = large
 imhere = 2600
 DO  i = 1,ncont
   IF (z(j) < 0) GO TO 2610
   IF (ii   == 1) GO TO 2510
   ii = 1
   CALL page1
   WRITE  (nout,2500) ufm
   2500 FORMAT (a23,' 209, THE FOLLOWING CONTINUATION INPUT CARDS HAVE ',  &
       'NO PARENTS',//)
   nogo = -1
   2510 CALL page2 (1)
   l = z(j)
   IF (l > 10000000) GO TO 2540
   2520 m = l + 19
   WRITE  (nout,2530) (z(k),k=l,m)
   2530 FORMAT (10X,20A4)
   GO TO 2610
   
   2540 REC = l/10000000
   l   = l - REC*10000000
   IF (REC-recx < 0.0) THEN
     GO TO  2550
   ELSE IF (REC-recx == 0.0) THEN
     GO TO  2520
   ELSE
     GO TO  2600
   END IF
   2550 CALL REWIND (tapecc)
   IF (REC == 1) GO TO 2580
   skip = REC - 1
   2560 DO  k = 1,skip
     CALL fwdrec (*3020,tapecc)
   END DO
   2580 CALL READ (*3020,*2620,tapecc,z(ib),nzib,1,LEN)
   recx = REC
   GO TO 2520
   2600 skip = REC - recx - 1
   IF (skip < 0.0) THEN
     GO TO  2620
   ELSE IF (skip == 0.0) THEN
     GO TO  2580
   ELSE
     GO TO  2560
   END IF
   2610 j = j + bsize
 END DO
 GO TO 2700
 2620 CALL mesage (-2,tapecc,NAME)
 
!     CLOSE CONTINUAION CARD FILE TAPECC, IF IT WAS OPENED
!     DISABLE FREE-FIELD INPUT OPTION IN XREAD.
 
 2700 IF (tapecc > 0) CALL CLOSE (tapecc,rew)
 ffflag = 0
 wasff  = 0
 IF (nogo /= -3) GO TO 2730
 WRITE  (nout,2710) ufm
 2710 FORMAT (a23,' 3008, CONTINUATION CARDS WERE NOT ADDED TO SORTED ',  &
     'BULKDATA DECK DUE TO INSUFFICIENT CORE CONDITION.')
 IF (cpflg /= 0) WRITE (nout,2720)
 2720 FORMAT (5X,'THE NPTP FILE OR TAPE GENERATED IN THIS RUN IS NOT ',  &
     'SUITABLE FOR RESTART')
 CALL mesage (-61,0,0)
 2730 IF (nogo /= 0) nogo = 1
 IF (.NOT. debug) GO TO 3200
 
!     DEBUG NPTP ECHO
 
 imhere = 2730
 WRITE (nout,140) imhere,ffflag,wasff
 CALL OPEN (*3100,nptp,z(ibuf1),rdrew)
 2740 CALL skpfil (nptp,+1)
 CALL READ (*2770,*2770,nptp,buf(1),2,1,j)
 IF (buf(1) /= bulkda(1) .OR. buf(2) /= bulkda(2)) GO TO 2740
 2750 CALL READ (*2770,*2770,nptp,buf(1),20,1,j)
 WRITE  (nout,2760) (buf(j),j=1,10),(buf(j),j=17,20)
 2760 FORMAT (' ==NPTP==>',5(1X,2A4),'...',2(1X,2A4))
 GO TO 2750
 2770 CALL CLOSE (nptp,rew)
 GO TO 3200
 
 
!     INTERNAL ROUTINE TO SET RESTART BITS - CRDFLG
 
!     BITS SET ONLY IF JOB IS A RESTART RUN, AND
!       1. ALL NEW BULK DATA CARDS,   EXCEPT CONTINUATION CARDS
!       2. ALL DELETED CARDS IN OPTP, EXCEPT CONTINUATION CARDS
!       3. THE PARENTS OF THE CONTINUATION CARDS IN 1 AND 2
 
 2800 kard1 = buf(1)
 kard2 = buf(2)
 IF (kard1 /= param(1) .OR. kard2 /= param(2)) GO TO 2810
 kard1 = buf(3)
 kard2 = buf(4)
 2810 imhere = 2810
 IF (debug) WRITE (nout,2820) imhere,from,nogo,kard1,kard2
 2820 FORMAT (/,' *** IMHERE',i5,', FROM',i5,', NOGO=',i3,3X,2A4)
 IF (nogo /= 0) GO TO 2850
 k = numx1*2
 DO  i = 1,k,2
   IF (kard1 /= icards(i) .OR. kard2 /= icards(i+1)) CYCLE
   j = i/2
   m = (j/31) + 1
   n = MOD(j,31) + 2
   ibits(m) = orf(ibits(m),itwo(n))
   IF (debug) WRITE (nout,2830) kard1,kard2
   2830 FORMAT (5X,'BITS SET SUCCESSFULLY FOR ',2A4)
   EXIT
 END DO
 2850 GO TO crdflg, (200,810,880,2040)
 
!     ERRORS
 
 2900 tape = tape1
 GO TO  2960
 2910 tape = tape2
 GO TO  2960
 2920 tape = tape3
 GO TO  2960
 2930 tape = filea
 GO TO  2960
 2940 tape = filex
 GO TO  2960
 2950 tape = tapecc
 IF (tapecc <= maxscr) GO TO 2960
 2951 WRITE  (nout,2955) sfm
 2955 FORMAT (a25,' 212, NUMBER OF AVAILABLE SCRATCH FILES EXEEDED.',5X,  &
     'RE-RUN JOB WITH MORE CORE')
 GO TO  3140
 2960 WRITE  (nout,2970) sfm,tape
 2970 FORMAT (a25,' 210, COULD NOT OPEN SCRATCH FILE',i5)
 GO TO  3140
 2980 WRITE  (nout,2990) sfm
 2990 FORMAT (a25,' 211, ILLEGAL EOR ON SCRATCH')
 GO TO  3140
 3000 WRITE  (nout,3010) sfm,tape
 3010 FORMAT (a25,' 212, ILLEGAL EOF ON SCRATCH',i5)
 GO TO  3140
 3020 WRITE  (nout,3030)
 3030 FORMAT (//26X,'212, TAPECC ERROR')
 tape = tapecc
 GO TO  3000
 3040 WRITE  (nout,3050) sfm
 3050 FORMAT (a25,' 213, ILLEGAL EOF ON OPTP')
 GO TO  3140
 3060 WRITE  (nout,3070) sfm,imhere
 3070 FORMAT (a25,' 213X, ILLEGAL DATA ON OPTP.  IMHERE =',i7)
 nogo = 1
 GO TO  810
 3080 WRITE  (nout,3090) sfm
 3090 FORMAT (a25,' 214, OPTP COULD NOT BE OPENED')
 GO TO  3140
 3100 WRITE  (nout,3110) sfm
 3110 FORMAT (a25,' 215, NPTP COULD NOT BE OPENED')
 GO TO  3140
 3120 WRITE  (nout,3130) sfm,imhere
 3130 FORMAT (a25,' 219, MISSING ENDDATA CARD.  IMHERE =',i7)
 nogo = 1
 GO TO  350
 3140 WRITE  (nout,3150) imhere
 3150 FORMAT (5X,'IMHERE =',i6)
 CALL mesage (-37,0,NAME)
 
!     TURN OFF XSORT FLAG AND FREE-FIELD FLAG
 
 3200 ixsort = 0
 RETURN
END SUBROUTINE xsort2
