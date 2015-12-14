SUBROUTINE nsinfo (jump)
     
!     THIS ROUTINE READS AND PROCESSES DATA IN THE NASINFO FILE
 
!     JUMP = 2, NSINFO IS CALLED BY NASCAR TO OPEN NASINFO FILE AND
!               PROCESS THE SYSTEM PRESET PARAMETERS IN THE 2ND SECTION
!               OF THE FILE, AND THE BCD WORDS (USED ONLY BY NUMTYP
!               SUBROUTINE) IN THE 3RD SECTION
!     JUMP = 3, NSINFO IS CALLED BY TTLPGE TO PROCESS THE INSTALLATION-
!               CENTER-TO-USER MESSAGES STORED IN THE 4TH SECTION OF
!               THE NASINFO FILE
!     JUMP = 4, NSINFO IS CALLED BY XCSA TO ECHO DIAG 48 MESSAGE STORED
!               IN THE 5TH SECITON OF THE NASINFO FILE.
 
!     SINCE DIAG48 MAY NOT BE CALLED, NASINFO FILE IS CLOSED BY XCSA
 
!     WRITTEN BY G.CHAN/UNISYS    6/1990
 
 
 INTEGER, INTENT(OUT)                     :: jump
 IMPLICIT INTEGER (a-z)
!WKBR 8/94 SUN INTEGER         NAME(2),NTAB(5),CARDX(4),CARD(20),DIAG48(4)
 INTEGER :: NAME(2),cardx(4),card(20),diag48(4)
 REAL :: time
!WKBR CHARACTER*167   IFILE
 CHARACTER (LEN=144) :: ifile
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /machin/ mach
 COMMON /system/ sys(100)
 COMMON /ntime / nt,time(16)
 COMMON /output/ dum(64),pghdg3(32)
 COMMON /numtpx/ nbcd,bcd(1)
 COMMON /BLANK / iblnk(60)
 EQUIVALENCE     (cardx(1),card(1))
 EQUIVALENCE     (sys(1),sysbuf), (sys( 2),nout), (sys( 9),nlpp  ),  &
     (sys(14),mxlns), (sys(19),ehco), (sys(31),hicore),  &
     (sys(35),lprus), (sys(37),lu  ), (sys(36),nprus ),  &
     (sys(20),pltop), (sys(92),dict), (sys(76),nosbe ),  &
     (sys(77),bndit), (sys(91),lpch)
 DATA    ttpg     / 0            /
!WKBR 8/94 SUN DATA    EQU   , R  ,  S  ,  BNK      ,   EQUALS   , NAME         /
!WKBR 8/94 SUN1        1H=   , 1HR,  1HS,  4H       ,   4H====   , 4HNSIN,2HFO  /
 DATA    equ   , s  ,  bnk      ,   equals   , NAME         /  &
     1H=   , 1HS,  4H       ,   4H====   , 4HNSIN,2HFO  /
 
 DATA    relse , tpg,  pop,  tim,  mxl,   bsz      , s3s   ,skp3  /  &
     4HELEA, 3HTPG,3HPOP,3HTIM,3HMXL, 3HBSZ    , 3HS3S ,0     /
DATA    lpp   , hic,  bnd,  ech,  nos,   pru,  npr,  pch,  END   /  &
3HLPP , 3HHIC,3HBND,3HECH,3HNOS, 3HPRU,3HNPR,3HPCH,3HEND /
DATA    s88   , s89,  s90,  s92,  s94,   s96,  s97,  s98,  s99   /  &
    3HS88 , 3HS89,3HS90,3HS92,3HS94, 3HS96,3HS97,3HS98,3HS99 /
DATA    diag48                         , dd   ,dic  ,cod  ,key   /  &
    4H d i, 4H a g, 4H   4, 2H 8   , 3H$. ,3HDIC,3HCOD,3HKEY /

SELECT CASE ( jump )
  CASE (    1)
    GO TO 550
  CASE (    2)
    GO TO 200
  CASE (    3)
    GO TO 350
  CASE (    4)
    GO TO 400
END SELECT

!     JUMP = 2
!     ========

!     OPEN NASINFO FILE, AND SET LU, THE 37TH WORD OF /SYSTEM/

!     CURRENTLY 'NASINFO' IS USED FOR ALL MACHINES OF TYPE 5 AND HIGHER

200  lu  = 99
CALL nasopn (*280, lu, ifile)

!     SEARCH FOR FIRST EQUAL-LINE

210  READ (lu,220,ERR=275,END=275) cardx
220  FORMAT (20A4)
IF (card(1) /= equals .AND. card(2) /= equals) GO TO 210

!     READ AND PROCESS THE 2ND SECTION OF NASINFO FILE

230  READ (lu,235,END=500) symbol,EQ,value
235  FORMAT (a4,a1,i7)
IF (symbol == bnk) GO TO 230
IF (EQ     /= equ) GO TO 520
IF (symbol == tim) GO TO 250
IF (symbol == END) GO TO 290
IF (value  == -99) GO TO 230
IF (symbol == s3s) GO TO 240
IF (symbol == bsz) sysbuf  = value
IF (symbol == lpp) nlpp    = value
IF (symbol == hic) hicore  = value
IF (symbol == mxl) mxlns   = value
IF (symbol == tpg) ttpg    = value
IF (symbol == ech) echo    = value
IF (symbol == pch) lpch    = value
IF (symbol == dic) dict    = value
IF (symbol == bnd) bndit   = value
IF (symbol == pop) pltop   = value
IF (symbol == pru) lprus   = value
IF (symbol == npr) nprus   = value
IF (symbol == nos) nosbe   = value
IF (symbol == cod) code    = value
IF (symbol == key) key     = value
symb1 = khrfn1(bnk,1,symbol,1)
IF (symb1  /=   s) GO TO 230
IF (symbol == s88) sys(88) = value
IF (symbol == s89) sys(89) = value
IF (symbol == s90) sys(90) = value
IF (symbol == s92) sys(92) = value
IF (symbol == s94) sys(94) = value
IF (symbol == s96) sys(96) = value
IF (symbol == s97) sys(97) = value
IF (symbol == s98) sys(98) = value
IF (symbol == s99) sys(99) = value
GO TO 230

!     SKIP JUMP 3 PRINTOUT

240  skp3 = 1
GO TO 230

!     READ IN 16 GINO TIME CONSTANTS (NT=16)

250  IF (value /= nt) GO TO 270
READ (lu,260,END=500) time
260  FORMAT (12X,8F7.2, /12X,8F7.2)
GO TO 230
270  READ (lu,235,END=500) symbol
READ (lu,235,END=500) symbol
GO TO 230

!     NASINFO DOES NOT EXIST (or IS WRITE-PROTECTED), SET LU TO ZERO

275  CLOSE (UNIT=lu)
CALL mesage (2,0,NAME)
280  WRITE (nout, 285) ifile
285  FORMAT ('0*** USER WARNING MESSAGE, UNABLE TO OPEN ',  &
    'THE FOLLOWING NASINFO FILE -- '//
!WKBR*        1X, A167/)  &
1X, a44/)
lu = 0
GO TO 550

!     READ PASS THE 2ND EQUAL-LINE. CONTINUE INTO 3RD SECTION

290  READ (lu,220,END=500) cardx
IF (card(1) /= equals .AND. card(2) /= equals) GO TO 290


!     THIS 3RD SECTION CONTAINS BCD WORDS WHICH ARE REALLY REAL NUMBERS.
!     (THE BINARY REPRESENTATIONS OF SOME REAL NUMBERS AND THEIR
!     CORRESPONDING BCD WORDS ARE EXACTLY THE SAME. SUBROUTINE NUMTYP
!     MAY IDENTIFY THEM AS TYPE BCD. ANY WORD ON THE BCD LIST WILL BE
!     REVERTED BACK TO AS TYPE REAL. THE LIST IS MACHINE DEPENDENT)

!     SKIP FIRST 5 COMMENT LINES

READ (lu,220,END=500)
READ (lu,220)
READ (lu,220)
READ (lu,220)
READ (lu,220)

300  READ (lu,305,END=500) machx,nbcd
305  FORMAT (i2,i3)
IF (machx == mach) GO TO 320
IF (nbcd == 0) GO TO 300
DO  i = 1,nbcd,19
  READ (lu,325)
END DO
GO TO 300
320  IF (nbcd == 0) GO TO 340
jb = 1
DO  i = 1,nbcd,19
  je = jb + 18
  READ (lu,325) (bcd(j),j=jb,je)
  325  FORMAT (5X,19(a4,1X))
END DO

!     READ PASS THE 3RD EQUAL-LINE, THEN RETURN

340  READ (lu,220,END=500) cardx
IF (card(1) /= equals .AND. card(2) /= equals) GO TO 340
IF (ttpg /= 0) jump = ttpg
GO TO 550

!     JUMP = 3
!     ========

!     READ AND ECHO OUT INSTALLATION-CENTER-TO-USER MESSAGES, SAVED IN
!     THE 4TH SECTION OF NASINFO FILE
!     TERMINATE MESSAGES BY THE LAST EQUAL-LINE.

!     IN THIS MESSAGE SECTION ONLY, SKIP INPUT LINE IF A '$.  ' SYMBOL
!     IS IN FIRST 4 COLUMNS.

350  IF (lu == 0 .OR. skp3 == 1) GO TO 550
CALL page1
360  READ (lu,220,END=500) card
IF (card(1) ==     dd) GO TO 360
IF (card(1) /= equals) GO TO 380
IF (card(2) == equals) GO TO 550
CALL page1
WRITE  (nout,370)
370  FORMAT (//)
GO TO 360
380  WRITE  (nout,390) card
390  FORMAT (25X,20A4)
GO TO 360

!     JUMP = 4
!     ========

!     PROCESS DIAG48 MESSAGE, SAVED IN THE 5TH SECTION OF NASINFO FILE

400  CALL sswtch (20,l20)
IF (lu == 0) GO TO 480
DO  i = 10,20
  pghdg3(i) = bnk
END DO
pghdg3(6) = diag48(1)
pghdg3(7) = diag48(2)
pghdg3(8) = diag48(3)
pghdg3(9) = diag48(4)
line  = nlpp + 1
count = 0

!     READ AND PRINT RELEASE NEWS
!     PRINT LAST TWO YEARS OF NEWS ONLY, IF DIAG 20 IS ON
!     (MECHANISM - GEAR TO THE 'nn RELEASE' LINE AND '========' LINES)

one = 1
IF (l20 == 1) one = 0
420  READ (lu,220,END=540) card
IF (card(1) == equals .AND. card(2) == equals) GO TO 470
IF (one == -1) GO TO 440
IF (card(2) /= relse  .OR.  card(4) /=   bnk) GO TO 440
count = count + one
IF (count <= 2) GO TO 440
430  READ (lu,220,END=540) cardx
IF (card(1) /= equals .OR. card(2) /= equals) GO TO 430
GO TO 470
440  IF (line < nlpp) GO TO 460
CALL page1
IF (line == nlpp) GO TO 450
line = 3
GO TO 460
450  WRITE (nout,370)
line = 5
460  line = line + 1
WRITE (nout,390) card
GO TO 420

!     READ AND PRINT THE REST OF SECTION 5

470  IF (one == -1) GO TO 540
one  = -1
line = nlpp + 1
GO TO 420

480  WRITE  (nout,490) uim
490  FORMAT (a29,', DIAG48 MESSAGES ARE NOT AVAILABLE DUE TO ABSENCE ',  &
    'OF THE NASINFO FILE')
GO TO 540

!     ERROR

500  WRITE  (nout,510) sfm
510  FORMAT (a25,' 3002, EOF ENCOUNTERED WHILE READING NASINFO FILE')
STOP 'JOB TERMINATED IN SUBROUTINE NSINFO'
520  WRITE  (nout,530) symbol,EQ,value
530  FORMAT ('0*** ERROR IN NASINFO FILE - LINE - ',a4,a1,i7)
GO TO 230

540  IF (l20 == 0) GO TO 550
CLOSE (UNIT=lu)
CALL pexit
550  RETURN
END SUBROUTINE nsinfo
