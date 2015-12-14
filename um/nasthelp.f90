PROGRAM nasthelp
 
!     THIS PROGRAM PROVIDES ON-LINE SCREEN HELP FOR NASTRAN USER'S
!     MANUAL INFORMATION. THE COMPLETE MANUAL IS STORED IN THE
!     FOLLOWING ASCII TEXT FILES, WHICH ARE ALSO ACCESSIBLE TO ANY
!     SYSTEM EDITOR:
 
!          MANUAL SECTION                                FILE NAME
!          -------------------                           ---------
!       1. EXECUTIVE CONTROL                             EXEC.TXT
!       2. CASE CONTROL                                  CASE.TXT
!       3. INPUT BULK DATA                               BULK.TXT
!       4. PLOTTING                                      PLOT.TXT
!       5. DMAP                                          DMAP.TXT
!       6. SUBSTRUCTURE                                  SUBS.TXT
!       7. ERROR MESSAGES                                MSSG.TXT
!       8. NASTRAN DICTIONARY                            DICT.TXT
!       9. INTRODUCTION & GENERAL INFORMATION            INTR.TXT
!      10. USER'S MASTER FILE AND USER GENERATED INPUT   UMFL.TXT
!      11. RIGID FORMATS                                 RFMT.TXT
 
 
!     IN ADDITION,
!     IF FILE SELECTION IS FOLLOWED BY ',C', THIS PROGRAM WILL ONLY
!     CHECK THE INPUT FILE FOR OUT-OF-ORDER ITEMS, DUPLICATE ITEMS,
!     AND HEADER 12 CHARACTERS (EXEC, CASE, BULK, PLOT, DMAP AND
!     SUBS.TXT FILES ONLY)
 
!     IF FILE SELECTION IS FOLLOWED BY ',P', THIS PROGRAM WILL PRINT
!     THE ENTIRE CONTENTS OF THE FILE WITH PROPER CARRIAGE CONTROL AND
!     PAGING
 
!     DESIGN REQUIREMENTS FOR MANUAL TEXT FILES
!     (1) A POUND SIGN (#) ON COLUMN 1, MUST PRECEED EACH ITEM
!     (2) '=PAGE=' IN FIRST 6 COLUMNS OF A LINE IS A PAGE MARK
!     (3) EACH ITEM MUST BEGIN WITH ONE OF THE FOLLOWING WORDS, 12 CHAR.
!         EACH               111111
!                 ..123456789012345..(COLUMN)
!                   Executive Co
!                   Case Control
!                   Input Data C
!                   Structure Pl
!                   X-Y Output D
!                   Name:
!                   Substructure Co (special
!                   Substructure Mo  15
!                   Substructure Op  chars.)
!         (REVISED 4/93, DUE TO CHANGES IN THE .TXT FILES, HEADER WORDS
!         IN (3) ARE NO LONGER USED)
!     (4) SEARCH BY KEY NASTRAN WORD AFTER EACH HEADER WORDS IN (3)
!     (5) KEY NASTRAN WORDS MUST BE IN ALPHA-NUMERIC SORT
!     (6) USE ',C' OPTION IN FILE SELECTION FOR DATA CARDS CHECK
 
!     SOME OF THE ABOVE COMMENTS MAY NO LONGER BE TRUE (1993)
 
!     FORTRAN FILE ASSIGNMENTS -
 
!                  FORTRAN
!     FILE NAME    UNIT NO.   STATUS            FILE CONTENTS
!     -----------  -------   ---------  --------------------------------
!     SYS$INPUT        5      INPUT     KEYBOARD INPUT
!     SYS$OUTPUT       6      OUTPUT    TERMINAL OUTPUT
!     COV.TAB          4      INPUT     CONTAINS MACHINE DEPENDENT ASCII
!                           (OPTIONAL)  SPECIAL SYMBOL CONVERSION TABLE
!     NASTRAN MANUAL   3      INPUT     EXEC,CASE,BULK,PLOT,SUBS,DMAP,
!        .TXT FILES                     MSSG,DICT,INTR UMFL and RFMT.TXT
!     USER GIVEN       2      OUTPUT    OUTPUT PRINT FILE
!      FILE NAME            (OPTIONAL)
 
!     PROGRAM FLAGS USED:
!     FL    = 1  THRU 11, FOR 11 DIFFERENT MANUAL.TXT FILES
!     SEC   = 1, MEANS SEARCH BY SECTION ALLOWED, ZERO OTHERWISE
!     HD12  = 0, NO HEADER LINE ON TEXT. OTHEREWISE,
!                HD12(FL) IS THE APPROPIATE HEADER LINE FOR FILE FL
!     BASE  = N, SKIP N WORDS ON HEADER LINE WHEN SEARCHING KEY WORD
!           = 0, SET TO ZERO DUE TO CHANGES IN 1993 USER'S MANUAL
!     MIDPT = AN INTEGER OF A CHARACTER SYMBOL INDICATING THE MID POINT
!                ON TEXT FILE
!     MDPT  = M, NO. OF LINES TO SKIP TO MID POINT OF TEXT FILE
!           = 0, MEANS NO SKIPPING
!     J4    =    2ND ALTERNATE BASE FOR PLOT.TXT
!     J5,J6 =    2ND AND 3RD ALTERNATE BASE FOR SUBS.TXT
 
!     WRITTEN BY GORDON CHAN/UNISYS       3/1992
!     REVISED FOR NEW .TXT FILES FORMAT   4/1993
!     General cleanup and comments
!     added by Reg Mitchell, GSFC         8/1994
!     Last modification                   9/9/94
 
 IMPLICIT INTEGER (a-z)
 LOGICAL :: check,PRINT,first,ety,popen,debug
 INTEGER :: ival(256),numsub(256),nval(256,10),base(11)
 CHARACTER (LEN=1) :: ka,kb,kc,kd,ke,kf,kg,kh,ki,kj,kk,kl,km,kn,ko,kp,kq,    &
                      kr,ks,kt,ku,kv,kw,kx,ky,kz,bnk1,lb1,up1,lt1,mns,plus,  &
                      cma,dot,num,qm, lla,llz, lc,ic,jc,jx1,yes,no, a1(80),  &
                      a11,k1(8)
 CHARACTER (LEN=3) :: pag3, a3
 CHARACTER (LEN=4) :: a4,key4,key42,bnk4,new4,was4,appr,app,solu,sol,        &
                      stp4,stp4l,exit4,exit4l,quit4,quit4l,nxtp,nxts,nxtb,   &
                      file,filex(11)
 CHARACTER (LEN=6) :: key6, a6, pag6, hd6
 CHARACTER (LEN=7) :: dev_dir
 CHARACTER (LEN=8) :: fname, key8, dbgo, dbgf
 CHARACTER (LEN=9) :: date9
 CHARACTER (LEN=12) :: a12, hd12(6)
 CHARACTER (LEN=16) :: mach
 CHARACTER (LEN=35) :: tctf
 CHARACTER (LEN=42) :: ou12
 CHARACTER (LEN=44) :: k44, ka44
 CHARACTER (LEN=48) :: a48, b48
 CHARACTER (LEN=79) :: a79
 CHARACTER (LEN=80) :: a80
 
 COMMON /khr/ ka,kb,kc,kd,ke,kf,kg,kh,ki,kj,kk,kl,km,kn,ko,kp,kq,  &
              kr,ks,kt,ku,kv,kw,kx,ky,kz,bnk1,lb1,up1,lt1,mns,plus,&
              cma,dot,num(10)
              
 EQUIVALENCE  (ka44,ka),(yes,ky),(no,kn),(pag3,pag6),(FILE,fname),  &
              (a1(1),a11,a4,a3,a6,a12,a48,a79,a80,jx1),(a1(3),jc),  &
              (k1(1),key4,key6,key8),(key42,k1(5)),(a1(2),ic)
              
 DATA    in,out,tb,ou / 5, 6, 4, 2  /, nlp,ety / 21, .true. /,  &
         base  /  0,  0,  0,  0, 0,  0, 0, 0, 0, 0, 0       /,  &
         j4,j5,j6, t1,       t2,    t3 , lla, llz, new4     /   &
         21,26,41, 16777216, 65536, 256, 'a', 'z', '|   '   /,  &
         k44 /'ABCDEFGHIJKLMNOPQRSTUVWXYZ #^<-+,.1234567890'/,  &
         bnk4,stp4,pag6,b48  / '   ', 'STOP','=PAGE=', ' '  /,  &
         nxtp,nxts,nxtb,qm   / 'P  ', 'S   ','B   ', '?'    /,  &
         appr,app, solu,sol  / 'APPR','APP ','SOLU','SOL '  /,  &
         filex / 'EXEC','CASE','BULK','PLOT','DMAP','SUBS'  ,   &
         'MSSG','DICT','INTR','UMFL','RFMT'/,fname/'XXXX.TXT'   &
         /,           mach / '  UNIX VERSION  '  /,  date9 / 'AUG. 1994' /,  &
         first,check,print,debug    / .true.,  3*.false.    /,  &
         lu,ivff / 3, 12 /,     hd6 / 'Name: '/,   hd12     /   &
         'Executive Co', 'Case Control', 'Input Data C'     ,   &
         'Structure Pl',  'X-Y Output D', 'Substructure'    /,  &
         tctf / ' or terminate current text file(^):'       /,  &
         dbgo,dbgf / 'DEBUG ON', 'DEBUG OF'/,stp4l /'stop'/,    &
         exit4,exit4l / 'EXIT','exit' /, quit4,quit4l / 'QUIT','quit' /
 DATA dev_dir/'DEV_DIR'/
 
 complf(i) = NOT(i)
 
 ka44   = k44
 num1   = ICHAR(num(1))
 num9   = ICHAR(num(9))
 la     = ICHAR(lla)
 lz     = ICHAR(llz)
 ba     = ICHAR(ka )
 aa     = ba - la
 numcov = 0
 popen  = .false.
 
!     OPEN THE SPECIAL CHARACTER CONVERSION FILE (UNIT 2).
 
 OPEN (UNIT=tb,FILE='COV.TAB',ACCESS='SEQUENTIAL',FORM='FORMATTED',  &
     STATUS='OLD',ERR=130)
 
!     COV.TAB FILE BEGINS WITH A HEADER RECORD, THEN FOLLOWED BY RECORDS
!     OF 3 INTEGER WORDS, IN 3I4 FORMAT, WHICH ARE:
!     INCOMING SYMBOL, NO. OF BYTE, AND OUTGOING CORRESPONDING SYMBOL
 
!        1   1 124
!        2   1 124
!        :   :  :
!      256   1 124
 
 READ (tb,100,END=130)
 DO  i = 1,256
   READ (tb,100,END=130) ival(i),n,(nval(i,j),j=1,n)
   IF (n > 10) WRITE (out,110) i,n
   100 FORMAT (12I4)
   110 FORMAT (' *** Error in COV.TAB    I,N =',2I5)
   numsub(i) = n
 END DO
 i = 257
 CLOSE (UNIT=tb)
!     SPECIAL CHARACTER TABLE HAS BEEN READ
 130 numcov = i - 1
 
!     PRINT PROGRAM HEADER
 
 WRITE  (out,140)
 140 FORMAT (/////////)
 WRITE  (out,150) mach,date9
 150 FORMAT (34X,4H****, /32X,1H*,6X,1H*, /31X,1H*,8X,1H*, /31X,  &
     18H*  n a s t h e l p, /31X,1H*,8X,1H*, /32X,1H*,6X,1H*,  &
     /34X,4H****, ///15X,a16,10X,17HSYSTEM release - ,a9)
 WRITE  (out,160)
 160 FORMAT (//,' Is your screen capable of MORE THAN 80 columns? ',  &
     '(Y or N (default))')
 READ   (in,170) lc
 170 FORMAT (3A1)
 ety = .false.
 IF (lc == bnk1) GO TO 180
 IF (ICHAR(lc) >= la .AND. ICHAR(lc) <= lz) lc = CHAR(ICHAR(lc)+aa)
 IF (lc == yes) ety = .true.
 180 WRITE  (out,190) nlp
 190 FORMAT (//,' Enter LINES PER PAGE (default is',i3,') ')
 READ   (in,200,ERR=180) j
 200 FORMAT (i2)
 IF (j > 1) nlp = j
 nlp4 = nlp - 4
 WRITE  (out,210)
 210 FORMAT (/,' A text line marked by | in column 1 indicates that ',  &
     'this line contains updated', /,' material since the ',  &
     'June 1986 NASTRAN Users'' MANUAL')
 IF (ety) WRITE (out,215)
 215 FORMAT (/,' If your screen loses character in column 80, it is ',  &
     'because your terminal lacks',/,' 80-column capability.')
 WRITE  (out,220)
 220 FORMAT (//,' NASTHELP accepts both UPPER and lower case input')
!    1         /,' The terms STOP and QUIT are interchangable')
!     READY TO BEGIN READING A MANUAL
 GO TO 240
!     CLOSE CURRENT MANUAL IF STILL OPEN
 230 CLOSE (UNIT=lu)
 
!     PROCESS REQUEST FOR USER MANUAL SELECTION
 
 240 WRITE  (out,250)
 250 FORMAT (/,' Enter letter for desired part of NASTRAN',  &
     ' User''S MANUAL',//, '  Introduction(I)     Don''T KNOW(?)',/,  &
     '  Executive(E)        Case control(C)     Bulkdata(B)',/,  &
     '  DMAP(D)             Rigid Formats(R)    Plotting(P)',/,  &
     '  Messages(M)         Substructures(S)    UMF/UGI(U)',/,  &
     '  Dictionary(T)       Stop/quit(STOP)')
 READ (in,170) lc,ic,jc
!     CONVERT TO UPPER CASE IF NECESSARY
 IF (ICHAR(lc) >= la .AND. ICHAR(lc) <= lz) lc = CHAR(ICHAR(lc)+aa)
 IF (ICHAR(ic) >= la .AND. ICHAR(ic) <= lz) ic = CHAR(ICHAR(ic)+aa)
 IF (ICHAR(jc) >= la .AND. ICHAR(jc) <= lz) jc = CHAR(ICHAR(jc)+aa)
!     CHECK FOR A COMMA IN SECOND POSITION
 IF (ic /= cma) GO TO 260
 IF (jc ==  kc) check = .true.
 IF (jc ==  kp) PRINT = .true.
 260 IF (lc == ks .AND. jc == ko) GO TO 2210
 first = .true.
 last  = -1
 sec   =  1
 mdpt  =  0
 mq    =  0
 IF (lc == ke) GO TO 310
 IF (lc == kc) GO TO 320
 IF (lc == kb) GO TO 330
 IF (lc == kp) GO TO 340
 IF (lc == kd) GO TO 350
 IF (lc == ks) GO TO 360
 IF (lc == km) GO TO 370
 IF (lc == kt) GO TO 380
 IF (lc == ki) GO TO 390
 IF (lc == ku) GO TO 400
 IF (lc == kr) GO TO 405
 IF (lc == qm) GO TO 280
 WRITE  (out,270)
 270 FORMAT (/,' *** SELECTION error')
 GO TO 240
 
 280 mq = 1
 GO TO 310
 
!     GET NEW MANUAL REQUEST
 300 CLOSE (UNIT=lu)
 GO TO 240
 
!     SET FL TO THE REQUESTED NASTRAN MANUAL (EXEC=1, CASE=2, ETC.)
!     DEFINE MID-POINT IN FILE FOR SKIPPING, AND SET HEADER SEARCH
!     MIDPT IS THE FIRST LETTER OF KEY WORDS,
!     MDPT  IS NO. OF RECORDS TO BE SKIPPED.
 
!     USE ANY SYSTEM EDITOR TO LOCATE THE MID-POINT OF FILE
 
 310 fl   = 1
 GO TO  410
 320 fl   = 2
 GO TO  410
 330 fl   = 3
 midpt= ICHAR(kn)
 mdpt = 10741
 GO TO  410
 340 fl   = 4
 GO TO  410
 350 fl   = 5
 GO TO  410
 360 fl   = 6
 GO TO  410
 370 fl   = 7
 sec  = 0
 midpt= ICHAR(num(3))
 mdpt = 3550
 GO TO  410
 380 fl   = 8
 sec  = 0
 GO TO  410
 390 fl   = 9
 IF (nlp-4 == nlp4) nlp = nlp - 1
 GO TO  410
 400 fl   = 10
 GO TO  410
 405 fl   = 11
 GO TO  410
 
!     OPEN THE REQUESTED NASTRAN MANUAL FILE
 
 410 FILE = filex(fl)
 was1 = complf(0)
 was2 = was1
 OPEN (UNIT=lu,FILE=fname,ACCESS='SEQUENTIAL',FORM='FORMATTED',  &
     STATUS='OLD',ERR=412)
 GO TO 415
 412 WRITE  (out,413) fname
 413 FORMAT (//1X,a8,' file DOES NOT EXIST')
 IF (jc /= ku) GO TO 240
 GO TO 2210
 
 415 CONTINUE
 IF (check) GO TO 1700
 IF (PRINT) GO TO 2000
 sount = 0
!     SEARCH KEY WORD OR NUMBER AS APPROPRIATE, BASED ON MANUAL FLAG.
!     1=EXEC,2=CASE,3=BULK,4=PLOT,5=DMAP,6=SUBS,7=MSSG,8=DICT,
!     9=INTR,10=UMFL,11=RFMT
 420 SELECT CASE ( fl )
   CASE (    1)
     GO TO 650
   CASE (    2)
     GO TO 650
   CASE (    3)
     GO TO 650
   CASE (    4)
     GO TO 650
   CASE (    5)
     GO TO 650
   CASE (    6)
     GO TO 650
   CASE (    7)
     GO TO 430
   CASE (    8)
     GO TO 540
   CASE (    9)
     GO TO 610
   CASE (   10)
     GO TO 610
   CASE (   11)
     GO TO 650
 END SELECT
 
!     MESSAGE SEARCH IN FILE MSSG.TXT
 
 430 WRITE  (out,440) tctf
 440 FORMAT (/,' Enter MESSAGE NUMBER (up to 4 digits) or STOP,',a35)
 READ   (in,450) a4
 450 FORMAT (a4)
 IF (a11 == up1) GO TO 300
 IF (a4 == stp4) GO TO 2200
 IF (a4 == stp4l) GO TO 2200
 IF (a4 == exit4) GO TO 2200
 IF (a4 == exit4l) GO TO 2200
 IF (a4 == quit4) GO TO 2200
 IF (a4 == quit4l) GO TO 2200
 REWIND lu
 count = 0
 j = 0
 460 IF (a1(3) /= bnk1) GO TO 470
 j = j + 1
 IF (j > 2) GO TO 430
 a1(3) = a1(2)
 a1(2) = a1(1)
 a11   = bnk1
 GO TO 460
 470 IF (a1(4) == bnk1) GO TO 490
 IF (mdpt == 0 .OR. ICHAR(a11) < midpt) GO TO 490
 DO  i = 1,mdpt
   READ (lu,450,END=500)
 END DO
 count = mdpt
 490 READ (lu,450,END=500) key4
 count = count + 1
 IF (key4 == bnk4 .OR. key4 == new4 .OR. key4 /= a4) GO TO 490
 BACKSPACE lu
 kount = 2
 GO TO 1200
 500 WRITE  (out,510)
 510 FORMAT (' *** No Such MESSAGE NO. ***',/)
 GO TO 430
 
!     DICTIONARY SEARCH IN FILE DICT.TXT
 
 540 WRITE  (out,550) tctf
 550 FORMAT (/,' Enter DICTIONARY word or STOP,',a35)
 READ   (in,560) a6
 560 FORMAT (a6)
 IF (a11 == up1) GO TO 300
 DO  j = 1,6
   IF (ICHAR(a1(j)) >= la .AND. ICHAR(a1(j)) <= lz)  &
       a1(j) = CHAR(ICHAR(a1(j))+aa)
 END DO
 IF (a4 == stp4) GO TO 300
 IF (a4 == quit4) GO TO 300
 IF (a4 == stp4l) GO TO 300
 IF (a4 == quit4l) GO TO 300
 REWIND lu
 580 READ (lu,560,END=590) key6
 IF (key4 == bnk4 .OR. key4 == new4 .OR. key6 /= a6) GO TO 580
 BACKSPACE lu
 kount = 1
 GO TO 1200
 590 WRITE  (out,600)
 600 FORMAT (' *** No such term in NASTRAN Dictionary')
 GO TO 540
 
!     INTRODUCTION OR USER-MASTER-FILE SEARCH OF FILES INTR OR UMFL.TXT
 
!     SEARCH BY SECTION, SUBSECTION, AND PAGE ONLY
!     SECTION AND SUBSECTION MUST BE PRECEEDED BY A BLANK LINE
 
 610 WRITE  (out,620) tctf
 620 FORMAT (/,' Enter section(S), sub-section(B), page(P), ',  &
     'or stop(STOP),', /,a35)
 GO TO 670
 
 630 pass = pass + 1
 IF (pass >= 2) WRITE (out,640) fname
 640 FORMAT (17X,'No such WORD in ',a8,' file')
 REWIND lu
 count = 0
 IF (pass == 1) GO TO 700
 
!     GENERAL SEARCH FOR FILE TYPES = EXEC, CASE, BULK, PLOT, DMAP,
!     SUBS OR RFMT.TXT
!     KEY WORD or SECTION SEARCH
 
 650 WRITE  (out,660) tctf
 660 FORMAT (/,' Enter NASTRAN KEY WORD, STOP, next page(P),',  &
     ' next section(S),',/,' next sub-section(B), ',a35)
 670 READ   (in,680) key8
 680 FORMAT (a8)
 pass = 0
 last = count
 IF (k1(1) == up1) GO TO 230
!     NOT A ^ CHARACTER, CONVERT TO UPPER CASE AND PROCESS
 DO  iloop = 1,8
   IF (ICHAR(k1(iloop)) >= la .AND. ICHAR(k1(iloop)) <= lz)  &
       k1(iloop) = CHAR(ICHAR(k1(iloop))+aa)
 END DO
 IF (key4 == stp4) GO TO 1600
 IF (key4 == quit4) GO TO 1600
 IF (key4 == bnk4) GO TO 650
 IF (key4 == nxtp) GO TO 980
 IF (key4 == nxts) GO TO 750
 IF (key4 == nxtb) GO TO 730
 IF (key8 == dbgo) debug = .true.
 IF (key8 == dbgf) debug = .false.
 IF (key8 == dbgo .OR. key8 == dbgf) GO TO 650
 IF (fl == 9 .OR. fl == 10) GO TO 610
 700 jdx = 9
 IF (sec == 1) GO TO 900
 IF (k1(1) /= ks .AND. k1(1) /= kb) GO TO 900
 
 710 WRITE  (out,720) fname
 720 FORMAT (/,' *** Search by SECTION is not practical on this ',a8, ' file')
 GO TO 650
 
 730 key4 = nxtb
 IF (sount > 0) GO TO 755
 WRITE  (out,740)
 740 FORMAT (/,' *** SUBSECTION is requested without first request of',  &
     ' SECTION ***')
 GO TO 650
 750 key4  = nxts
 sount = 0
 IF (count <= 1) GO TO 1200
 755 a4   = stp4
 760 was4 = a4
 READ (lu,770,END=880) a12
 770 FORMAT (a12)
 count = count + 1
 IF (pass == 2 .AND. count == last) GO TO 1620
 IF (a4 == bnk4 .OR. a4 == new4) GO TO 760
 IF (was4 /= bnk4) GO TO 760
 i = ICHAR(a11)
 IF (a11 == bnk1) i = ICHAR(a1(2))
 IF (i < num1 .OR. i > num9) GO TO 760
 ndot = 0
 DO  i = 2,11
   IF (a1(i) /= dot) CYCLE
   IF (a1(i+1) /= bnk1 .AND. a1(i+1) /= dot) ndot = ndot + 1
 END DO
 IF (ndot-1 < 0) THEN
   GO TO   760
 ELSE IF (ndot-1 == 0) THEN
   GO TO   790
 ELSE
   GO TO   810
 END IF
 790 IF (key4 == nxtb) GO TO 820
 sount = count - 1
 800 count = count - 1
 BACKSPACE lu
 GO TO 1200
 810 IF (key4 == nxts) GO TO 760
 GO TO 800
 820 WRITE  (out,830) tctf
 830 FORMAT (/,' *** End of SECTION ***', /,' return to Key(K), ',  &
     'return to beginning of section(R), next section(N)', /,' stop(STOP),',a35)
 READ (in,170) ic
 IF (ICHAR(ic) >= la .AND. ICHAR(ic) <= lz) ic = CHAR(ICHAR(ic)+aa)
 IF (ic == up1) GO TO 230
 IF (ic == ks ) GO TO 2200
 IF (ic == kk ) GO TO 650
 IF (ic == kn ) GO TO 860
 IF (ic /= kr ) GO TO 820
 IF (sount <= 1) GO TO 870
 j = count - sount + 2
 840 DO  i = 1,j
   BACKSPACE lu
   count = count - 1
 END DO
 IF (count < 0) count = 0
 key4 = nxts
 GO TO 1200
 860 j = 1
 GO TO 840
 870 REWIND lu
 key4 = nxts
 GO TO 1200
 880 WRITE  (out,890)
 890 FORMAT (' *** End of File ***')
 REWIND lu
 sount = 0
 GO TO 650
 
!     KEY WORD SEARCH - FIRST SEARCH HEADING THEN KEY WORD
 
 900 jdx = jdx - 1
 IF (k1(jdx) == bnk1) GO TO 900
 IF (jdx <= 0) GO TO 630
 IF (fl  /= 1) GO TO 910
 
!     SOME KEY WORDS IN EXECUTIVE CONTROL SECTION MAY BE ABBREVIATED.
!     4 BYTES ARE USED FOR ALL EXECUTIVE CONTROL KEY WORDS
 
 IF (jdx > 4) jdx = 4
 IF (key4 == appr) key4 = app
 IF (key4 == solu) key4 = sol
 
!     IF KEY IS LESS THAN 4 LETTERS, ADD A BLANK AT THE END SO THAT
!     'SOF' IS NOT 'SOFIN', 'SOFOUT', etc.
 
 910 IF (jdx >= 4) GO TO 920
 jdx = jdx + 1
 k1(jdx) = bnk1
 
!     T1  = 2**8,  T2 = 2**16,  T3 = 2**24
!     IS0 = FIRST CHARACTER OF THE 8-BYTE KEY WORD IN NUMERIC VALUE
!     IS1 = FIRST  HALF OF THE 8-BYTE KEY WORD IN NUMERIC VALUE
!     IS2 = SECOND HALF OF THE 8-BYTE KEY WORD IN NUMERIC VALUE
 
!     THAT IS, WE WILL USE NUMERIC VALUE FOR KEY WORD SEARCH
 
 920 is0 = ICHAR(k1(1))
 is1 = is0*t1 + ICHAR(k1(2))*t2 + ICHAR(k1(3))*t3 + ICHAR(k1(4))
 is2 = ICHAR(k1(5))*t1 + ICHAR(k1(6))*t2 + ICHAR(k1(7))*t3 + ICHAR(k1(8))
 
!     COMPARE PRESENT KEY WORD AND PREVIOUS KEY AND DETERMINE WE NEED
!     TO REWIND FILE OR NOT
 
!     IF TEXT FILE IS NOT PRESORTED, WE NEED TO REWIND FILE ON EACH NEW
!     KEY WORD.  (USER'S MANUAL IS SORTED)
 
 IF (is1-was1 < 0) THEN
   GO TO   940
 ELSE IF (is1-was1 == 0) THEN
   GO TO   930
 ELSE
   GO TO   980
 END IF
 930 IF (is2-was2 < 0) THEN
   GO TO   940
 ELSE IF (is2-was2 == 0) THEN
   GO TO   950
 ELSE
   GO TO   980
 END IF
 940 REWIND lu
 count = 0
 GO TO 980
 950 WRITE (out,960)
 960 FORMAT (' Same KEY WORD as before. Continue? (Y,N) ')
 READ (in,170) ic
 IF (ICHAR(ic) >= la .AND. ICHAR(ic) <= lz) ic = CHAR(ICHAR(ic)+aa)
 IF (ic == no) GO TO 650
 REWIND lu
 count = 0
 
!     IF KEY WORD IS BEYOND MID-POINT, SKIP HALF OF THE RECORDS IN FILE
 
 IF (mdpt == 0 .OR. ICHAR(k1(1)) < midpt) GO TO 980
 l = mdpt - count + 1
 DO  j = 1,l
   READ (lu,170)
 END DO
 count = count + l
 
!     IVFF IS PAGE MARK.  PAG6 IS '=PAGE='
!     LOOK FOR PAGE MARK OR '=PA' FIRST
 
 980 READ (lu,170,END=630) jx1,ic,jc
 count = count + 1
 IF (pass == 2 .AND. count == last) GO TO 1620
 ivjx1 = ICHAR(jx1)
 IF (jx1 /= lb1 .AND. ivjx1 /= ivff .AND. a3 /= pag3) GO TO 980
 IF (debug) WRITE (out,990) jx1,ic,jc
 990 FORMAT (40X,'@980 Just read- ',8A1)
 IF (key4 == nxtp .AND. a3 == pag3) GO TO 1240
 j = 0
 1000 READ (lu,1250,END=1450) a80
 IF (debug) WRITE (out,1005) (a1(i),i=1,8)
 1005 FORMAT (36X,'@1005 Just read- ',8A1)
 count = count + 1
 IF (pass == 2 .AND. count == last) GO TO 1620
 j = j + 1
 IF (j >= 7) GO TO 980
 IF (a4 == bnk4 .OR. a4 == new4) GO TO 1000
 IF (ICHAR(a11) == ivff) GO TO 1000
 
!     KEY WORD HEADING SEARCH
 
!     ******************************************************
!     *   HEADER WORDS WERE REMOVED IN 1993 USER'S MANUAL  *
 j = 0
 IF (j == 0) GO TO 1120
!     ******************************************************
 
 SELECT CASE ( fl )
   CASE (    1)
     GO TO 1010
   CASE (    2)
     GO TO 1010
   CASE (    3)
     GO TO 1010
   CASE (    4)
     GO TO 1010
   CASE (    5)
     GO TO 1030
   CASE (    6)
     GO TO 1010
   CASE (    7)
     GO TO 1100
   CASE (    8)
     GO TO 1100
   CASE (    9)
     GO TO 1100
   CASE (   10)
     GO TO 1100
   CASE (   11)
     GO TO 1050
 END SELECT
 1010 IF (debug) WRITE (out,1020) a12,hd12(fl)
 1020 FORMAT (50X,a12,'==> ',a12)
 IF (a12 == hd12(fl)) GO TO 1050
 IF (fl == 4 .AND. a12 == hd12(5)) GO TO 1040
 GO TO 1000
 1030 IF (a6 /= hd6) GO TO 1000
 GO TO 1050
 
 1040 j = j4
 GO TO 1060
 1050 j = base(fl)
 1060 IF (fl /= 6) GO TO 1070
 IF (a1(14) ==  km) j = j5
 IF (a1(14) ==  ko) j = j6
 1070 IF (a1(j) /= bnk1) j = j - 1
 IF (debug) WRITE (out,1080) j
 1080 FORMAT (45X,'@1080  BASE J =',i3)
 1090 IF (a1(j+1) /= bnk1) GO TO 1120
 j = j + 1
 GO TO 1090
 
 1100 WRITE  (out,1110) fl
 1110 FORMAT (/,' *** SHOULD NOT BE HERE.  FL =',i3)
 GO TO 240
 
!     KEY WORD SEARCH
 
 1120 IF (debug) WRITE (out,1130) (a1(j+i),i=1,jdx),lt1,lt1, (k1(i),i=1,jdx)
 1130 FORMAT (30X,'@1130 - ',18A1)
 DO  i = 1,jdx
   IF (a1(j+i) /= k1(i)) GO TO 980
 END DO
 
!     KEY WORD FOUND ON FILE
 
 was1  = is1
 was2  = is2
 kount = 6
 IF (first) kount = 8
 WRITE  (out,1150)
 1150 FORMAT (//)
 IF (     ety) WRITE (out,1290) a80
 IF (.NOT.ety) WRITE (out,1300) a79
 
!     RECORD FOUND.  READ AND PRINT ON SCREEN
!     ALLOW UP TO 4 BLANK LINES PRINTED ON SCREEN
 
 1200 jount = count - 4
 bline = 0
 IF (numcov == 0) GO TO 1240
 DO  j = 1,80
   DO  i = 1,numcov
     IF (ICHAR(a1(j)) == ival(i)) a1(j) = CHAR(nval(i,1))
   END DO
 END DO
 GO TO 1240
 
 1220 WRITE  (out,1230)
 1230 FORMAT (///)
 kount = kount + 3
 
 1240 READ (lu,1250,END=1450) a80
 1250 FORMAT (a80)
 count = count + 1
 IF (ICHAR(a11) == ivff .OR. a6 == pag6) GO TO 1220
 IF (a48 /= b48) GO TO 1260
 IF (bline > 4) GO TO 1240
 bline = bline + 1
 WRITE (out,170) a11
 GO TO 1310
 1260 bline = 0
 IF (numcov == 0) GO TO 1280
 DO  j = 1,80
   DO  i = 1,numcov
     IF (ICHAR(a1(j)) == ival(i)) a1(j) = CHAR(nval(i,1))
   END DO
 END DO
 1280 IF (a11  ==  lb1) GO TO 1470
 IF (     ety) WRITE (out,1290) a80
 IF (.NOT.ety) WRITE (out,1300) a79
 1290 FORMAT (1X,a80)
 1300 FORMAT (1X,a79)
 1310 kount = kount + 1
 IF (MOD(kount,nlp) /= 0) GO TO 1240
 SELECT CASE ( fl )
   CASE (    1)
     GO TO 1320
   CASE (    2)
     GO TO 1320
   CASE (    3)
     GO TO 1320
   CASE (    4)
     GO TO 1320
   CASE (    5)
     GO TO 1320
   CASE (    6)
     GO TO 1320
   CASE (    7)
     GO TO 420
   CASE (    8)
     GO TO 420
   CASE (    9)
     GO TO 1320
   CASE (   10)
     GO TO 1320
   CASE (   11)
     GO TO 1320
 END SELECT
 1320 IF (.NOT.first) GO TO 1350
 first = .false.
 WRITE (out,1332) nlp4
 1332 FORMAT(' (Y,N,STOP,1,2,...,',i2,',-n,P,S,B,^,PRINT,HELP or <CR>)')
 GO TO 1350
 1330 WRITE (out,1340) nlp4,nlp,nlp4
 1340 FORMAT (' (Y,N,STOP,1,2,...,',i2,',-n,P,S,B,^,PRINT,HELP or <CR>)'  &
     /11X,'Y or <CR> = yes more', /11X,'N         = no more on this item',  &
     /11X,'STOP      = terminate NASTHELP',  &
     /11X,'1,2,...,n = keep bottom n lines on next page. (',i2, ' max)',  &
     /11X,'-n        = back up n+',i2,' lines',  &
     /11X,'P,S,B     = go to next page, next section, or next', ' sub-section',  &
     /11X,'^         = terminate current text file',  &
     /11X,'PRINT     = print text, up to last line on screen',  &
     /11X,'HELP      = echo options of MORE', /,' ...more? ')
 GO TO 1370
 1350 WRITE  (out,1360)
 1360 FORMAT (' ...more? ')
 1370 READ (in,170) ic,jc,lc
 IF (ic == up1) GO TO 230
 IF (ICHAR(lc) >= la .AND. ICHAR(lc) <= lz) lc = CHAR(ICHAR(lc)+aa)
 IF (ICHAR(ic) >= la .AND. ICHAR(ic) <= lz) ic = CHAR(ICHAR(ic)+aa)
 IF (ICHAR(jc) >= la .AND. ICHAR(jc) <= lz) jc = CHAR(ICHAR(jc)+aa)
 IF (ic == kh) GO TO 1330
 IF (ic == no) GO TO 420
!     Check for STop, EXit ot QUit
 IF (ic == ks .AND. jc == kt) GO TO 1600
 IF (ic == ke .AND. jc == kx) GO TO 1600
 IF (ic == kq .AND. jc == ku) GO TO 1600
!     Check for PRint request.
 IF (ic == kp .AND. jc == kr) GO TO 2000
 IF (sec == 0 .AND. (ic == ks .OR. ic == kb)) GO TO 710
 IF (ic == ks) GO TO 750
 IF (ic == kb) GO TO 1380
 IF (ic /= kp) GO TO 1390
 key4 = nxtp
 GO TO 980
 1380 IF (sount > 0) GO TO 730
 WRITE (out,740)
 GO TO 1350
 1390 kount = 0
 IF (ic == bnk1 .OR. ic == yes) GO TO 1240
 i =  0
 j = -1
 l = -1
 DO  k = 1,10
   IF (ic == num(k)) i = MOD(k,10)
   IF (jc == num(k)) j = MOD(k,10)
   IF (lc == num(k)) l = MOD(k,10)
 END DO
 IF (j+l == -2) ijl = i
 IF (l == -1 .AND. j /= -1) ijl = i*10  + j
 IF (l /= -1 .AND. j /= -1) ijl = i*100 + j*10 + l
 IF (ic == mns) GO TO 1410
 ijl   = MIN0(nlp4,ijl)
 kount = ijl + 1
 GO TO 1240
 1410 ijl   = ijl + nlp
 DO  l = 1,ijl
   BACKSPACE lu
 END DO
 count = count - ijl
 kount = 0
 IF (count < 0) count = 0
 GO TO 1240
 
 1450 IF (lc == ku) GO TO 1610
 WRITE  (out,1460)
 1460 FORMAT (/,' ...EOF. <CR> to continue')
 GO TO 1500
 
 1470 WRITE  (out,1480)
 1480 FORMAT (/,' ...End of Description')
 1500 IF (kount < nlp4) GO TO 650
 READ (in,170) ic
 IF (ICHAR(ic) >= la .AND. ICHAR(ic) <= lz) ic = CHAR(ICHAR(ic)+aa)
 IF (ic /= ks) GO TO 650
 
 1600 IF (mq == 0) GO TO 2200
 1610 CLOSE (lu)
 kount = 0
 SELECT CASE ( fl )
   CASE (    1)
     GO TO 320
   CASE (    2)
     GO TO 330
   CASE (    3)
     GO TO 1620
 END SELECT
 1620 WRITE  (out,1630)
 1630 FORMAT (/,'*** No such KEY WORD in EXEC, CASE and BULK.TXT files')
 GO TO 240
 
!     CHECK INPUT TEXT FORMATS FOR
!     EXEC, CASE, BULK, PLOT, DMAP SUBS and RFMT.TXT FILES
 
 1700 IF (fl <= 6 .OR. fl == 11) GO TO 1720
 WRITE  (out,1710) filex(fl)
 1710 FORMAT (//,' *** Data CHECK OPTION not valid for ',a4,'.TXT file')
 GO TO 240
 1720 READ (lu,1730,END=2200) a48
 1730 FORMAT (a48)
 IF (a12 == hd12(fl)) GO TO 1760
 IF (fl == 4 .AND. a12 == hd12(5)) GO TO 1760
 IF (a11 /= lb1) GO TO 1720
 1740 READ (lu,1730,END=2200) a48
 IF (a4 == bnk4 .OR. a4 == new4) GO TO 1740
 IF (a12 == hd12(fl)) GO TO 1790
 IF (fl == 4 .AND. a12 == hd12(5)) GO TO 1780
 WRITE  (out,1750) a48
 1750 FORMAT (1X,a48,' <== HEADER 12 CHAR. ERROR')
 GO TO 1720
 1760 WRITE  (out,1770) a48
 1770 FORMAT (1X,a48,' <== NO PRECEEDING # SYMBOL')
 GO TO 1720
 
 1780 j = j4
 GO TO 1800
 1790 j = base(fl)
 1800 IF (fl /= 6) GO TO 1810
 IF (a1(14) == km) j = j5
 IF (a1(14) == ko) j = j6
 1810 IF (a1(j  ) /= bnk1) j = j - 1
 1820 IF (a1(j+1) /= bnk1) GO TO 1830
 j = j + 1
 GO TO 1820
 1830 is1 = ICHAR(a1(j+1))*t1 + ICHAR(a1(j+2))*t2 + ICHAR(a1(j+3))*t3 &
          + ICHAR(a1(j+4))
 is2 = ICHAR(a1(j+5))*t1 + ICHAR(a1(j+6))*t2 + ICHAR(a1(j+7))*t3  &
     + ICHAR(a1(j+8))
 IF (is1-was1 < 0) THEN
   GO TO  1850
 ELSE IF (is1-was1 == 0) THEN
   GO TO  1840
 ELSE
   GO TO  1890
 END IF
 1840 IF (is2-was2 < 0) THEN
   GO TO  1850
 ELSE IF (is2-was2 == 0) THEN
   GO TO  1870
 ELSE
   GO TO  1890
 END IF
 1850 WRITE  (out,1860) a48
 1860 FORMAT (1X,a48,' <== OUT OF ORDER')
 GO TO  1910
 1870 WRITE  (out,1880) a48
 1880 FORMAT (1X,a48,' <== DUPLICATE')
 GO TO  1720
 1890 WRITE  (out,1900) a48
 1900 FORMAT (1X,a48)
 1910 was1 = is1
 was2 = is2
 GO TO 1720
 
!     PROCESS REQUEST TO PRINT MANUAL DATA FOUND
 
 2000 IF (popen) GO TO 2040
 WRITE  (out,2010)
 2010 FORMAT (/,1H-,' Enter OUTPUT FILE Name (assume default dir): ')
 READ   (in,770) a12
 IF (a4 /= bnk4 .AND. a4 /= stp4 .AND. a4 /= stp4l) GO TO 2030
 IF (a4 /= exit4 .AND. a4 /= exit4l) GO TO 2030
 IF (a4 /= quit4 .AND. a4 /= quit4l) GO TO 2030
 WRITE  (out,2020)
 2020 FORMAT (/,' Output PRINT Aborted')
 GO TO 2200
 2030 ou12 = a12
 IF (popen) GO TO 2040
 OPEN (UNIT=ou,FILE=ou12,STATUS='NEW',ACCESS='SEQUENTIAL',ERR=2140,  &
     FORM='FORMATTED')
 popen = .true.
 GO TO 2060
 2040 WRITE  (ou,2050)
 2050 FORMAT (1H1)
 2060 j = count - jount
 IF (j <= 0) GO TO 2120
 DO  i = 1,j
   BACKSPACE lu
 END DO
 DO  i = 1,j
   READ (lu,1250,END=2120) a80
   IF (a11 ==  lb1) CYCLE
   IF (a6  == pag6) GO TO 2090
   WRITE  (ou,2080) a80
   2080 FORMAT (1X,a80)
   CYCLE
   2090 WRITE  (ou,2050)
 END DO
 
 2120 WRITE  (out,2130) j
 2130 FORMAT (/,i9,' lines printed')
 GO TO 420
 
 2140 WRITE  (out,413) ou12
 GO TO 2000
 
!     END OF JOB PROCESSING
 
 2200 CLOSE  (lu)
 kount = 0
 2210 WRITE  (out,2220)
 2220 FORMAT (//,'  *** NASTHELP is done.  Have a good run! ***',// )
 IF (.NOT.popen) GO TO 2240
 CLOSE  (ou)
 WRITE  (out,2230) ou12
 2230 FORMAT (/,'  *** Don''T FORGET THE PRINT FILE IN ',a42)
 popen = .false.
 2240 CONTINUE
 
END PROGRAM nasthelp
