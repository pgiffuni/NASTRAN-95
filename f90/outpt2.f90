SUBROUTINE outpt2
     
!     COPY DATA BLOCK(S) ONTO FORTRAN UNIT.
 
!     CALL TO THIS MODULE IS
 
!     OUTPUT2    IN1,IN2,IN3,IN4,IN5/ /V,N,P1/V,N,P2/V,N,P3/
!                                      V,N,P4/V,N,P5/V,N,P6 $
 
!             P1 = 0, NO ACTION TAKEN BEFORE WRITE
!                     (DEFAULT P1=0)
!                =+N, SKIP FORWARD N DATA BLOCKS BEFORE WRITE
!                =-1, BEFORE WRITE, FORTRAN TAPE IS REWOUND AND A
!                     HEADER RECORD (RECORD NUMBER 0) ADDED TO TAPE
!                =-3, THE NAMES OF ALL DATA BLOCKS ON FORTRAN TAPE
!                     ARE PRINTED AND WRITE OCCURS AT THE END OF TAPE
!                =-9, WRITE A NULL FILE, ENDFILE AND REWIND FORTRAN
!                     TAPE.
 
!             P2 =    THE FORTRAN UNIT NO. ON WHICH THE DATA BLOCKS WILL
!                     BE WRITTEN. (DEFAULT P2=11).
 
!             P3 =    TAPE ID CODE FOR FORTRAN TAPE, AN ALPHANUMERIC
!                     VARIABLE WHOSE VALUE WILL BE WRITTEN ON A FORTRAN
!                     TAPE.
!                     THE WRITING OF THIS ITEM IS DEPENDENT ON THE
!                     VALUE OF P1 AS FOLLOWS.
!                          *P1*             *TAPE ID WRITTEN*
!                           +N                     NO
!                            0                     NO
!                           -1                    YES
!                           -3                     NO (WARNING CHECK)
!                    (DEFAULT P3 = XXXXXXXX).
 
!             P4 = 0, FORTRAN WRITTEN RECORD SIZE IS UNLIMITTED
!                     (DEFAULT FOR ALL MACHINES, EXECPT IBM)
!                =-N, MAXIMUM FORTRAN WRITTEN RECORD SIZE IS N TIMES
!                     THE SYSTEM BUFFER SIZE, N*BUFFSIZE
!                = N, MAXIMUM FORTRAN WRITTEN RECORD SIZE IS N WORDS.
!                -    IN ALL CASES, THE MAXIMUM FORTRAN WRITTEN RECORD
!                     SIZE SHOLD BE .GE. BUFFSIZE, AND .LE. AVAILABLE
!                     CORE
!                IBM, IF P4=0, AND SINCE IBM CAN NOT HANDLE UNLIMITED
!                     RECORD SIZE, RECORD SIZE P4 OF 1024 WORDS IS USED
 
!             P5 = 0  FOR NON-SPARSE, AND NON-ZERO FOR SPARSE MATRIX
!                     OUTPUT
!                = 0, KEY-WORD RECORD CONTAINS EFFECTIVELY ONE SINGLE
!                     WORD OF DATA (THIS IS THE ORIGINAL COSMIC/OUTPT2)
!                = NOT 0, KEY-WORD RECORD CONTAINS 2 WORDS, THUS ALLOW
!                     SPARSE MATRIX TO BE COPIED OUT.
!                     FIRST KEY WORD:
!                        >0, DEFINES THE LENGTH OF NEXT DATA RECORD
!                        =0, END-OF-FILE
!                        <0, END-OF-RECORD WITH ANOTHER RECORD TO FOLLOW
!                     SECOND KEY WORD:
!                        =0, TABLE DATA, OR P5 SPARSE MATRIX OPTION NOT
!                            REQUESTED
!                        >0, ROW-BASE FOR NEXT RECORD. FOR EXAMPLE:
!                            KEYS = 10,200 INDICATE NEXT DATA RECORD IS
!                            FOR ROW(200+1) THRU ROW(200+10)
!                            i.e. (ROW(KEY2+J),J=1,KEY1)
 
!             P6 = BLANK (DEFAULT)
!                = *MSC*,    OUTPUT2 WILL ISSUE RECORDS IN MSC/OUTPUT2
!                            FORMAT WHICH IS SLIGHTLY DIFFERENT FROM
!                            COSMIC/OUTPUT2.
!                            (P5 OPTION IS NOT AVAILABLE)
 
!     NOTES ABOUT P5
!             (1) P5 IS IGNORED IN TABLE DATA
!             (2) POSSIBLY, NON-ZERO ROW ELEMENT MAY START AT 2ND HALF
!                 OF A COMPLEX WORD
!             (3) UP TO 3 ZEROS MAY BE IMBEDDED IN NON-ZERO STRING
!             (4) THE CHOICE OF 2 KEY WORDS IN ONE KEY RECORD OVER 2 KEY
!                 WORDS IN TWO RECORDS (AS IN MSC/NASTRAN), IS NOT TO
!                 MAKE THE ORIGINAL COSMIC OUTPT2/INPTT2 OBSOLETE.
!                 (i.e. WE DON'T FOLLOW OTHER PEOPLE BLINDLY SO TO MAKE
!                 OURSELVES OBSOLETE)
!             (5) ALTHOUGH OUTPT2 ALWAYS WRITES 2 KEY WORDS OUT IN A
!                 RECORD. ONE MAY CHOOSE TO READ BACK ONE OR BOTH KEYS.
 
!     REVISED  11/90 BY G.CHAN/UNISYS TO INCLUDE P4 AND P5 PARAMETERS
!     LAST REVISED  2/93 BY G.CHAN    TO INCLUDE P6 PARAMETER
 
 IMPLICIT INTEGER (a-z)
 LOGICAL :: sparse,dp
 CHARACTER (LEN=6) :: mt,matrix,table
 DIMENSION       dx(3),trl(8),NAME(2),subnam(2),inp(3),namex(2),  &
     idhdr(7),idhdrx(7),p3x(2),tapcod(2)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
!WKBNB
 CHARACTER (LEN=80) :: dsnames
 COMMON /dsname/ dsnames(80)
!WKBNE
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /BLANK / p1,p2,p3(2),p4,p5,p6(2)  &
     /system/ buffsz,nout,dum6(6),nlpp,dum2(2),line,dum(2),d(3) /zzzzzz/ x(1)  &
     /machin/ mach
 COMMON /unpakx/ itype,irow,nrow,incr
 DATA    subnam/ 4HOUTP, 4HT2  /, matrix,table /'MATRIX',' TABLE'/
 DATA    inp   / 1HT, 1H1, 1H2 /, msc   / 4HMSC   /
 DATA    zero  , mone,mtwo,mtre,mnin    / 0,-1,-2,-3,-9 /
 DATA    idhdr / 4HNAST,4HRAN ,4HFORT,4H tap,4HE id,4H cod,4HE - /
!WKBI
 DATA    ifirst/0/
 
!     CHECK P2 AND P4 PARAMETERS
 
!WKBI 3/95 SPR94016
 lcor   = korsz(x) - buffsz
 IF (p2 >= 11 .AND. p2 <= 21 ) GO TO 20
 j = 11
 WRITE  (nout,10) uwm,p2,j,inp(i)
 10 FORMAT (a25,' FROM OUTPUT2 MODULE. UNACCEPTABLE FORTRAN UNIT',i3,  &
     ' WAS CHANGED TO',i3,' (INP',a1,1H))
 p2 = j
 20 IF (p4 < 0.0) THEN
   GO TO    25
 ELSE IF (p4 == 0.0) THEN
   GO TO    30
 ELSE
   GO TO    35
 END IF
 25 lrec = -p4*buffsz
 GO TO 40
 30 lrec = lcor
 IF (mach  ==   2) lrec = 1024
 IF (p6(1) == msc) lrec = 2*buffsz
 GO TO 40
 35 lrec = p4
 40 IF (lrec >   lcor) lrec = lcor
 IF (lrec < buffsz) lrec = buffsz
 IF (p4 /= 0) WRITE (nout,50) uim,lrec
 50 FORMAT (a29,' 4116, MAXIMUM FORTRAN RECORD SIZE USED IN OUTPUT2 ',  &
     'WAS',i8,' WORDS')
 p4 = lrec
!WKBNB
 IF ( ifirst /= 0 ) GO TO 51
 CLOSE ( UNIT=p2 )
 OPEN ( UNIT=p2, FILE=dsnames(p2), FORM='UNFORMATTED', STATUS='UNKNOWN' )
 ifirst = 1
 51    CONTINUE
!WKBNE
 IF (p6(1) == msc) CALL outmsc (*1000,*420)
 
 sparse = .false.
 IF (p5 /= 0) sparse = .true.
 endfil = 0
 endrec = 0
!WKBD 3/95 SPR94016      LCOR   = KORSZ(X) - BUFFSZ
 icrq   =-lcor
 IF (lcor <= 0) GO TO 890
 inbuf  = lcor + 1
 tapcod(1) = p3(1)
 tapcod(2) = p3(2)
 out = p2
 IF (p1 == mnin) GO TO 410
 IF (p1 < mtre .OR. p1 == mtwo) GO TO 810
 
 IF (p1 == mtre) GO TO 500
 IF (p1 <= zero) GO TO 80
 
 i = 1
 60 READ (out) key
 keyx = 2
 IF (key /= keyx) GO TO 900
 READ (out) namex
 READ (out) key
 IF (key >= 0) GO TO 920
 ASSIGN 70 TO ret
 nskip = 1
 GO TO 700
 70 i = i + 1
 IF (i <= p1) GO TO 60
 
 80 IF (p1 /= mone) GO TO 90
 
!     REWIND OUTPUT TAPE. (P1 = -1)
 
 REWIND out
 key = 3
 WRITE (out) key,zero
 WRITE (out) d
 key = 7
 WRITE (out) key,zero
 WRITE (out) idhdr
 key = 2
 WRITE (out) key,zero
 WRITE (out) p3
 endrec = endrec - 1
 WRITE (out) endrec,zero
 WRITE (out) endfil,zero
 endrec = 0
 
 90 DO  i = 1,5
   INPUT  = 100 + i
   trl(1) = INPUT
   CALL rdtrl (trl)
   IF (trl(1) <= 0) CYCLE
   CALL fname (INPUT,NAME)
   
!     OPEN INPUT DATA BLOCK TO READ WITH REWIND.
   
   CALL OPEN (*800,INPUT,x(inbuf),0)
   CALL skprec (INPUT,1)
   trl(8) = 1
   CALL rectyp (INPUT,irec1)
   IF (irec1 /= 0) GO TO 100
   trl(8) = 0
   CALL READ (*100,*100,INPUT,x(1),1,1,nf)
   CALL rectyp (INPUT,irec2)
   IF (irec2 == 0) GO TO 100
   trl(8) = 2
   100 CALL REWIND (INPUT)
   key = 2
   WRITE (out) key,zero
   WRITE (out) NAME
   endrec = endrec - 1
   WRITE (out) endrec,zero
   key = 8
   WRITE (out) key,zero
   WRITE (out) trl
   endrec = endrec - 1
   WRITE (out) endrec,zero
   INDEX = 0
   
!     COPY CONTENTS OF INPUT DATA BLOCK ONTO FILE.
!     (OR THE HEADER RECORD OF A MATRIX DATA BLOCK)
   
!     COMMENTS FROM G.CHAN/UNISYS  2/93
!     THE WRITES IN LOOP 110 AND 120 SEEM DATA TYPE (S.P. OR D.P.)
!     INCENSITIVE. THE D.P. DATA IN KELM, MELM AND BELM TABLES SHOULD
!     WORK OK.
   
   110 CALL READ (*310,*120,INPUT,x(1),lrec,0,nf)
   WRITE (out) lrec,zero
   WRITE (out) (x(l),l=1,lrec)
   GO TO 110
   
   120 WRITE (out) nf,zero
   WRITE (out) (x(l),l=1,nf)
   endrec = endrec - 1
   WRITE (out) endrec,zero
   IF (trl(8) == 0) GO TO 110
   IF (trl(8) == 1) GO TO 130
   IF (INDEX  > 0) GO TO 130
   INDEX = 1
   GO TO 110
   
!     COPY STRING FORMATTED MATRIX
   
   130 IF (trl(8) == 2 .AND. INDEX == 2) GO TO 140
   INDEX = 2
   nwds  = trl(5)
   dp    = .false.
   IF (nwds == 2 .OR. nwds == 4) dp = .true.
   dsp   = 1
   IF (dp) dsp = 2
   IF (nwds == 3) nwds = 2
!         NWDS=1,SP  -  =2,DP,CS  -  =4,CDP
   
   incr = 1
   nwds = trl(3)*nwds
   
!     CHECK FOR NULL MATRIX
   
   IF (trl(2) == 0 .OR. trl(3) == 0) GO TO 310
   
!     NWDS HAS NUMBER WORDS NEEDED PER COLUMN
   
   icrq = nwds - lcor
   IF (nwds > lcor) GO TO 890
   itype = trl(5)
   irow  = 1
   nrow  = trl(3)
   ncol  = trl(2)
   IF (trl(8) == 2) ncol = 1
   140 DO  l = 1,ncol
     CALL unpack (*180,INPUT,x)
     IF (sparse) GO TO 200
     150 DO  kb = 1,nwds,lrec
       ke = kb + lrec - 1
       IF (ke > nwds) ke = nwds
       kbe = ke - kb + 1
       WRITE (out) kbe,zero
       WRITE (out) (x(k),k=kb,ke)
     END DO
     
     170 endrec = endrec - 1
     WRITE (out) endrec,zero
     CYCLE
     180 IF (sparse) GO TO 170
     DO  k = 1,nwds
       x(k) = 0
     END DO
     GO TO 150
     
!     SPARSE MASTRIX OUT
     
     200 j12 = -1
     DO  j = 1,nwds,dsp
       IF (j12  >=  +1) GO TO 220
       IF (x(j) /= 0.0) GO TO 210
       IF (dp) IF (x(j+1)) 210,260,210
       CYCLE
       210 j12 = +1
       k2  = j - 1
       CYCLE
       220 IF (x(j) /= 0.0) CYCLE
       IF (dp) IF (x(j+1)) 260,230,260
       230 IF (j12 == -1) CALL mesage (-37,0,subnam)
       j12 = j12 + 1
       
!     ALLOW UP TO 3 IMBEDDED ZEROS
       
       IF (j12 <= 3) CYCLE
       IF (x(j-1) /= 0.0 .OR. x(j-2) /= 0.0) CYCLE
       j12 = -1
       k1  = j - k2
       IF (k1 > lrec) GO TO 240
       WRITE (out) k1,k2
       WRITE (out) (x(k2+k),k=1,k1)
       CYCLE
       240 ke = j
       kb = k2 + 1
       DO  kk = kb,ke,lrec
         k2 = kk - 1
         k1 = k2 + lrec
         IF (k1 > ke) k1 = ke
         WRITE (out) k1,k2
         WRITE (out) (x(k2+k),k=1,k1)
       END DO
       260 CONTINUE
     END DO
     
     IF (j12 == -1) GO TO 290
     j12 = -1
     k1  = nwds - k2
     IF (k1 >= lrec) GO TO 270
     WRITE (out) k1,k2
     WRITE (out) (x(k2+k),k=1,k1)
     GO TO 290
     270 ke = nwds
     kb = k2 + 1
     DO  kk = kb,ke,lrec
       k2 = kk - 1
       k1 = k2 + lrec
       IF (k1 > ke) k1 = ke
       WRITE (out) k1,k2
       WRITE (out) (x(k2+k),k=1,k1)
     END DO
     
     290 endrec = endrec - 1
     WRITE (out) endrec,zero
     
   END DO
   
   IF (trl(8) == 2) GO TO 110
   
!     CLOSE INPUT DATA BLOCK WITH REWIND
   
   310 CALL CLOSE (INPUT,1)
   
   WRITE (out) endfil,zero
   endrec = 0
   CALL page2 (-4)
   mt = matrix
   IF (trl(8) == 0) mt = table
   WRITE  (nout,320) uim,mt,NAME,out,(trl(ii),ii=2,7)
   320 FORMAT (a29,' 4114, ',a6,' DATA BLOCK ',2A4,  &
       ' WRITTEN ON FORTRAN UNIT',i4, /5X,'TRAILR =',5I6,i9)
   IF (sparse .AND. trl(8) /= 0) WRITE (nout,330)
   330 FORMAT (1H+,55X,'(SPARSE MATRIX)')
   
 END DO
 GO TO 1000
 
!     FINAL CALL TO OUTPUT2. (P1 = -9)
 
 410 WRITE (out) endfil,zero
 420 endrec = 0
 ENDFILE out
 REWIND out
 WRITE  (nout,430) uim
 430 FORMAT (a29,'. OUTPUT2 MODULE WROTE AN E-O-F RECORD, A SYSTEM ',  &
     'E-O-F MARK, AND REWOUND THE OUTPUT TAPE. (P1=-9)')
 GO TO 1000
 
!     OBTAIN LIST OF DATA BLOCKS ON FORTRAN TAPE.  (P1 = -3)
 
 500 REWIND out
 READ (out) key
 keyx = 3
 IF (key /= keyx) GO TO 900
 READ (out) dx
 READ (out) key
 keyx = 7
 IF (key /= keyx) GO TO 900
 READ (out) idhdrx
 DO  kf = 1,7
   IF (idhdrx(kf) /= idhdr(kf)) GO TO 830
 END DO
 READ (out) key
 keyx = 2
 IF (key /= keyx) GO TO 900
 READ (out) p3x
 IF (p3x(1) /= p3(1) .OR. p3x(2) /= p3(2)) GO TO 850
 520 ASSIGN 530 TO ret
 nskip = 1
 GO TO 700
 530 kf = 0
 540 CALL page1
 line = line + 8
 WRITE  (nout,550) out
 550 FORMAT (1H0, 50X, 30HFILE contents on fortran UNIT , i2,  &
     /51X, 32(1H-), ///54X, 4HFILE, 18X, 4HNAME/1H0)
 560 READ (out) key
 IF (key < 0) THEN
   GO TO   870
 ELSE IF (key == 0) THEN
   GO TO   600
 END IF
 570 READ (out) namex
 ASSIGN 580 TO ret
 nskip = 1
 GO TO 700
 580 kf   = kf + 1
 line = line + 1
 WRITE  (nout,590) kf,namex
 590 FORMAT (53X,i5,18X,2A4)
 IF (line-nlpp < 0) THEN
   GO TO   560
 ELSE
   GO TO   540
 END IF
 600 ASSIGN 90 TO ret
 nskip = -(kf+1)
 GO TO 700
 
!     SIMULATION OF SKPFIL (OUT,NSKIP)
 
 700 IF (nskip < 0) THEN
   GO TO   720
 ELSE IF (nskip == 0) THEN
   GO TO   710
 ELSE
   GO TO   730
 END IF
 710 GO TO ret, (70,90,530,580)
 720 REWIND out
 
!     NSKIP IS THE NEGATIVE OF THE NUMBER OF FILES TO BE SKIPPED
 
 nskip = -nskip
 730 DO  ns = 1,nskip
   740 READ (out) key
   IF (key < 0) THEN
     GO TO   740
   ELSE IF (key == 0) THEN
     GO TO   760
   END IF
   750 CONTINUE
!     ICRQ = KEY - LCOR
!     IF (KEY .GT. LCOR) GO TO 9917
   READ (out) l
   GO TO 740
   760 CONTINUE
 END DO
 GO TO 710
 
 
!     ERRORS
 
 800 mm = -1
 GO TO 950
 
 810 WRITE  (nout,820) ufm,p1
 820 FORMAT (a23,' 4120, MODULE OUTPUT2 - ILLEGAL VALUE FOR FIRST ',  &
     'PARAMETER =',i20)
 line = line + 2
 GO TO 940
 830 WRITE  (nout,840) ufm,(idhdrx(kf),kf=1,7)
 840 FORMAT (a23,' 4130, MODULE OUTPUT2 - ILLEGAL TAPE CODE HEADER = ', 7A4)
 line = line + 2
 GO TO 940
 850 WRITE  (nout,860) uwm,p3x,p3
 860 FORMAT (a25,' 4131, FORTRAN TAPE ID CODE -',2A4,'- DOES NOT MATCH'  &
     ,       ' THIRD OUTPUT2 DMAP PARAMETER -',2A4,2H-.)
 line = line + 2
 GO TO 520
 870 WRITE  (nout,880) sfm
 880 FORMAT (a25,' 4115, MODULE OUTPUT2 - SHORT RECORD.')
 line = line + 2
 GO TO 940
 890 mm = -8
 INPUT = icrq
 GO TO 950
 900 WRITE  (nout,930) sfm,key
 WRITE  (nout,910) keyx
 910 FORMAT (10X,17HEXPECTED value = ,i10,1H.)
 line = line + 3
 GO TO 940
 920 WRITE  (nout,930) sfm,key
 930 FORMAT (a25,' 2190, ILLEGAL VALUE FOR KEY =',i10,1H.)
 line = line + 2
 GO TO 940
 
 940 mm = -37
 950 CALL mesage (mm,INPUT,subnam)
 
 1000 RETURN
END SUBROUTINE outpt2
