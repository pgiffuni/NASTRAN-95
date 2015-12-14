SUBROUTINE table5 (*,in,out,trl,ibuf,wrt,lfn,fn)
     
!     THIS ROUTINE IS CALLED ONLY BY OUTPT5 TO COPY A TABLE FILE IN 'IN'
!     TO AN OUPUT FILE 'OUT', BY FORTRAN WRITE, FORMATTED OR UNFORMATTED
 
!     IN,OUT = INPUT AND OUTPUT FILE, INTEGERS
!     TRL    = TRAILER OF INPUT FILE, INTEGERS
!     P4     = 0, OUTPUT FILE IS TO BE WRITTEN UNFORMATTED, BINARY, INT.
!            = 1, OUTPUT FILE IS TO BE WRITTEN FORMATTED, INTEGER
!     TI     = ARRAY TO OVERRIDE DATA TYPE OUTPUT. INTEGERS
!              SEE RULES BELOW.
!     Z,IBUF = OPEN CORE AND GINO BUFFER POINTER, INTEGER
!     WRT,LFN= ARE COMMUNICATION FLAGS BETWEEN TABLE5 AND OUTPT5
!     FN     = ARRAY FOR INPUT FILE NAME
 
!     THE FOLLOWING CONVENTIONS ARE USED FOR FORMATTED TAPE -
 
!       A   '/'+A4  FORMAT FOR BCD WORD               ( 5 BYTES)
!       AN  'I'+I9  FORMAT FOR INTEGER                (10 BYTES)
!       A 'R'+E14.7 FORMAT FOR S.P. REAL NUMBER.      (15 BYTES)
!       A 'D'+D14.7 FORMAT FOR D.P. REAL NUMBER.      (15 BYTES)
!       A 'X'+4 BLANKS IS A FILLER, AT END OF A LINE  ( 5 BYTES)
 
!       EACH RECORD IS PRECEEDED BY L5 (IN I10 FORMAT) WHERE L5 IS THE
!       TOTAL NO. OF CHARACTERS OF THIS CURRENT RECORD DIVIDED BY 5.
 
!       EACH RECORD IS WRITTEN IN MULTIPLE LINES OF 130 CHARACTERS EACH.
!       (131 CHARACTERS TO BE EXACTLY - 130 PLUS A BLANK)
 
!       ONE OR TWO FILLERS MAY ATTACH TO THE END OF A LINE TO MAKE UP
!       130 CHARACTERS. THAT IS, INTEGER AND S.P.REAL NUMBER AT THE END
!       OF A LINE WILL NOT BE SPLITTED BETWEEN TWO LINES
 
!       IF A ZERO IS PRECEEDED BY A F.P. REAL NUMBER, IT WILL BE WRITTEN
!       OUT AS A REAL ZERO (0.0), INTEGER ZERO (0) OTHERWISE.
 
!       DUE TO THE FACTS THAT FLOATING POINT ZEROS ARE ALWAYS TREATED AS
!       INTEGERS, DOUBLE PRECISION CAN NOT BE DETECTED, AND OCCATIONALLY
!       AUTOMATIC DATA TYPE CHECKING MAY ERR, THE USER CAN OVERRIDE THE
!       OUTPUT DATA FORMAT BY DEFINING TI ARRAYS WITH THE FOLLOWING
!       RULES -
 
!          EACH TI PARAMETER MUST HOLD 9 DIGITS, FROM LEFT TO RIGHT.
!               ZEROS-FILLED IF NECCESSARY.
!               TOTALLY THERE ARE 10 TI PARAMETERS. THEREFORE, THERE ARE
!               UP TO 90 CONTINUOUS DIGITS CAN BE USED.
!               (DEFAULT IS 90 ZEROS)
!          EACH DIGIT HOLDS VALUE FROM 0 THROUGH 9, VALUE
!               0 MEANS DATA TYPE WILL BE SET AUTOMATICALLY BY TABLE5
!               1 MEANS DATA TYPE IS INTEGER
!               2 MEANS DATA TYPE IS REAL, SINGLE PRECISION
!               3 MEANS DATA TYPE IS BCD WORD (4 BYTES PER WORD)
!               4 MEANS DATA TYPE IS REAL, DOUGLE PRECISION
!             5-9 HAS SPECIAL MEANING. IT MEANS THERE ARE (5-9) VALUES
!                 OF DATA TYPE DEFINED BY THE NEXT VALUE FOLLOWING.
!          EACH DIGIT IN TI, EXCEPT 5 THRU 9, DEFINES THE CORRESPODING
!               DATA TYPE IN THE TABLE BLOCK DATA, STARTING FROM THE
!               FIRST DATA WORD AND CONTINUE TO THE LAST.
!          IF TI(1) IS NEGATIVE, INTERMEDIATE STEPS IN FORMAT GENERATION
!               ARE PRINTED OUT.
!     E.G.
!     TABLE- 3  4  3.4  5.0E-3  TESTING  .6D+7  9  G  3.2  8  0.  0  4
!            12 13  14  15  28  61   88   14   44 .7D+7
!     TI   - TI(1) =-112233413, TI(2) = 212516140  OR
!            TI(1) = 604000025, TI(2) = 060400000 (7TH AND 24 WORDS ARE
!                                            D.P. AND 12TH WORD IS REAL)
!     NOTE - 2 BCD WORDS IN 'TESTING',
!            ALL OTHERS ARE 1 COMPUTER WORD PER DATA ENTRY
!            TI(2), THE LAST TI USED HERE, MUST FILL UP WITH ZEROS TO
!               MAKE UP A 9-DIGIT WORD.
 
!     TO READ THE OUTPUT FILE, USE TABLE-V SUBROUTINE AS REFERENCE
 
!     NOTE - THE FORMATTED OUTPUT FILE CAN BE VIEWED AND/OR EDITTED BY
!            THE SYSTEM EDITOR
 
!     WRITTEN BY G.CHAN/UNISYS,  1989
 
!  $MIXED_FORMATS
 
 
 , INTENT(IN OUT)                         :: *
 INTEGER, INTENT(OUT)                     :: in
 INTEGER, INTENT(IN)                      :: out
 INTEGER, INTENT(OUT)                     :: trl(7)
 INTEGER, INTENT(IN)                      :: ibuf
 INTEGER, INTENT(OUT)                     :: wrt
 INTEGER, INTENT(OUT)                     :: lfn
 INTEGER, INTENT(OUT)                     :: fn(3,1)
 IMPLICIT INTEGER (a-z)
 LOGICAL :: debug,tion,dp
 INTEGER :: NAME(2),sub(2)
 REAL :: temp(2),rz(1)
 DOUBLE PRECISION :: dtemp
 CHARACTER (LEN=10) :: FMT(30),fmti,fmtr,fmtd,fmtb,fmtx,lpren,rpren, lpri10
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm,uwm,uim
 COMMON /system/  sysbuf,nout
 COMMON /BLANK /  dummy(4),p4,ti(1)
 COMMON /zzzzzz/  z(1)
!WKBI 7/94
 COMMON /machin/  mach
 EQUIVALENCE      (z(1),rz(1)) ,  (dtemp,temp(1))
 DATA    sub   /  4HTABL,4HE5  /, debug   /   .false.      /
 DATA    fmti,    fmtr         / '1HI,I9,' ,  '1HR,E14.7,' /
 DATA    fmtb,    fmtd         / '1H/,A4,' ,  '1HD,D14.7,' /
 DATA    fmtx,    lpri10       / '1HX,4X,' ,  '(I10,'      /
 DATA    lpren,   rpren, del   / '(', '1X)',   4H),.)      /
 DATA    END,     tble         /  4H*END,      4HTBLE      /
 
 debug = .false.
 IF (ti(1) < 0) debug =.true.
 ti(1) = IABS(ti(1))
 tion  = .false.
 DO  l = 1,10
   IF (ti(l) /= 0) tion=.true.
 END DO
 IF (debug) CALL page1
 IF (debug) WRITE (nout,20)
 20 FORMAT (///5X,'*** IN TABLE5/OUTPUT5 ***')
 kore  = ibuf - 2
 
!     OPEN INPUT FILE, AND READ FILE NAME IN THE FILE HEADER RECORD
!     WRITE ONE HEADER RECORD, IN OUTPT5 MATRIX HEADER FORMAT, TO
!     OUTPUT TAPE
 
 CALL OPEN (*810,in,z(ibuf),0)
 CALL READ (*820,*830,in,NAME,2,1,kk)
 IF (debug) WRITE (nout,30) NAME
 30 FORMAT (/5X,'PROCESSING...',2A4,/)
 i = 0
 j = 1
 trl(7) = 0
 IF (p4 == 0) WRITE (out   ) i,j,j,dtemp,(trl(k),k=2,7),NAME
 IF (p4 == 1) WRITE (out,40) i,j,j,dtemp,(trl(k),k=2,7),NAME
 40 FORMAT (3I8,/,d26.17,6I8,2A4)
 
 50 IF (p4 == 1) GO TO 100
 
!     UNFORMATED WRITE
 
 j = 2
 60 CALL READ  (*700,*70,in,z(j),kore,1,kk)
 j = 0
 GO TO 840
 70 IF (j == 1) GO TO 80
 j = 1
 z(1) = kk
 80 CALL WRITE (out,z(1),kk,1)
 GO TO 60
 
!     FORMATTED WRITE
 
 100 j = 2
 CALL READ (*700,*110,in,z(j),kore,1,kk)
 j = 0
 GO TO 840
 
!     SET UP USER DIRECTED TI TABLE IN Z(KK2) THRU Z(KK3)
 
 110 IF (debug) WRITE (nout,120) (ti(j),j=1,10)
 120 FORMAT (//5X,'TI PARAMETERS =',/4X,10(1X,i9))
 kk1 = kk  + 2
 kk2 = kk1 + kk
 kk3 = kk2 + kk
 j   = kore - kk3 - 9
 IF (j < 0) GO TO 840
 DO  k = kk1,kk3
   z(k) = 0
 END DO
 IF (.NOT.tion) GO TO 260
 k  = kk1 - 9
 ll = 0
 l  = -1
 150 IF (l >= 0) GO TO 170
 l  = 8
 ll = ll + 1
 k  = k  + 9
 IF (k >= kk2 .OR. ll > 10) GO TO 200
 til= ti(ll)
 IF (til > 0) GO TO 170
 l  = -1
 GO TO 150
 170 til10 = til/10
 z(k+l)= til - til10*10
 til   = til10
 l  = l - 1
 GO TO 150
 
 200 k  = kk2 - 1
 IF (debug) WRITE (nout,210) (z(j),j=kk1,k)
 210 FORMAT (//5X,'DIGITIZED TI PARAMTERS =',/,(3X,25I3))
 i  = kk2
 DO  j = kk1,k
   jz = z(j)
   IF (jz <= 4) GO TO 230
   ji = jz + i - 1
   jj = z(j+1)
   IF (jj > 4) GO TO 860
   DO  l = i,ji
     z(l) = jj
   END DO
   i  = ji + 1
   z(j+1) = -1
   CYCLE
   230 IF (jz == -1) CYCLE
   z(i) = jz
   i  = i + 1
 END DO
 i  = kk3 - 1
 IF (debug) WRITE (nout,250) (z(j),j=kk2,i)
 250 FORMAT (//,5X,'DECODED TI PARAMETERS =',/,(3X,25I3))
 
!     COUNT HOW MANY 5-BYTE WORDS TO BE GENERATED, FILLERS INCLUDED
 
 260 kk2 = kk2 - 1
 k   = kk1
 pjj = 1
 l5  = 10
 
 IF (debug) CALL page1
 DO  i = 1,kk
   k   = k + 1
   pjj = jj
   IF (tion) GO TO 290
   280 jj  = numtyp(z(i+1)) + 1
   GO TO 300
   290 jj  = z(kk2+i) + 1
   IF (jj == 1) GO TO 280
   300 SELECT CASE ( jj )
     CASE (    1)
       GO TO 310
     CASE (    2)
       GO TO 320
     CASE (    3)
       GO TO 340
     CASE (    4)
       GO TO 380
     CASE (    5)
       GO TO 340
   END SELECT
!              0,  I,  R,  B,  D
   
!     ZERO
   
   310 jj  = 3
   IF (pjj == 3 .OR. pjj == 5) GO TO 340
   jj  = 2
   
!     INTEGER
   
   320 IF (MOD(l5,130) <= 120) GO TO 330
   z(k)= 6
   k   = k  + 1
   l5  = l5 + 5
   330 z(k)= jj
   l5  = l5 + 10
   CYCLE
   
!     REAL, S.P. OR D.P.
   
   340 j   = MOD(l5,130)
   IF (j-120 < 0) THEN
     GO TO   370
   ELSE IF (j-120 == 0) THEN
     GO TO   350
   ELSE
     GO TO   360
   END IF
   350 l5  = l5 + 5
   z(k)= 6
   k   = k  + 1
   360 l5  = l5 + 5
   z(k)= 6
   k   = k + 1
   370 z(k)= jj
   l5  = l5 + 15
   CYCLE
   
!     BCD
   
   380 z(k)= jj
   l5  = l5 + 5
   
 END DO
 
!     NOW WE FORM THE FORMAT
 
 dp  = .false.
 kk  = k
 z(1) = (l5-10)/5
 FMT(1) = lpri10
 
 l5  = 10
 l   = 1
 i   = 1
 ib  = 1
 k   = kk1
 500 IF (l5 < 130) GO TO 540
 l   = l + 1
 FMT(l) = rpren
 IF (.NOT.debug) GO TO 520
 CALL page2 (-5)
 WRITE  (nout,510) (FMT(j),j=1,l)
 510 FORMAT (/,' DYNAMICALLY GENERATED FORMAT =',/,(1X,7A10))
!WKBD 7/94   520 WRITE  (OUT,FMT,ERR=530) (RZ(J),J=IB,I)
!WKBNB 7/94
 520 IF ( mach /= 5 .AND. mach /= 2 ) GO TO 525
 WRITE  (out,FMT,ERR=530) (rz(j),j=ib,i)
 GO TO 530
 525 isave = nout
 nout  = out
 CALL forwrt ( FMT, rz(ib), i-ib+1)
 nout  = isave
!WKBNE 7/94
 530 ib  = i + 1
 l5  = 0
 l   = 1
 FMT(1) = lpren
 
 540 k   = k + 1
 IF (k > kk) GO TO 650
 i   = i + 1
 l   = l + 1
 j   = z(k)
 SELECT CASE ( j )
   CASE (    1)
     GO TO 600
   CASE (    2)
     GO TO 600
   CASE (    3)
     GO TO 610
   CASE (    4)
     GO TO 620
   CASE (    5)
     GO TO 630
   CASE (    6)
     GO TO 640
 END SELECT
!              0,  I,  R,  B,  D, FL
 600 FMT(l) = fmti
 l5  = l5 + 10
 GO TO 500
 
!     S.P. REAL NUMBERS
 
 610 FMT(l) = fmtr
 l5  = l5 + 15
 GO TO 500
 
 620 FMT(l) = fmtb
 l5  = l5 + 5
 GO TO 500
 
!     D.P. NUMBERS
 
 630 FMT(l) = fmtd
 l5     = l5 + 15
 temp(1)= rz(l  )
 temp(2)= rz(l+1)
 z(l  ) = SNGL(dtemp)
 z(l+1) = del
 dp     =.true.
 GO TO 500
 
!     FILLER
 
 640 FMT(l) = fmtx
 l5 = l5 + 5
 i  = i  - 1
 GO TO 500
 
 650 l  = l + 1
 FMT(l) = rpren
 IF (.NOT.debug) GO TO 660
 CALL page2 (-5)
 WRITE (nout,510) (FMT(j),j=1,l)
 
!     REMOVED SECOND HALVES OF ALL D.P. NUMBERS IF THEY ARE PRESENT
!     THEN WRITE THE ARRAY OUT WITH THE GENERATED FORMAT
 
 660 IF (.NOT.dp) GO TO 680
 k   = ib - 1
 DO  j = ib,i
   IF (z(j) == del) CYCLE
   k   = k + 1
   z(k)= z(j)
 END DO
 i   = k
!WKBD 7/94  680 WRITE (OUT,FMT,ERR=690) (RZ(J),J=IB,I)
!WKBNB 7/94
 680 IF ( mach /= 2 .AND. mach /= 5 ) GO TO 685
 WRITE (out,FMT,ERR=690) (rz(j),j=ib,i)
 GO TO 690
 685 isave = nout
 nout  = out
 CALL forwrt ( FMT, rz(ib), i-ib+1)
 nout  = isave
!WKBNE 7/94
 
!     RETURN TO PROCESS ANOTHER RECORD ON INPUT FILE
 
 690 debug = .false.
 GO TO 50
 
!     ALL DONE. SET WRT FLAG, UPDATE LFN AND FN, AND CLOSE INPUT FILE
!     AND ECHO USER MESSAGES
 
 700 wrt = 1
 IF (lfn < 0) lfn = 0
 lfn = lfn + 1
 fn(1,lfn) = NAME(1)
 fn(2,lfn) = NAME(2)
 fn(3,lfn) = tble
 CALL CLOSE (in,1)
 IF (p4 == 1) GO TO 730
 CALL page2 (-7)
WRITE  (out) i,END
WRITE  (nout,710) uim,NAME
710 FORMAT (a29,' FROM OUTPUT5 MODULE, SUCCESSFUL TABLE-DATA ',  &
    'TRANSFERED FROM INPUT FILE ',2A4,' TO OUTPUT TAPE', //5X,  &
    'A HEADER RECORD WAS FIRST WRITTEN, THEN FOLLOWED BY')
WRITE  (nout,720)
720 FORMAT (5X,'FORTRAN UNFORMATTED (BINARY) WRITE')
GO TO 950
730 i  = 1
WRITE  (out,740) i,END
740 FORMAT (1X,i9,1X,a4)
CALL page2 (-13)
WRITE  (nout,710) uim,NAME
WRITE  (nout,750)
750 FORMAT (5X,'FORTRAN FORMATTED WRITE, 130 CHARACTERS PER LINE -',  &
    /10X,'(''/'',A4 FOR BCD WORD       ( 5 BYTES)',  &
    /11X,'''I'',I9 FOR INTEGER        (10 BYTES)',  &
    /11X,'''R'',E14.7 FOR S.P. REAL   (15 BYTES)',  &
    /11X,'''D'',D14.7 FOR D.P. NUMBER (15 BYTES)',  &
    /11X,'''X    '', FOR FILLER       ( 5 BYTES)')
GO TO 950

!     ERROR

810 j = 1
GO TO 850
820 j = 2
GO TO 850
830 j = 3
GO TO 850
840 in= j
j = 8
850 CALL mesage (j,in,sub)
GO TO 880
860 WRITE  (nout,870) uwm,ji,jj
870 FORMAT (a25,', OUTPTT5 MODULE PARAMETER ERROR.  WRONG INDEX ',  &
    'VALUES',2I3)
880 CALL fname (in,NAME)
WRITE  (nout,890) NAME
890 FORMAT (/5X,'TABLE DATA BLOCK ',2A4,' WAS NOT COPIED TO OUTPUT', ' TAPE')
900 CALL fwdrec (*950,in)
GO TO 900

950 RETURN 1
END SUBROUTINE table5
