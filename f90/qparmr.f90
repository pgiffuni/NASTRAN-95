SUBROUTINE qparmr
     
!     MODULE PARAMR PERFORMS THE FOLLOW OP ON PARAMETERS IN SINGLE
!     PRECISION
!     (COMPANION MODULE PARAMD AND SUBROUTINE QPARMD)
 
!     DMAP
!     PARAMR  / /C,N,OP/ V,N,OUTR/V,N,IN1R/V,N,IN2R/
!                        V,N,OUTC/V,N,IN1C/V,N,IN2C/V,N,FLAG $
 
!         OP        COMPUTE
!         --        -------------------------------------------
!      BY DEFAULT   FLAG = 0
!      1  ADD       OUTR = IN1R + IN2R
!      2  SUB       OUTR = IN1R - IN2R
!      3  MPY       OUTR = IN1R * IN2R
!      4  DIV       OUTR = IN1R / IN2R (IF IN2R = 0, FLAG IS SET TO +1)
!      5  NOP       RETURN
!      6  SQRT      OUTR = SQRT(IN1R)
!      7  SIN       OUTR = SIN(IN1R) WHERE IN1R IS IN RADIANS
!      8  COS       OUTR = COS(IN1R) WHERE IN1R IS IN RADIANS
!      9  ABS       OUTR = ABS(IN1R)
!     10  EXP       OUTR = EXP(IN1R)
!     11  TAN       OUTR = TAN(IN1R) WHERE IN1R IS IN RADIANS
!     12  ADDC      OUTC = IN1C + IN2C
!     13  SUBC      OUTC = IN1C - IN2C
!     14  MPYC      OUTC = IN1C * IN2C
!     15  DIVC      OUTC = IN1C / IN2C (IF IN2C = 0, FLAG IS SET TO +1)
!     16  COMPLEX   OUTC = (IN1R,IN2R)
!     17  CSQRT     OUTC = CSQRT(IN1C)
!     18  NORM      OUTR = SQRT(OUTC(1)**2 + OUTC(2)**2)
!     19  REAL      IN1R = OUTC(1),   IN2R = OUTC(2)
!     20  POWER     OUTR = IN1R**IN2R
!     21  CONJ      OUTC = CONJG(IN1C)
!     22  EQ        FLAG =-1 IF IN1R COMPARES WITH IN2R
!     23  GT        FLAG =-1 IF IN1R IS GT IN2R
!     24  GE        FLAG =-1 IF IN1R IS GE IN2R
!     25  LT        FLAG =-1 IF IN1R IS LT IN2R
!     26  LE        FLAG =-1 IF IN1R IS LE IN2R
!     27  NE        FLAG =-1 IF IN1R IS NE IN2R
!     28  LOG       OUTR = ALOG10(IN1R)
!     29  LN        OUTR = ALOG(IN1R)
!     30  FIX       FLAG = OUTR
!     31  FLOAT     OUTR = FLOAT(FLAG)
 
!     NEW OP CODE ADDED IN THIS NEW VERSION, 12/1988 -
 
!     32  ERR       IF FLAG IS 0, SYSTEM NOGO FLAG IS SET TO ZERO
!                   IF FLAG IS NON-ZERO, JOB TERMINATED IF ANY PREVIOUS
!                      PARAMR (OR PARAMD) CONTAINS NON-FATAL ERROR(S)
 
 LOGICAL :: prt
 INTEGER :: op,opcode(50),flag,ivps(1),NAME(2),il(8),ilx(8), nam(2),blnk
 REAL :: in1r,in2r,in1c,in2c
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm,uwm,uim
 COMMON /BLANK /  op(2),outr,in1r,in2r,outc(2),in1c(2),in2c(2),flag
 COMMON /xvps  /  vps(1)
 COMMON /ilxxr /  il1,il2,il3,il4,il5,il6,il7,il8
 COMMON /system/  ibuf,nout,nogo,dummy(33),ksys37
 EQUIVALENCE      (vps(1),ivps(1)), (il,il1)
 DATA NAME     /  4HQPAR,4HMR  /      ,ifirst / 15  /
 DATA opcode   /  4HADD ,4HSUB ,4HMPY ,4HDIV ,4HNOP ,  &
     4HSQRT,4HSIN ,4HCOS ,4HABS ,4HEXP , 4HTAN ,4HADDC,4HSUBC,4HMPYC,4HDIVC,  &
     4HCOMP,4HCSQR,4HNORM,4HREAL,4HPOWE, 4HCONJ,4HEQ  ,4HGT  ,4HGE  ,4HLT  ,  &
     4HLE  ,4HNE  ,4HLOG ,4HLN  ,4HFIX , 4HFLOA,4HERR ,4H    ,4H    ,4H    ,  &
     4H    ,4H    ,4H    ,4H    ,4H    , 4H    ,4H    ,4H    ,4H    ,4H    ,  &
     4H    ,4H    ,4H    ,4H    ,4H    /
 DATA ilx      /  4H1ST ,4H2ND ,4H3RD ,4H4TH ,4H5TH ,  &
     4H6TH ,4H7TH ,4H8TH               /
 DATA parm,nam /  4HPARM,4H/par,3HAMR/,blnk  /4H    /
 
!     SUPPRESSED ALL INPUT/OUTPUT CHECK MESSAGES IF DIAG 37 IS ON
 
 CALL sswtch (37,i)
 prt = i == 0
 IF (prt) nam(1) = blnk
 IF (prt) nam(2) = blnk
 
!     COMPUTE VPS INDEXES AND PARAMETER NAMES
 
 DO  i = 2,8
   CALL fndpar (-i,il(i))
 END DO
 IF (.NOT.prt) GO TO 4
 CALL page2 (ifirst)
 ifirst = 6
 WRITE  (nout,3) uim,op
 3 FORMAT (a29,' FROM PARAMR MODULE - OP CODE = ',2A4, /5X,  &
     '(ALL PARAMR MESSAGES CAN BE SUPPRESED BY DIAG 37)')
 
!     BRANCH ON OPERATION CODE
 
 4 iflag = flag
 flag  = 0
 ierr  = 0
 
 DO  iop = 1,32
   IF (op(1) == opcode(iop)) GO TO  &
       (  10,  20,  30,  40,  50,  60,  70,  80,  90, 100,  &
       110, 120, 130, 140, 150, 160, 170, 180, 190, 200,  &
       210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320    ), iop
 END DO
 WRITE  (nout,6) op(1),nam
 6 FORMAT (22X,'UNRECOGNIZABLE OP CODE = ',a4,'  (INPUT ERROR) ',2A4)
 CALL mesage (-7,0,NAME)
 
! *******
!     REAL NUMBER FUNCTIONS
! *******
 
!     ADD
 
 10 outr = in1r + in2r
 GO TO 600
 
!     SUBTRACT
 
 20 outr = in1r - in2r
 GO TO 600
 
!     MULTIPLY
 
 30 outr = in1r*in2r
 GO TO 600
 
!     DIVIDE
 
 40 outr = 0.0
 IF (in2r == 0.d0) GO TO 45
 outr = in1r/in2r
 GO TO 600
 45 WRITE  (nout,47) nam
 47 FORMAT (5X,'ERROR - DIVIDED BY ZERO  ',2A4)
 ierr = 1
 flag =+1
 IF (il8 <= 0) GO TO 730
 ivps(il8) = flag
 i = il8 - 3
 WRITE  (nout,48) ivps(i),ivps(i+1),flag,nam
 48 FORMAT (22X,2A4,2H =,i10,'   (OUTPUT)  ',2A4)
 GO TO 730
 
!     NOP
 
 50 RETURN
 
!     SQUARE ROOT
 
 60 IF (in1r < 0.0) GO TO 65
 outr = SQRT(in1r)
 GO TO 650
 65 WRITE  (nout,67) nam
 67 FORMAT (5X,'ERROR - OPERATING ON A NEGATIVE NUMBER  ',2A4)
 outr = 0.0
 ierr = 1
 GO TO 650
 
!     SINE
 
 70 outr = SIN(in1r)
 GO TO 650
 
!     COSINE
 
 80 outr = COS(in1r)
 GO TO 650
 
!     ABSOLUTE VALUE
 
 90 outr = ABS(in1r)
 GO TO 650
 
!     EXPONENTIAL
 
 100 outr = EXP(in1r)
 GO TO 650
 
!     TANGENT
 
 110 outr = TAN(in1r)
 GO TO 650
 
!     NORM
 
 180 outr = SQRT(outc(1)**2 + outc(2)**2)
 GO TO 690
 
!     POWER
 
 200 outr  = in1r**in2r
 GO TO 600
 
!     LOG
 
 280 IF (in1r < 0.0) GO TO 65
 outr = ALOG10(in1r)
 GO TO 650
 
!     NATURAL LOG
 
 290 IF (in1r < 0.0) GO TO 65
 outr = ALOG(in1r)
 GO TO 650
 
!     FLOAT
 
 310 outr = iflag
 GO TO 670
 
!     ERR
 
 320 IF (iflag /= 0 .AND. ksys37 /= 0) GO TO 970
 ksys37 = 0
 nogo   = 0
 IF (prt) WRITE (nout,325)
 325 FORMAT (5X,'SYSTEM NOGO FLAG IS RESET TO INTEGER ZERO',/)
 GO TO 990
 
! *******
!     COMPLEX FUNCTIONS
! *******
 
!     ADD COMPLEX
 
 120 outc(1) = in1c(1) + in2c(1)
 outc(2) = in1c(2) + in2c(2)
 GO TO 730
 
!     SUBTRACT COMPLEX
 
 130 outc(1) = in1c(1) - in2c(1)
 outc(2) = in1c(2) - in2c(2)
 GO TO 730
 
!     MULTIPLY COMPLEX
 
 140 outc(1) = in1c(1)*in2c(1) - in1c(2)*in2c(2)
 outc(2) = in1c(1)*in2c(2) + in1c(2)*in2c(1)
 GO TO 730
 
!     DIVIDE COMPLEX
 
 150 denom = in2c(1)**2 + in2c(2)**2
 IF (denom == 0.0) GO TO 155
 outc(1) = (in1c(1)*in2c(1) + in1c(2)*in2c(2))/denom
 outc(2) = (in1c(2)*in2c(1) - in1c(1)*in2c(2))/denom
 GO TO 730
 155 outc(1) = 0.0
 outc(2) = 0.0
 GO TO 45
 
!     COMPLEX
 
 160 outc(1) = in1r
 outc(2) = in2r
 GO TO 710
 
!     COMPLEX SQUARE ROOT
 
 170 outc(1) = (in1c(1)**2 + in1c(2)**2)**0.25  &
     *COS(0.5*ATAN2(in1c(2),in1c(1)))
 outc(2) = (in1c(1)**2 + in1c(2)**2)**0.25 *SIN(0.5*ATAN2(in1c(2),in1c(1)))
 GO TO 760
 
!     CONJUGATE
 
 210 outc(1) = in1c(1)
 outc(2) =-in1c(2)
 GO TO 760
 
!     REAL
 
 190 in1r = outc(1)
 in2r = outc(2)
 GO TO 770
 
!     EQUAL
 
 220 IF (in1r == in2r) flag = -1
 GO TO 660
 
!     GREATER THAN
 
 230 IF (in1r > in2r) flag = -1
 GO TO 660
 
!     GREATER THAN OR EQUAL
 
 240 IF (in1r >= in2r) flag = -1
 GO TO 660
 
!     LESS THAN
 
 250 IF (in1r < in2r) flag = -1
 GO TO 660
 
!     LESS THAN OR EQUAL
 
 260 IF (in1r <= in2r) flag = -1
 GO TO 660
 
!     NOT EQUAL
 
 270 IF (in1r /= in2r) flag = -1
 GO TO 660
 
!     FIX
 
 300 flag = outr
 GO TO 720
 
! ---------------------------------------------------
 
!     INPUT PARAMETER ECHO
 
 600 ASSIGN 620 TO irtn3
 ASSIGN 800 TO irtn4
 610 IF (.NOT.prt) GO TO 615
 i = il3 - 3
 IF (il3 <= 0) WRITE (nout,640) ilx(3),parm,in1r
 IF (il3 > 0) WRITE (nout,640) ivps(i),ivps(i+1),in1r
 615 IF (il3 == 0) ierr = 1
 GO TO irtn3, (620,800)
 620 IF (.NOT.prt) GO TO 645
 j = il4 - 3
 IF (il4 <= 0) WRITE (nout,640) ilx(4),parm,in2r
 IF (il4 > 0) WRITE (nout,640) ivps(j),ivps(j+1),in2r
 640 FORMAT (22X,2A4,3H = ,e13.6,'  (INPUT)')
 645 IF (il4 == 0) ierr = 1
 GO TO irtn4, (800,880,910)
 
 650 ASSIGN 800 TO irtn3
 GO TO 610
 
 660 ASSIGN 620 TO irtn3
 ASSIGN 910 TO irtn4
 GO TO 610
 
 670 IF (.NOT.prt) GO TO 685
 i = il8 - 3
 IF (il8 <= 0) WRITE (nout,680) ilx(8),parm,iflag
 IF (il8 > 0) WRITE (nout,680) ivps(i),ivps(i+1),iflag
 680 FORMAT (22X,2A4,2H =,i10,'   (INPUT)')
 685 IF (il8 == 0) ierr = 1
 GO TO 800
 
 690 IF (.NOT.prt) GO TO 705
 i = il5 - 3
 IF (il5 <= 0) WRITE (nout,700) ilx(5),parm,outc
 IF (il5 > 0) WRITE (nout,700) ivps(i),ivps(i+1),outc
 700 FORMAT (22X,2A4,4H = (,e13.6,1H,,e13.6,')   (INPUT)')
 705 IF (il5 == 0) ierr = 1
 GO TO 800
 
 710 ASSIGN 620 TO irtn3
 ASSIGN 880 TO irtn4
 GO TO 610
 
 720 IF (.NOT.prt) GO TO 725
 i = il2 - 3
 IF (il2 <= 0) WRITE (nout,640) ilx(2),parm,outr
 IF (il2 > 0) WRITE (nout,640) ivps(i),ivps(i+1),outr
 725 IF (il2 == 0) ierr = 1
 GO TO 910
 
 730 ASSIGN 750 TO irtn6
 740 IF (.NOT.prt) GO TO 745
 i = il6 - 3
 IF (il6 <= 0) WRITE (nout,700) ilx(6),parm,in1c
 IF (il6 > 0) WRITE (nout,700) ivps(i),ivps(i+1),in1c
 745 IF (il6 == 0) ierr = 1
 GO TO irtn6, (750,880)
 750 IF (.NOT.prt) GO TO 755
 j = il7 - 3
 IF (il7 <= 0) WRITE (nout,700) ilx(7),parm,in2c
 IF (il7 > 0) WRITE (nout,700) ivps(j),ivps(j+1),in2c
 755 IF (il7 == 0) ierr = 1
 GO TO 880
 
 760 ASSIGN 880 TO irtn6
 GO TO 740
 
 770 IF (.NOT.prt) GO TO 775
 i = il5 - 3
 IF (il5 <= 0) WRITE (nout,700) ilx(5),parm,outc
 IF (il5 > 0) WRITE (nout,700) ivps(i),ivps(i+1),outc
 775 IF (il5 == 0) ierr = 1
 GO TO 840
 
!     OUTPUT PARAMETER CHECK
 
!     SECOND PARAMETER - OUTR
 
 800 IF (il2 > 0) GO TO 820
 WRITE  (nout,810) ilx(2),nam
 810 FORMAT (22X,a4,'PARAMETER IS MISSING  (OUTPUT ERROR)  ',2A4)
 ierr = 1
 GO TO 950
 820 IF (ierr == 0) vps(il2) = outr
 i = il2 - 3
 IF (prt) WRITE (nout,830) ivps(i),ivps(i+1),vps(il2)
 830 FORMAT (22X,2A4,3H = ,e13.6,'  (OUTPUT)')
 GO TO 950
 
!     THIRD AND FOURTH PARAMETERS - INR1, INR2
 
 840 IF (il3 > 0) GO TO 850
 WRITE (nout,810) ilx(3),nam
 ierr = 1
 GO TO 860
 850 IF (ierr == 0) vps(il3) = in1r
 i = il3 - 3
 IF (prt) WRITE (nout,830) ivps(i),ivps(i+1),vps(il3)
 860 IF (il4 > 0) GO TO 870
 WRITE (nout,810) ilx(4),nam
 ierr = 1
 GO TO 950
 870 IF (ierr == 0) vps(il4) = in2r
 j = il4 - 3
 IF (prt) WRITE (nout,830) ivps(j),ivps(j+1),vps(il4)
 GO TO 950
 
!     FIFTH PARAMETER - OUTC
 
 880 IF (il5 > 0) GO TO 890
 WRITE (nout,810) ilx(5),nam
 ierr = 1
 GO TO 950
 890 IF (ierr == 1) GO TO 895
 vps(il5  ) = outc(1)
 vps(il5+1) = outc(2)
 895 i = il5 - 3
 IF (prt) WRITE (nout,900) ivps(i),ivps(i+1),vps(il5),vps(il5+1)
 900 FORMAT (22X,2A4,4H = (,e13.6,1H,,e13.6,')   (OUTPUT)')
 GO TO 950
 
!     EIGHTH PARAMETER - FLAG
 
 910 IF (il8 >  0) GO TO 920
 WRITE (nout,810) ilx(8),nam
 ierr = 1
 GO TO 950
 920 IF (ierr == 0) ivps(il8) = flag
 i = il8 - 3
 IF (prt) WRITE (nout,930) ivps(i),ivps(i+1),ivps(il8)
 930 FORMAT (22X,2A4,2H =,i10,6X,'(OUTPUT)')
 
 950 IF (ierr  == 0) GO TO 990
 WRITE  (nout,960) uwm,nam
 960 FORMAT (a25,' - I/O ERROR, OUTPUT NOT SAVED. OUTPUT DEFAULT ',  &
     'VALUE REMAINS ',2A4,/)
 GO TO 990
 970 WRITE  (nout,980)
 980 FORMAT (5X,'JOB TERMINATED DUE TO PREVIOUS ERROR(S)',/)
 CALL pexit
 990 IF (ksys37 == 0) ksys37 = ierr
 RETURN
 
END SUBROUTINE qparmr
