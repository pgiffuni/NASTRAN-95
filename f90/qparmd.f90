SUBROUTINE qparmd
     
!     MODULE PARAMD PERFORMS THE FOLLOW OP ON PARAMETERS IN DOUBLE
!     PRECISION
!     (REFERENCE - MODULE PARAMR AND SUBROUTINE QPARMR)
 
!     DMAP
!     PARAMD  / /C,N,OP/ V,N,OUTD/V,N,IND1/V,N,IND2/
!                        V,N,OUTC/V,N,INC1/V,N,INC2/V,N,FLAG $
 
!         OP        COMPUTE
!         --        -------------------------------------------
!      BY DEFAULT   FLAG = 0
!      1  ADD       OUTD = IND1 + IND2
!      2  SUB       OUTD = IND1 - IND2
!      3  MPY       OUTD = IND1 * IND2
!      4  DIV       OUTD = IND1 / IND2 (IF IND2 = 0, FLAG IS SET TO +1)
!      5  NOP       RETURN
!      6  SQRT      OUTD = DSQRT(IND1)
!      7  SIN       OUTD = DSIN(IND1) WHERE IND1 IS IN RADIANS
!      8  COS       OUTD = DCOS(IND1) WHERE IND1 IS IN RADIANS
!      9  ABS       OUTD = DABS(IND1)
!     10  EXP       OUTD = DEXP(IND1)
!     11  TAN       OUTD = DTAN(IND1) WHERE IND1 IS IN RADIANS
!     12  ADDC      OUTC = INC1 + INC2
!     13  SUBC      OUTC = INC1 - INC2
!     14  MPYC      OUTC = INC1 * INC2
!     15  DIVC      OUTC = INC1 / INC2 (IF INC2 = 0, FLAG IS SET TO +1)
!     16  COMPLEX   OUTC = (IND1,IND2)
!     17  CSQRT     OUTC = DCSQRT(INC1)
!     18  NORM      OUTD = DSQRT(OUTC(1)**2 + OUTC(2)**2)
!     19  REAL      IND1 = OUTC(1),  IND2 = OUTC(2)
!     20  POWER     OUTD = IND1**IND2
!     21  CONJ      OUTC = DCONJG(INC1)
!     22  EQ        FLAG =-1 IF IND1 COMPARES WITH IND2
!     23  GT        FLAG =-1 IF IND1 IS GT IND2
!     24  GE        FLAG =-1 IF IND1 IS GE IND2
!     25  LT        FLAG =-1 IF IND1 IS LT IND2
!     26  LE        FLAG =-1 IF IND1 IS LE IND2
!     27  NE        FLAG =-1 IF IND1 IS NE IND2
!     28  LOG       OUTD = DLOG10(IND1)
!     29  LN        OUTD = DLOG(IND1)
!     30  FIX       FLAG = OUTD
!     31  FLOAT     OUTD = FLOAT(FLAG)
 
!     NEW OP CODE ADDED IN THIS NEW VERSION, 12/1988 -
 
!     32  ERR       IF FLAG IS 0, SYSTEM NOGO FLAG IS SET TO ZERO
!                   IF FLAG IS NON-ZERO, JOB TERMINATED IF ANY PREVIOUS
!                      PARAMD (OR PARAMR) CONTAINS NON-FATAL ERROR(S)
 
 LOGICAL :: prt
 INTEGER :: op,opcode(50),flag,ivps(1),NAME(2),il(8),ilx(8), nam(2),blnk
 REAL :: temp(2)
 DOUBLE PRECISION :: outd,ind1,ind2,outc,inc1,inc2,denom,tempd
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm,uwm,uim
 COMMON /BLANK /  op(2),outd,ind1,ind2,outc(2),inc1(2),inc2(2),flag
 COMMON /xvps  /  vps(1)
 COMMON /ilxxd /  il1,il2,il3,il4,il5,il6,il7,il8
 COMMON /system/  ibuf,nout,nogo,dummy(33),ksys37
 EQUIVALENCE      (vps(1),ivps(1)), (tempd,temp(1)), (il,il1)
 DATA NAME     /  4HQPAR,4HMD  /      ,ifirst / 15  /
 DATA opcode   /  4HADD ,4HSUB ,4HMPY ,4HDIV ,4HNOP ,  &
     4HSQRT,4HSIN ,4HCOS ,4HABS ,4HEXP , 4HTAN ,4HADDC,4HSUBC,4HMPYC,4HDIVC,  &
     4HCOMP,4HCSQR,4HNORM,4HREAL,4HPOWE, 4HCONJ,4HEQ  ,4HGT  ,4HGE  ,4HLT  ,  &
     4HLE  ,4HNE  ,4HLOG ,4HLN  ,4HFIX , 4HFLOA,4HERR ,4H    ,4H    ,4H    ,  &
     4H    ,4H    ,4H    ,4H    ,4H    , 4H    ,4H    ,4H    ,4H    ,4H    ,  &
     4H    ,4H    ,4H    ,4H    ,4H    /
 DATA ilx      /  4H1ST ,4H2ND ,4H3RD ,4H4TH ,4H5TH ,  &
     4H6TH ,4H7TH ,4H8TH               /
 DATA parm,nam /  4HPARM,4H/par,3HAMD/,blnk  /4H    /
 
!     SUPPRESS ALL PARAMETER CHECK MESSAGES IF DIAG 37 IS ON
 
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
 3 FORMAT (a29,' FROM PARAMD MODULE - OP CODE = ',2A4, /5X,  &
     '(ALL PARAMD MESSAGES CAN BE SUPPRESSED BY DIAG 37)')
 
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
 6 FORMAT (22X,'UNRECOGNIZABLE OP CODE = ',a4,' (INPUT ERROR)  ',2A4)
 CALL mesage (-7,0,NAME)
 
! *******
!     D.P. REAL NUMBER FUNCTIONS
! *******
 
!     ADD
 
 10 outd = ind1 + ind2
 GO TO 600
 
!     SUBTRACT
 
 20 outd = ind1 - ind2
 GO TO 600
 
!     MULTIPLY
 
 30 outd = ind1*ind2
 GO TO 600
 
!     DIVIDE
 
 40 outd = 0.d+0
 IF (ind2 == 0.d+0) GO TO 45
 outd = ind1/ind2
 GO TO 600
 45 WRITE  (nout,47) nam
 47 FORMAT (5X,'ERROR - DIVIDED BY ZERO  ',2A4)
 ierr = 1
 flag =+1
 IF (il8 <= 0) GO TO 730
 ivps(il8) = flag
 i = il8 - 3
 IF (prt) WRITE (nout,48) ivps(i),ivps(i+1),flag
 48 FORMAT (22X,2A4,2H =,i10,'   (OUTPUT)')
 GO TO 730
 
!     NOP
 
 50 RETURN
 
!     SQUARE ROOT
 
 60 IF (ind1 < 0.d+0) GO TO 65
 outd = DSQRT(ind1)
 GO TO 650
 65 WRITE  (nout,67) nam
 67 FORMAT (5X,'ERROR - OPERATING ON A NEGATIVE NUMBER  ',2A4)
 outd = 0.d+0
 ierr = 1
 GO TO 650
 
!     SINE
 
 70 outd = DSIN(ind1)
 GO TO 650
 
!     COSINE
 
 80 outd = DCOS(ind1)
 GO TO 650
 
!     ABSOLUTE VALUE
 
 90 outd = DABS(ind1)
 GO TO 650
 
!     EXPONENTIAL
 
 100 outd = DEXP(ind1)
 GO TO 650
 
!     TANGENT
 
 110 outd = DTAN(ind1)
 GO TO 650
 
!     NORM
 
 180 outd = DSQRT(outc(1)**2 + outc(2)**2)
 GO TO 690
 
!     POWER
 
 200 outd  = ind1**ind2
 GO TO 600
 
!     LOG
 
 280 IF (ind1 < 0.d+0) GO TO 65
 outd = DLOG10(ind1)
 GO TO 650
 
!     NATURAL LOG
 
 290 IF (ind1 < 0.d+0) GO TO 65
 outd = DLOG(ind1)
 GO TO 650
 
!     FLOAT
 
 310 outd = iflag
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
 
 120 outc(1) = inc1(1) + inc2(1)
 outc(2) = inc1(2) + inc2(2)
 GO TO 730
 
!     SUBTRACT COMPLEX
 
 130 outc(1) = inc1(1) - inc2(1)
 outc(2) = inc1(2) - inc2(2)
 GO TO 730
 
!     MULTIPLY COMPLEX
 
 140 outc(1) = inc1(1)*inc2(1) - inc1(2)*inc2(2)
 outc(2) = inc1(1)*inc2(2) + inc1(2)*inc2(1)
 GO TO 730
 
!     DIVIDE COMPLEX
 
 150 denom = inc2(1)**2 + inc2(2)**2
 IF (denom == 0.d+0) GO TO 155
 outc(1) = (inc1(1)*inc2(1) + inc1(2)*inc2(2))/denom
 outc(2) = (inc1(2)*inc2(1) - inc1(1)*inc2(2))/denom
 GO TO 730
 155 outc(1) = 0.d+0
 outc(2) = 0.d+0
 GO TO 45
 
!     COMPLEX
 
 160 outc(1) = ind1
 outc(2) = ind2
 GO TO 710
 
!     COMPLEX SQUARE ROOT
 
 170 outc(1) = (inc1(1)**2 + inc1(2)**2)**0.25D0  &
     *DCOS(0.5D0*DATAN2(inc1(2),inc1(1)))
 outc(2) = (inc1(1)**2 + inc1(2)**2)**0.25D0  &
     *DSIN(0.5D0*DATAN2(inc1(2),inc1(1)))
 GO TO 760
 
!     CONJUGATE
 
 210 outc(1) = inc1(1)
 outc(2) =-inc1(2)
 GO TO 760
 
!     REAL
 
 190 ind1 = outc(1)
 ind2 = outc(2)
 GO TO 770
 
!     EQUAL
 
 220 IF (ind1 == ind2) flag = -1
 GO TO 660
 
!     GREATER THAN
 
 230 IF (ind1 > ind2) flag = -1
 GO TO 660
 
!     GREATER THAN OR EQUAL
 
 240 IF (ind1 >= ind2) flag = -1
 GO TO 660
 
!     LESS THAN
 
 250 IF (ind1 < ind2) flag = -1
 GO TO 660
 
!     LESS THAN OR EQUAL
 
 260 IF (ind1 <= ind2) flag = -1
 GO TO 660
 
!     NOT EQUAL
 
 270 IF (ind1 /= ind2) flag = -1
 GO TO 660
 
!     FIX
 
 300 flag = outd
 GO TO 720
 
! ---------------------------------------------------
 
!     INPUT PARAMETER ECHO
 
 600 ASSIGN 620 TO irtn3
 ASSIGN 800 TO irtn4
 610 IF (.NOT.prt) GO TO 615
 i = il3 - 3
 IF (il3 <= 0) WRITE (nout,640) ilx(3),parm,ind1
 IF (il3 > 0) WRITE (nout,640) ivps(i),ivps(i+1),ind1
 615 IF (il3 == 0) ierr = 1
 GO TO irtn3, (620,800)
 620 IF (.NOT.prt) GO TO 645
 j = il4 - 3
 IF (il4 <= 0) WRITE (nout,640) ilx(4),parm,ind2
 IF (il4 > 0) WRITE (nout,640) ivps(j),ivps(j+1),ind2
 640 FORMAT (22X,2A4,3H = ,d15.8,'  (INPUT)')
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
 700 FORMAT (22X,2A4,4H = (,d15.8,1H,,d15.8,')   (INPUT)')
 705 IF (il5 == 0) ierr = 1
 GO TO 800
 
 710 ASSIGN 620 TO irtn3
 ASSIGN 880 TO irtn4
 GO TO 610
 
 720 IF (.NOT.prt) GO TO 725
 i = il2 - 3
 IF (il2 <= 0) WRITE (nout,640) ilx(2),parm,outd
 IF (il2 > 2) WRITE (nout,640) ivps(i),ivps(i+1),outd
 725 IF (il2 == 0) ierr = 1
 GO TO 910
 
 730 ASSIGN 750 TO irtn6
 740 IF (.NOT.prt) GO TO 745
 i = il6 - 3
 IF (il6 <= 0) WRITE (nout,700) ilx(6),parm,inc1
 IF (il6 > 0) WRITE (nout,700) ivps(i),ivps(i+1),inc1
 IF (il6 == 0) ierr = 1
 745 GO TO irtn6, (750,880)
 750 IF (.NOT.prt) GO TO 755
 j = il7 - 3
 IF (il7 <= 0) WRITE (nout,700) ilx(7),parm,inc2
 IF (il7 > 0) WRITE (nout,700) ivps(j),ivps(j+1),inc2
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
 
!     SECOND PARAMETER - OUTD
 
 800 IF (il2  >  0) GO TO 820
 WRITE  (nout,810) ilx(2),nam
 810 FORMAT (22X,a4,'PARAMETER IS MISSING    (OUTPUT ERROR)  ',2A4)
 ierr = 1
 GO TO 950
 820 IF (ierr == 0) GO TO 825
 temp(1) = vps(il2  )
 temp(2) = vps(il2+1)
 outd  = tempd
 825 tempd = outd
 vps(il2  ) = temp(1)
 vps(il2+1) = temp(2)
 i = il2 - 3
 IF (prt) WRITE (nout,830) ivps(i),ivps(i+1),outd
 830 FORMAT (22X,2A4,3H = ,d15.8,'  (OUTPUT)')
 GO TO 950
 
!     THIRD AND FOURTH PARAMETERS - IND1, IND2
 
 840 IF (il3 >  0) GO TO 850
 WRITE (nout,810) ilx(3),nam
 ierr = 1
 GO TO 860
 850 IF (ierr == 0) GO TO 855
 temp(1) = vps(il3  )
 temp(2) = vps(il3+1)
 ind1  = tempd
 855 tempd = ind1
 vps(il3  ) = temp(1)
 vps(il3+1) = temp(2)
 i = il3 - 3
 IF (prt) WRITE (nout,830) ivps(i),ivps(i+1),ind1
 860 IF (il4 >  0) GO TO 870
 WRITE (nout,810) ilx(4),nam
 ierr = 1
 GO TO 950
 870 IF (ierr == 0) GO TO 875
 temp(1) = vps(il4  )
 temp(2) = vps(il4+1)
 ind2  = tempd
 875 tempd = ind2
 vps(il4  ) = temp(1)
 vps(il4+1) = temp(2)
 j = il4 - 3
 IF (prt) WRITE (nout,830) ivps(j),ivps(j+1),ind2
 GO TO 950
 
!     FIFTH PARAMETER - OUTC
 
 880 IF (il5 >  0) GO TO 890
 WRITE (nout,810) ilx(5),nam
 ierr = 1
 GO TO 950
 890 IF (ierr == 0) GO TO 895
 temp(1) = vps(il5  )
 temp(2) = vps(il5+1)
 outc(1) = tempd
 temp(1) = vps(il5+2)
 temp(2) = vps(il5+3)
 outc(2) = tempd
 895 tempd = outc(1)
 vps(il5  ) = temp(1)
 vps(il5+1) = temp(2)
 tempd = outc(2)
 vps(il5+2) = temp(1)
 vps(il5+3) = temp(2)
 i = il5 - 3
 IF (prt) WRITE (nout,900) ivps(i),ivps(i+1),outc
 900 FORMAT (22X,2A4,4H = (,d15.8,1H,,d15.8,')   (OUTPUT)')
 GO TO 950
 
!     EIGHTH PARAMETER - FLAG
 
 910 IF (il8 >  0) GO TO 920
 WRITE (nout,810) ilx(8),nam
 ierr = 1
 GO TO 950
 920 IF (ierr == 0) ivps(il8) = flag
 i = il8 - 3
 IF (prt) WRITE (nout,930) ivps(i),ivps(i+1),ivps(il8)
 930 FORMAT (22X,2A4,2H =,i12,6X,'(OUTPUT)')
 
 950 IF (ierr ==  0) GO TO 990
 WRITE  (nout,960) uwm,nam
 960 FORMAT (a25,' - I/O ERROR, OUTPUT NOT SAVED. OUTPUT DEFAULT ',  &
     'VALUE REMAINS ',2A4,/)
 GO TO 990
 970 WRITE  (nout,980) nam
 980 FORMAT (5X,'JOB TERMINTATED DUE TO PREVIOUS ERROR(S)  ',2A4,/)
 CALL pexit
 990 IF (ksys37 == 0) ksys37 = ierr
 RETURN
 
END SUBROUTINE qparmd
