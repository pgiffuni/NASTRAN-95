      SUBROUTINE EQMCK
C
C     EQMCK CREATES AN OUTPUT FILE OF MPC CONSTRAINT FORCES AND AN
C     OVERALL TOTAL OF FORCES AND MOMENTS ON THE MODEL TO PROVIDE AN
C     EQUILIBRIUM CHECK.
C     VALID ONLY FOR STATICS AND REAL EIGENVALUE ANALYSIS.
C
C     DMAP CALLING SEQUENCE (DEFAULT PARAMETERS SHOWN)
C                                                        LAMA
C     EQMCK  CASECC,EQEXIN,GPL,BGPDT,SIL,USET,KGG,GM,UGV,PGG,QG,CSTM/
C            OQM1/V,Y,OPT=0/V,Y,GRDPNT=-1/V,N,NSKIP/V,Y,SUBNAM $
C     WHERE
C     OPT .EQ. 0, CREATE OQM
C         .LT. 0, CALCULATE ST1
C         .GT. 0, CALCULATE ST1 AND CREATES OQM
C     GRDPNT - POINT ABOUT WHICH EQUILIBRIUM IS CALCULATED.
C     NSKIP  - NO. RECORDS TO SKIP ON APPENDED FILES (1 OR GREATER),
C              NEGATIVE IF EIGENVALUE PROBLEM.
C     SUBNAM - RESERVED FOR FUTURE USE
C
      INTEGER         BGPDT,CSTM,CASECC,EQEXIN,GM,GPL,NAME(2),OQM,PGG,
     1                QG,SF(7),SIL,UGV,USET,PARM,KFIL(14),TRL,SFL(7)
      CHARACTER       UFM*23,UWM*25
      COMMON /XMSSG / UFM,UWM
      COMMON /BLANK / IOPT,IGRID,NSKIP
      COMMON /EQMK1 / K(21),KMPC,KLOAD,KSPC ,PARM(4),TRL(7)
      COMMON /SYSTEM/ KSYSTM(80)
      EQUIVALENCE     (KSYSTM(2),NOUT),(K(1),CASECC),(K(2),EQEXIN),
     1                (K(3),GPL),(K(4),BGPDT),(K(5),SIL),(K(6),USET),
     2                (K(7),KGG),(K(8),GM),(K(9),UGV),(K(10),PGG),
     3                (K(11),QG),(K(12),CSTM),(K(13),LAMA),(K(14),OQM),
     4                (K(15),SF(1))
      DATA     KFIL /
C          ... CASECC,EQEXIN,GPL,BGPDT,SIL,USET,KGG,GM ,UGV,PGG,QG,CSTM,
     1         101  , 102,   103,104,  105,106, 107,108,109,110,111,112,
C          ... LAMA , OQM .....
     2         110  , 201 /
      DATA     SFL  / 301,302,303,304,305,306,307 /
      DATA     NAME / 4HEQMC, 2HK  /
C
      OQM   = 0
      KMPC  = 0
      KSPC  = 0
      KLOAD = 0
      PARM(3) = NAME(1)
      PARM(4) = NAME(2)
      DO 5 I = 1,7
    5 SF(I) = SFL(I)
C
      DO 10 I = 1,11
      TRL(1) = KFIL(I)
      CALL RDTRL (TRL)
      K(I) = TRL(1)
   10 CONTINUE
      LAMA = K(10)
      CSTM = KFIL(12)
C
C     ALWAYS NECESSARY FILES
C
      PARM(2) = KFIL(1)
      IF (CASECC .LT. 0) GO TO 120
      PARM(2) = KFIL(2)
      IF (EQEXIN .LT. 0) GO TO 120
      PARM(2) = KFIL(13)
      IF (NSKIP.LT.0 .AND. LAMA.LT.0) GO TO 120
C
C     FILES FOR OQM
C
      L = 0
      IF (IOPT .LT. 0) GO TO 40
      IF (GPL.LT.0 .OR. SIL.LT.0 .OR. USET.LT.0) GO TO 20
      OQM = KFIL(14)
C
C     MPC CONSTRAINTS
C
   20 IF (GM.LT.0 .OR. UGV.LT.0 .OR. KGG.LT.0) GO TO 30
      KMPC = 1
   30 IF (KMPC.LE.0 .OR. IOPT.LT.0) OQM = -KFIL(14)
      IF (OQM  .GT. 0) GO TO 40
      IF (IOPT .LT. 0) GO TO 40
      CALL PAGE2 (2)
      WRITE  (NOUT,35) UWM,NAME
   35 FORMAT (A25,' 2370, MULTI-POINT CONSTRAINT FORCES NOT CALCULATED',
     1       ' IN ',A4,A2,' DUE TO MISSING INPUT FILE.')
      IF (IOPT .EQ. 0) GO TO 70
C
C     ST1 CALCULATION
C
CWKBD 11/93 SPR93007   40 IF (IOPT  .EQ. 0) GO TO 60
CWKBD 11/93 SPR93007      IF (BGPDT .LT. 0) GO TO 50
CWKBI 11/93 SPR93007
   40 CONTINUE
      IF (PGG.GE.0 .AND. NSKIP.GE.0) KLOAD  = 1
      IF (QG .GE. 0) KSPC = 1
      L = KSPC + KMPC + KLOAD
CWKBNB 11/93 SPR93007
      IF (IOPT  .EQ. 0) GO TO 60
      IF (BGPDT .LT. 0) GO TO 50
CWKBNE 11/93 SPR93007
      IF (L .GT. 0) GO TO 60
   50 CALL PAGE2 (2)
      WRITE (NOUT,110) UWM,NAME
      IF (IOPT .LT. 0) GO TO 70
      IOPT = 0
C
   60 CONTINUE
      IF (IOPT.LT.0 .AND. L  .EQ.0) GO TO 70
      IF (IOPT.EQ.0 .AND. OQM.LE.0) GO TO 70
      IF (IOPT.GT.0 .AND. L.EQ.0 .AND. OQM.LE.0) GO TO 70
C
C     CREATE MPC DATA AND OQM
C
      IF (KMPC.GT.0 .OR. (NSKIP.GT.1 .AND. L.GT.0)) CALL EQMCKM
      IF (L .EQ. 0) GO TO 70
C
C     CALCULATE D-T FOR ST1
C
      I = IGRID
      CALL EQMCKA (IGRID,BGPDT,CSTM,EQEXIN,SF(2),L)
      IF (IGRID .NE. 0) IGRID = I
      IF (L .EQ. 0) GO TO 140
C
C     CALCULATE AND OUTPUT ST1
C
      CALL EQMCKS
C
   70 RETURN
C
C     ERROR MESSAGES
C
  110 FORMAT (A25,' 2371, EQUILIBRIUM FORCES NOT CALCULATED IN ',A4,A2,
     1       ' DUE TO MISSING INPUT FILE.')
  120 PARM(1) = 1
      CALL MESAGE (PARM(1),PARM(2),PARM(3))
      GO TO 70
C
  140 CALL PAGE2 (2)
      WRITE  (NOUT,150) UWM,NAME
  150 FORMAT (A25,' 2372, ',A4,A2,' IS UNABLE TO CALCULATE RIGID BODY ',
     1       'TRANSFORMATION FOR SCALAR MODEL.')
      GO TO 70
      END
