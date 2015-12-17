SUBROUTINE eqmck
     
!     EQMCK CREATES AN OUTPUT FILE OF MPC CONSTRAINT FORCES AND AN
!     OVERALL TOTAL OF FORCES AND MOMENTS ON THE MODEL TO PROVIDE AN
!     EQUILIBRIUM CHECK.
!     VALID ONLY FOR STATICS AND REAL EIGENVALUE ANALYSIS.
 
!     DMAP CALLING SEQUENCE (DEFAULT PARAMETERS SHOWN)
!                                                        LAMA
!     EQMCK  CASECC,EQEXIN,GPL,BGPDT,SIL,USET,KGG,GM,UGV,PGG,QG,CSTM/
!            OQM1/V,Y,OPT=0/V,Y,GRDPNT=-1/V,N,NSKIP/V,Y,SUBNAM $
!     WHERE
!     OPT .EQ. 0, CREATE OQM
!         .LT. 0, CALCULATE ST1
!         .GT. 0, CALCULATE ST1 AND CREATES OQM
!     GRDPNT - POINT ABOUT WHICH EQUILIBRIUM IS CALCULATED.
!     NSKIP  - NO. RECORDS TO SKIP ON APPENDED FILES (1 OR GREATER),
!              NEGATIVE IF EIGENVALUE PROBLEM.
!     SUBNAM - RESERVED FOR FUTURE USE
 
 INTEGER :: bgpdt,cstm,casecc,eqexin,gm,gpl,NAME(2),oqm,pgg,  &
     qg,sf(7),sil,ugv,uset,parm,kfil(14),trl,sfl(7)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /BLANK / iopt,igrid,nskip
 COMMON /eqmk1 / k(21),kmpc,kload,kspc ,parm(4),trl(7)
 COMMON /system/ ksystm(80)
 EQUIVALENCE     (ksystm(2),nout),(k(1),casecc),(k(2),eqexin),  &
     (k(3),gpl),(k(4),bgpdt),(k(5),sil),(k(6),uset),  &
     (k(7),kgg),(k(8),gm),(k(9),ugv),(k(10),pgg),  &
     (k(11),qg),(k(12),cstm),(k(13),lama),(k(14),oqm), (k(15),sf(1))
 DATA     kfil / &
!          ... CASECC,EQEXIN,GPL,BGPDT,SIL,USET,KGG,GM ,UGV,PGG,QG,CSTM,  &
 101  , 102,   103,104,  105,106, 107,108,109,110,111,112, &
!          ... LAMA , OQM .....  &
 110  , 201 / 
 DATA     sfl  / 301,302,303,304,305,306,307 /
 DATA     NAME / 4HEQMC, 2HK  /
 
 oqm   = 0
 kmpc  = 0
 kspc  = 0
 kload = 0
 parm(3) = NAME(1)
 parm(4) = NAME(2)
 DO  i = 1,7
   sf(i) = sfl(i)
 END DO
 
 DO  i = 1,11
   trl(1) = kfil(i)
   CALL rdtrl (trl)
   k(i) = trl(1)
 END DO
 lama = k(10)
 cstm = kfil(12)
 
!     ALWAYS NECESSARY FILES
 
 parm(2) = kfil(1)
 IF (casecc < 0) GO TO 120
 parm(2) = kfil(2)
 IF (eqexin < 0) GO TO 120
 parm(2) = kfil(13)
 IF (nskip < 0 .AND. lama < 0) GO TO 120
 
!     FILES FOR OQM
 
 l = 0
 IF (iopt < 0) GO TO 40
 IF (gpl < 0 .OR. sil < 0 .OR. uset < 0) GO TO 20
 oqm = kfil(14)
 
!     MPC CONSTRAINTS
 
 20 IF (gm < 0 .OR. ugv < 0 .OR. kgg < 0) GO TO 30
 kmpc = 1
 30 IF (kmpc <= 0 .OR. iopt < 0) oqm = -kfil(14)
 IF (oqm  > 0) GO TO 40
 IF (iopt < 0) GO TO 40
 CALL page2 (2)
 WRITE  (nout,35) uwm,NAME
 35 FORMAT (a25,' 2370, MULTI-POINT CONSTRAINT FORCES NOT CALCULATED',  &
     ' IN ',a4,a2,' DUE TO MISSING INPUT FILE.')
 IF (iopt == 0) GO TO 70
 
!     ST1 CALCULATION
 
!WKBD 11/93 SPR93007   40 IF (IOPT  .EQ. 0) GO TO 60
!WKBD 11/93 SPR93007      IF (BGPDT .LT. 0) GO TO 50
!WKBI 11/93 SPR93007
 40 CONTINUE
 IF (pgg >= 0 .AND. nskip >= 0) kload  = 1
 IF (qg >= 0) kspc = 1
 l = kspc + kmpc + kload
!WKBNB 11/93 SPR93007
 IF (iopt  == 0) GO TO 60
 IF (bgpdt < 0) GO TO 50
!WKBNE 11/93 SPR93007
 IF (l > 0) GO TO 60
 50 CALL page2 (2)
 WRITE (nout,110) uwm,NAME
 IF (iopt < 0) GO TO 70
 iopt = 0
 
 60 CONTINUE
 IF (iopt < 0 .AND. l  == 0) GO TO 70
 IF (iopt == 0 .AND. oqm <= 0) GO TO 70
 IF (iopt > 0 .AND. l == 0 .AND. oqm <= 0) GO TO 70
 
!     CREATE MPC DATA AND OQM
 
 IF (kmpc > 0 .OR. (nskip > 1 .AND. l > 0)) CALL eqmckm
 IF (l == 0) GO TO 70
 
!     CALCULATE D-T FOR ST1
 
 i = igrid
 CALL eqmcka (igrid,bgpdt,cstm,eqexin,sf(2),l)
 IF (igrid /= 0) igrid = i
 IF (l == 0) GO TO 140
 
!     CALCULATE AND OUTPUT ST1
 
 CALL eqmcks
 
 70 RETURN
 
!     ERROR MESSAGES
 
 110 FORMAT (a25,' 2371, EQUILIBRIUM FORCES NOT CALCULATED IN ',a4,a2,  &
     ' DUE TO MISSING INPUT FILE.')
 120 parm(1) = 1
 CALL mesage (parm(1),parm(2),parm(3))
 GO TO 70
 
 140 CALL page2 (2)
 WRITE  (nout,150) uwm,NAME
 150 FORMAT (a25,' 2372, ',a4,a2,' IS UNABLE TO CALCULATE RIGID BODY ',  &
     'TRANSFORMATION FOR SCALAR MODEL.')
 GO TO 70
END SUBROUTINE eqmck
