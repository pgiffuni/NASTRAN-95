SUBROUTINE semint (debug1)
     
    !     SEMINT IS THE EXECUTION MONITOR FOR THE PREFACE.
    !     UMF IS NO LONGER SUPPORTED.
 
    !     FOR DEBUG PURPOSE, PRINT OUT GOES TO UNIT 6, NOT OUTTAP
 
    INTEGER, INTENT(IN OUT)                  :: debug1
    INTEGER :: axic,axif,outtap,plotf,hicore
    CHARACTER (LEN=1) :: ufm*23,uwm*25,uim*29,subr(13)*6

    COMMON /xmssg / ufm,uwm,uim
    COMMON /ifpx1 / ncds,t1(2,370)
    COMMON /machin/ mach,dummy4(4)
    COMMON /system/ system,outtap,nogo,intap,dumm15(15),plotf,  &
        dumm6(6),axic,dummy3(3),hicore,dummy6(6),  &
        axif,dumm30(30),isubs,isy70(7),isy77
    COMMON /xechox/ echo(4)
    COMMON /xxread/ inflag,insave

    DATA     bcd1 , bcd2,  bcd3,  bcd4  ,bcd5  ,bcd6,  bcd7   /  &
             4HXCSA,4HIFP1,4HXSOR,4HXGPI,4HGNFI,4HTTIO,4HTTLP /
    DATA     bcd8 , bcd9,  bcd10 ,bcd11                       /  &
            4HTTOT,4HSOLI,4HFLUI,1HD                         /
    DATA     subr / 'NASCAR','GNFIAT','TMTSIO','TMTSLP','XCSA  ',  &
           'TMTSOT','ASDMAP','IFP1  ','XSORT2','IFP   ', 'IFP3  ', &
           'XGPI  ','BANDIT'/

    insave = intap
 
    !     READ AND PROCESS THE NASTRAN CARD (IF PRESENT).
 
    !hgs disable nascar
    iby = 0
    IF(iby == 0) THEN
        IF (debug1 > 0) WRITE (6,10) subr(1)
10      FORMAT (/,' -LINK1 DEBUG- SEMINT CALLING ',a6,' NEXT',/)
        CALL nascar
    END IF
 
    !     DEFINE OPEN CORE FOR UNIVAC, VAX, AND UNIX
 
    IF (mach == 3 .OR. mach >= 5) CALL defcor
 
    !     GENERATE INITIAL FILE TABLES.
    !     COMPUTE NASTRAN TIMING CONSTANTS.
    !     READ EXECUTIVE CONTROL DECK AND SAVE NOGO FLAG.
    !     READ CASE CONTROL DECK, SORT BULK DATA AND EXECUTE
    !     INPUT FILE PROCESSOR UNLESS BULK DATA IS MISSING.
    !     IF CONICAL SHELL PROBLEM, EXECUTE IFP3.
 
    CALL conmsg (bcd5,1,1)
    IF (debug1 > 0) WRITE (6,10) subr(2)
    CALL gnfiat
 
    !     CALL THE TIME TEST ROUTINES TO COMPUTE THE NASTRAN
    !     TIMING CONSTANTS AND INITIALIZE COMMON /NTIME/
 
    !     GENERATE THE I/O TIMES AND
    !     CPU TIMES FOR VARIOUS TYPES OF LOOPS
 
    CALL conmsg (bcd6,1,0)
    IF (debug1 > 0) WRITE (6,10) subr(3)
    CALL tmtsio (*2000,debug1)
    CALL conmsg (bcd7,1,0)
    IF (debug1 > 0) WRITE (6,10) subr(4)
    CALL tmtslp
 
    !     PROCESS EXECUTIVE CONTROL CARDS
 
2000 CALL conmsg (bcd1,1,1)
    IF (debug1 > 0) WRITE (6,10) subr(5)
    CALL xcsa
 
    !     OUTPUT THE COMMON /NTIME/ ENTRIES IF DIAG 35 IS TURNED ON
 
    CALL sswtch (35,l35)
    IF (l35 == 0) GO TO 3000
    CALL conmsg (bcd8,1,0)
    IF (debug1 > 0) WRITE (6,10) subr(6)
    CALL tmtsot
 
    !     PROCESS SUBSTRUCTURING DMAP
 
3000 nogox = nogo
    nogo  = 0
    IF (debug1 > 0 .AND. isubs /= 0) WRITE (6,10) subr(7)
    IF (isubs /= 0) CALL asdmap
 
    !     PROCESS CASE CONTROL CARDS
 
    CALL conmsg (bcd2,1,1)
    IF (debug1 > 0) WRITE (6,10) subr(8)
    CALL ifp1
    nogo1 = nogo
    IF (nogo == -9) nogo = 1
    IF (nogo <  0) nogo = 0
    kaxif = 0
 
    !     REVERT TO OLD XSORT TO PROCESS BULKDATA CARDS IF DIAG 42 IS
    !     TURNED ON,  OTHERWISE, USE XSORT2 FOR SPEED AND EFFICIENCY
 
    CALL conmsg (bcd3,1,0)
    CALL sswtch (42,l42)
    IF (debug1 > 0) WRITE (6,10) subr(9)
    IF (l42  ==  1) CALL xsort
    IF (l42  ==  0) CALL xsort2
    IF (nogo == -2) GO TO 4000
 
    !     INPUT FILE PROCESSOR(S) TO CHECK EACH BULKDATA CARD
 
    IF (debug1 > 0) WRITE (6,10) subr(10)
    CALL ifp
    IF (debug1 > 0 .AND. axic /= 0) WRITE (6,10) subr(11)
    IF (axic /= 0) CALL ifp3
 
    !     SET KAXIF AS IFP4 WILL MODIFY AXIF
 
    kaxif = axif
    IF (kaxif == 1 .OR. kaxif == 3) CALL ifp4
    IF (kaxif == 2 .OR. kaxif == 3) CALL ifp5
 
    !     SUPPRESS NOGO FLAG IF USER REQUESTS UNDEFORMED STRUCTURE PLOT VIA
    !     NASTRAN PLOTOPT CARD
 
4000 IF (nogo == -2) nogo = 0
    IF (nogo == 0 .AND. nogo1 < 0) nogo = nogo1
    IF (nogo >= 1 .AND. nogo1 < 0) nogo = -9
    IF (nogo1 == 0) nogo1 = nogo
 
    !     NOGO FLAG CONDITIONS
    !     NOGOX.NE. 0, FATAL ERROR IN EXECUTIVE CONTROL
    !     NOGO .EQ.-9, FATAL ERROR IN BULKDATA AND IN PLOT COMMANDS
    !     NOGO .EQ. 0, NO FATAL ERROR DETECTED IN ENTIRE INPOUT DECK
    !     NOGO .GT. 0, FATAL ERROR IN BULKDATA, NO ERROR IN PLOT COMMANDS
    !     NOGO .LT. 0, NO ERROR IN BULKDATA, FATAL ERROR IN PLOT COMMANDS
 
    IF (nogox /= 0) GO TO 5500
    IF (nogo < 0) THEN
        GO TO  4100
    ELSE IF (nogo == 0) THEN
        GO TO  4300
    ELSE
        GO TO  4200
    END IF
4100 IF (nogo == -9 .AND. plotf /= 3) GO TO 5500
    IF (plotf <= 1) GO TO 4200
    nogo = 0
    GO TO 4300
4200 nogo = 1
 
    !     EXECUTE GENERAL PROBLEM INITIALIZATION IF DATA PERMITS.
 
4300 IF (nogo /= 0) CALL mesage (-61,0,0)
    CALL conmsg (bcd4,1,0)
    IF (debug1 > 0) WRITE (6,10) subr(12)
    CALL xgpi
 
    !     CALL BANDIT TO GENERATE GRID-POINT RE-SEQUENCE CARDS IF DATA
    !     PERMITS
 
    IF (nogo /= 0 .AND. nogo1 < 0) nogo = -9
    IF (nogo == 0 .AND. nogo1 /= 0) nogo = nogo1
    IF (isy77 < 0 .OR. nogo /= 0) GO TO 5100
    IF (axic /= 0  .OR. kaxif == 1 .OR. kaxif == 3) GO TO 5000
    IF (debug1 >  0) WRITE (6,10) subr(13)
    CALL bandit
    GO TO 5100
5000 WRITE (outtap,6100) uim
    bcdx = bcd10
    IF (axic /= 0) bcdx = bcd9
    WRITE (outtap,6200) bcdx,bcd11
    WRITE (outtap,6300)
 
    !     TERMINATE NASTRAN IF LINK 1 ONLY IS REQUESTED BY USER
 
5100 IF (isy77 == -2) CALL pexit
 
    !     EXIT ACCORDING TO PLOT OPTION REQUEST
    !     SET PLOTF TO NEGATIVE ONLY IF JOB IS TO BE TERMINATED AFTER PLOTS
    !     IN LINK2
 
    j = plotf + 1
    IF (nogo == 0) THEN
        SELECT CASE ( j )
            CASE (    1)
                GO TO 5800
            CASE (    2)
                GO TO 5800
            CASE (    3)
                GO TO 5700
            CASE (    4)
                GO TO 5700
            CASE (    5)
                GO TO 5800
            CASE (    6)
                GO TO 5800
        END SELECT
    END IF
    IF (nogo > 0) THEN
        SELECT CASE ( j )
            CASE (    1)
                GO TO 5300
            CASE (    2)
                GO TO 5300
            CASE (    3)
                GO TO 5600
            CASE (    4)
                GO TO 5600
            CASE (    5)
                GO TO 5600
            CASE (    6)
                GO TO 5600
        END SELECT
    END IF
    IF (nogo < 0) THEN
        SELECT CASE ( j )
            CASE (    1)
                GO TO 5500
            CASE (    2)
                GO TO 5500
            CASE (    3)
                GO TO 5500
            CASE (    4)
                GO TO 5600
            CASE (    5)
                GO TO 5600
            CASE (    6)
                GO TO 5200
        END SELECT
    END IF
    !                      PLOTF =   0,   1,   2,   3,   4,   5
 
5200 IF (nogo+9 == 0) THEN
        GO TO  5500
    ELSE
        GO TO  5800
    END IF
5300 IF (plotf > 1) WRITE (outtap,5400)
5400 FORMAT ('0*** ATTEMPT TO PLOT UNDEFORMED MODEL IS ABANDONED DUE',  &
        ' TO FATAL ERROR IN BULK DATA')
5500 CALL mesage (-61,0,0)
5600 WRITE  (outtap,5650) uwm
5650 FORMAT (a25,' - FATAL ERRORS ENCOUNTERED IN USER INPUT DECK,',  &
        /5X,'HOWEVER, NASTRAN WILL ATTEMPT TO PLOT THE UNDEFORMED',  &
        ' STRUCTURE AS REQUESTED BY USER')
5700 plotf = -plotf
5800 RETURN
 
6100 FORMAT (a29,' - GRID-POINT RESEQUENCING PROCESSOR BANDIT IS NOT',  &
        ' USED DUE TO')
6200 FORMAT (5X,'THE PRESENCE OF AXISYMMETRIC ',a4,a1,' DATA')
6300 FORMAT (1H0,10X,'**NO ERRORS FOUND - EXECUTE NASTRAN PROGRAM**')
 
END SUBROUTINE semint
