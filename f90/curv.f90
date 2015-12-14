SUBROUTINE curv
     
!     MAIN DRIVING ROUTINE OF MODULE -CURV-.
 
!     DMAP CALLING SEQUENCE.
 
!     CURV   OES1,MPT,CSTM,EST,SIL,GPL/OES1M,OES1G/P1/P2 $
 
 LOGICAL :: foes1g, eofos1, strain
 INTEGER :: subr(6), FILE, mcb(7)
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm, uwm, uim, sfm
 COMMON /BLANK / ip1, ip2
 COMMON /system/ isysbf, ioutpt
 COMMON /curvtb/ indexs(108)
 COMMON /zzzzzz/ iz(1)
 
!     COMMON /ZZCURV/ MUST BE AT THE LONGEST OF OVERLAYS WITH CURV1,
!     CURV2, AND CURV3.
 
 EQUIVALENCE     (indexs( 16),lmcsid), (indexs( 52),lcore),  &
     (indexs( 79),loc), (indexs( 80),FILE),  &
     (indexs( 81),imsg), (indexs(100),eofos1),  &
     (indexs(103),foes1g), (indexs(104),strain), (indexs(105),logerr)
 DATA     subr / 4HCURV,4H1   ,4HCURV,4H2   ,4HCURV,4H3     /
 
 
!     CHECK TO SEE IF COMPUTATIONS NEED TO BE DONE
 
 IF (ip1 < 0) RETURN
 
!     CHECK TO SEE IF THE INPUT FILE EXISTS
 
 mcb(1) = 101
 CALL rdtrl (mcb(1))
 IF (mcb(1) <= 0) RETURN
 
!     PERFORM INITIALIZATION AND CREATE ESTX ON SCRATCH FILE 1.
 
 DO  i = 1,107
   indexs(i) = 777777777
 END DO
 imsg = 0
 jsub = 1
 CALL curv1
 IF (imsg  == -8) GO TO 10001
 IF (imsg  <  0) GO TO 9000
 IF (lmcsid <= 0) GO TO 8000
 
!     CREATE OES1M FOR NEXT SUBCASE IF NOT AT EOF IN OES1.
 
 100 IF (eofos1) GO TO 4000
 jsub = 2
 CALL curv2
 IF (imsg == -8) GO TO 10001
 IF (imsg <  0) GO TO 9000
 
!     IF OES1G IS TO BE FORMED CALL CURV3 OVERLAY.  PROCESS CURRENT
!     SUBCASE
 
 IF (.NOT.foes1g) GO TO 100
 jsub = 3
 CALL curv3
 IF (imsg == -8) GO TO 10001
 IF (imsg <  0) GO TO 9000
 GO TO 100
 
!     EOF HIT IN OES1.  ALL THROUGH.
 
 4000 CONTINUE
 RETURN
 
!     NO NON-ZERO MATERIAL COORDINATE SYSTEM IDS ENCOUNTERED
 
 8000 CALL page2 (3)
 WRITE (ioutpt,8100) uwm
 IF (.NOT.strain) WRITE (ioutpt,8200)
 IF (     strain) WRITE (ioutpt,8300)
 8100 FORMAT (a25,' 3173, NO NON-ZERO MATERIAL COORDINATE SYSTEM IDS ',  &
     'ENCOUNTERED IN MODULE CURV.')
 8200 FORMAT (39H stresses in material coordinate system, 14H NOT computed.)
 8300 FORMAT (49H strains/curvatures in material coordinate system,  &
     14H NOT computed.)
 GO TO 4000
 
!     ERROR CONDITION IN CURV1, CURV2, OR CURV3.
 
 9000 IF (imsg /= -37) GO TO 9999
 WRITE  (ioutpt,9100) sfm,jsub,imsg,loc,jsub,FILE
 9100 FORMAT (a25,' 3174, SUBROUTINE CURV',i1,  &
     ' HAS RETURNED WITH ERROR CONDITION ',i4, /5X,  &
     'LOCATION CODE = ',i4,' IN SUBROUTINE CURV',i1, /5X, 'FILE NUMBER   = ',i4)
 WRITE  (ioutpt,9998) indexs
 9998 FORMAT (/5X,29H constants in COMMON /curvtb/ , /,(3X,4I15))
 
!     INSURE ALL FILES CLOSED
 
 9999 CONTINUE
 DO  i = 1,9
   DO  j = 100,300,100
     CALL CLOSE (i+j,1)
   END DO
 END DO
 10001 WRITE (ioutpt,9100) sfm,jsub,imsg,loc,jsub,FILE
 jsub = 2*jsub - 1
 IF (imsg == -8) FILE = lcore
 CALL mesage (imsg,FILE,subr(jsub))
 GO TO 4000
END SUBROUTINE curv
