SUBROUTINE shcsgd (*,cflag,ccsid,ctheta,pflag,pcsid,ptheta,  &
        necpt,tubd,csid,thetad,tumsd)
     
!     WITH ENTRY SHCSGS (*,CFLAG,CCSID,CTHETA,PFLAG,PCSID,PTHETA,
!    1                   NECPT,TUBS,CSID,THETAS,TUMSS)
 
 
!     'COORDINATE SYSTEM GENERATOR' ROUTINE FOR SHELL ELEMENTS.
 
!     THIS ROUTINE USES THE VALUES IN THE EST TABLE TO CREATE
!     APPROPRIATE MATERIAL/STRESS COORDINATE SYSTEM TRANSFORMATIONS.
 
!     INPUT:
!            CFLAG    - INDICATOR FLAG FROM CONNECTION
!            CCSID    - CSID  FROM CONNECTION
!            CTHETA   - ANGLE FROM CONNECTION
!            PFLAG    - INDICATOR FLAG FROM PROPERTY
!            PCSID    - CSID  FROM PROPERTY
!            PTHETA   - ANGLE FROM PROPERTY
!            NECPT    - ARRAY OF LENGTH 4, WORDS 2-4 ARE THE LOCATION
!                       WHERE THE TRANSFORMATION NEEDS TO BE CALCULATED
!            TUBD/S   - USER TO BASIC TRANSFORMATION
!     OUTPUT:
!            TUMSD/S  - USER TO MATERIAL/STRESS TRANSFORMATION
!            CSID     - CSID  USED FOR CALCULATIONS
!            THETAD/S - THETA USED FOR CALCULATIONS
 
!     NOTES:
!     1- IF CSID HAS BEEN SPECIFIED, SUBROUTINE TRANSD IS CALLED TO
!        CALCULATE [TBMS] (MATERIAL/STRESS TO BASIC TRANSFORMATION).
!        [TBMS] IS THEN PREMULTIPLIED BY [TUB] TO OBTAIN [TUMS].
!        THEN USING THE PROJECTION OF X-AXIS, AN ANGLE IS CALCULATED
!        UPON WHICH STEP 2 IS TAKEN.
!     2- IF THETA HAS BEEN SPECIFIED, INPLANE TRANSFORMATION IS USED TO
!        CALCULATE [TUMS] (MATERIAL/STRESS TO USER TRANSFORMATION).
!     3- IF THE CONNECTION VALUE IS LEFT BLANK, THE PROPERTY VALUE IS
!        USED.
!     4- NON-STANDARD RETURN IS TAKEN WHEN THE X-AXIS OF THE SPECIFIED
!        COORDINATE SYSTEM DOES NOT HAVE A PROJECTION ON THE X-Y PLANE
!        OF THE ELEMENT COORD. SYSTEM
 
 
 
 , INTENT(IN OUT)                         :: *
 INTEGER, INTENT(IN OUT)                  :: cflag
 INTEGER, INTENT(IN)                      :: ccsid
 REAL, INTENT(IN)                         :: ctheta
 INTEGER, INTENT(IN OUT)                  :: pflag
 INTEGER, INTENT(IN)                      :: pcsid
 REAL, INTENT(IN)                         :: ptheta
 INTEGER, INTENT(OUT)                     :: necpt(4)
 DOUBLE PRECISION, INTENT(IN)             :: tubd(9)
 INTEGER, INTENT(OUT)                     :: csid
 DOUBLE PRECISION, INTENT(OUT)            :: thetad
 DOUBLE PRECISION, INTENT(OUT)            :: tumsd(9)
 
 REAL :: tubs(9),tumss(9),tbmss(9),xms,yms,thetas,eps1s,  &
     pis,twopis,raddgs,degrds,flips
 DOUBLE PRECISION :: tbmsd(9),xmd,ymd, eps1d, pid,twopid,raddgd,degrdd,flipd
 COMMON /condas/  pis,twopis,raddgs,degrds
 COMMON /condad/  pid,twopid,raddgd,degrdd
 EQUIVALENCE      (tbmss(1),tbmsd(1))
 DATA    eps1d ,  eps1s /1.0D-7, 1.0E-7 /
 
 
!     DOUBLE PRECISION VERSION
 
 flipd = 1.0D0
 IF (cflag == 0) GO TO 130
 
!     DETERMINE THETA FROM THE PROJECTION OF THE X-AXIS OF THE MATERIAL/
!     STRESS COORD. SYSTEM, DETERMINED BASED ON CCSID, ONTO THE XY-PLANE
!     OF THE ELEMENT COORD. SYSTEM.
 
 csid = ccsid
 IF (ccsid > 0) GO TO 110
 
!     [TUMS] = [TUB]
 
 DO  i = 1,9
   tumsd(i) = tubd(i)
 END DO
 GO TO 120
 
!     [TUMS] = [TUB] [TBMS]
 
 110  necpt(1) = ccsid
 CALL transd (necpt,tbmsd)
 CALL gmmatd (tubd,3,3,0, tbmsd,3,3,0, tumsd)
 
 120  xmd = tumsd(1)
 ymd = tumsd(4)
 IF (DABS(xmd) <= eps1d .AND. DABS(ymd) <= eps1d) RETURN 1
 thetad = DATAN2(ymd,xmd)
 IF (tumsd(9) < 0.0D0) flipd = -1.0D0
 GO TO 190
 
 130  IF (ctheta == 0.0) GO TO 140
 
!     DETERMINE THETA FROM CTHETA
 
 thetad = DBLE(ctheta)*degrdd
 GO TO 190
 
!     DEFAULT IS CHOSEN, LOOK FOR VALUES OF PCSID AND/OR PTHETA ON THE
!     PSHELL CARD.
 
 140  IF (pflag == 0) GO TO 180
 
!     DETERMINE THETA FROM THE PROJECTION OF THE X-AXIS OF THE MATERIAL/
!     STRESS COORD. SYSTEM, DETERMINED BASED ON PCSID, ONTO THE XY-PLANE
!     OF THE ELEMENT COORD. SYSTEM.
 
 csid = pcsid
 IF (pcsid > 0) GO TO 160
 
!     [TUMS] = [TUB]
 
 DO  i = 1,9
   tumsd(i) = tubd(i)
 END DO
 GO TO 170
 
!     [TUMS] = [TUB] [TBMS]
 
 160  necpt(1) = pcsid
 CALL transd (necpt,tbmsd)
 CALL gmmatd (tubd,3,3,0, tbmsd,3,3,0, tumsd)
 
 170  xmd = tumsd(1)
 ymd = tumsd(4)
 IF (DABS(xmd) <= eps1d .AND. DABS(ymd) <= eps1d) RETURN 1
 thetad = DATAN2(ymd,xmd)
 IF (tumsd(9) < 0.0D0) flipd = -1.0D0
 GO TO 190
 
!     DETERMINE THETA FROM PTHETA
 
 180  thetad = DBLE(ptheta)*degrdd
 
!     IF THE Z-AXIS OF THE TARGET MATERIAL/STRESS COORD. SYSTEM WAS NOT
!     POINTING IN THE SAME GENERAL DIRECTION AS THE Z-AXIS OF THE USER
!     COORD. SYSTEM, FLIP THE Y- AND Z-AXES OF THE FINAL COORDINATE
!     SYSTEM TO ACCOUNT FOR IT.
 
 190  tumsd(1) = DCOS(thetad)
 tumsd(2) =-flipd*DSIN(thetad)
 tumsd(3) = 0.0D0
 tumsd(4) = DSIN(thetad)
 tumsd(5) = flipd*DCOS(thetad)
 tumsd(6) = 0.0D0
 tumsd(7) = 0.0D0
 tumsd(8) = 0.0D0
 tumsd(9) = flipd
 
 RETURN
 
 
 ENTRY shcsgs (*,cflag,ccsid,ctheta,pflag,pcsid,ptheta,  &
     necpt,tubs,csid,thetas,tumss)
!     ======================================================
 
!     SINGLE PRECISION VERSION
 
 flips = 1.0
 IF (cflag == 0) GO TO 230
 
!     DETERMINE THETA FROM THE PROJECTION OF THE X-AXIS OF THE MATERIAL/
!     STRESS COORD. SYSTEM, DETERMINED BASED ON CCSID, ONTO THE XY-PLANE
!     OF THE ELEMENT COORD. SYSTEM.
 
 csid = ccsid
 IF (ccsid > 0) GO TO 210
 
!     [TUMS] = [TUB]
 
 DO  i = 1,9
   tumss(i) = tubs(i)
 END DO
 GO TO 220
 
!     [TUMS] = [TUB] [TBMS]
 
 210  necpt(1) = ccsid
 CALL transs (necpt,tbmss)
 CALL gmmats (tubs,3,3,0, tbmss,3,3,0, tumss)
 
 220  xms = tumss(1)
 yms = tumss(4)
 IF (ABS(xms) <= eps1s .AND. ABS(yms) <= eps1s) RETURN 1
 thetas = ATAN2(yms,xms)
 IF (tumss(9) < 0.0) flips = -1.0
 GO TO 290
 
 230  IF (ctheta == 0.0) GO TO 240
 
!     DETERMINE THETA FROM CTHETA
 
 thetas = ctheta*degrds
 GO TO 290
 
!     DEFAULT IS CHOSEN, LOOK FOR VALUES OF PCSID AND/OR PTHETA ON THE
!     PSHELL CARD.
 
 240  IF (pflag == 0) GO TO 280
 
!     DETERMINE THETA FROM THE PROJECTION OF THE X-AXIS OF THE MATERIAL/
!     STRESS COORD. SYSTEM, DETERMINED BASED ON PCSID, ONTO THE XY-PLANE
!     OF THE ELEMENT COORD. SYSTEM.
 
 csid = pcsid
 IF (pcsid > 0) GO TO 260
 
!     [TUMS] = [TUB]
 
 DO  i = 1,9
   tumss(i) = tubs(i)
 END DO
 GO TO 270
 
!     [TUMS] = [TUB] [TBMS]
 
 260  necpt(1) = pcsid
 CALL transs (necpt,tbmss)
 CALL gmmats (tubs,3,3,0, tbmss,3,3,0, tumss)
 
 270  xms = tumss(1)
 yms = tumss(4)
 IF (ABS(xms) <= eps1s .AND. ABS(yms) <= eps1s) RETURN 1
 thetas = ATAN2(yms,xms)
 IF (tumss(9) < 0.0) flips = -1.0
 GO TO 290
 
!     DETERMINE THETA FROM PTHETA
 
 280  thetas = ptheta*degrds
 
!     IF THE Z-AXIS OF THE TARGET MATERIAL/STRESS COORD. SYSTEM WAS NOT
!     POINTING IN THE SAME GENERAL DIRECTION AS THE Z-AXIS OF THE USER
!     COORD. SYSTEM, FLIP THE Y- AND Z-AXES OF THE FINAL COORDINATE
!     SYSTEM TO ACCOUNT FOR IT.
 
 290  tumss(1) = COS(thetas)
 tumss(2) =-flips*SIN(thetas)
 tumss(3) = 0.0
 tumss(4) = SIN(thetas)
 tumss(5) = flips*COS(thetas)
 tumss(6) = 0.0
 tumss(7) = 0.0
 tumss(8) = 0.0
 tumss(9) = flips
 
 RETURN
END SUBROUTINE shcsgd
