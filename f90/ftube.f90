SUBROUTINE ftube
     
!     THIS IS THE FLUID TUBE ELEMENT IN HEAT TRANSFER.
!     IT COMPUTES AND OUTPUTS THE CONDUCTIVITY AND/OR CAPACITY MATRICES
!     OF THE ELEMENT.
 
!     - SINGLE AND DOUBLE PRECISION VERSION -
 
!     EST ENTRY FOR -FTUBE- ELEMENT.
!     ==============================
 
!     EST( 1) = ELEMENT ID
!     EST( 2) = SIL-A
!     EST( 3) = SIL-B
!     EST( 4) = HEAT CAPACITY/UNIT VOLUME = RHO C
!     EST( 5) = VOLUME FLOW RATE = VDOT         P
!     EST( 6) = DIAMETER AT A
!     EST( 7) = DIAMETER AT B = DIAMETER AT A IF NOT DEFINED.
!     EST( 8) = CSID-A  NOT USED
!     EST( 9) = XA
!     EST(10) = YA
!     EST(11) = ZA
!     EST(12) = CSID-B  NOT USED
!     EST(13) = XB
!     EST(14) = YB
!     EST(15) = ZB
!     EST(16) = AVG TEMP OF ELEMENT.  NOT USED.
 
 
 LOGICAL :: heat     ,error
 INTEGER :: dict(7)  ,estid    ,iest(1)
 REAL :: rk(4)    ,id1      ,id2      ,dict5
 DOUBLE PRECISION :: k(4)     ,length
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg /  ufm      ,uwm      ,uim
 COMMON /system/  sysbuf   ,ioutpt
 COMMON /emgprm/  dum15(15),kmb(3)   ,iprec    ,error    ,heat
 COMMON /emgdic/  dmmm(4)  ,estid
 COMMON /emgest/  est(16)
 COMMON /condas/  pi
 EQUIVALENCE      (iest(1) , est(1)) ,(rk(1)   , k(1))   , (dict(5) , dict5 )
 
 IF (.NOT.heat) GO TO 240
 dict(1) = estid
 dict(2) = 1
 dict(3) = 2
 dict(4) = 1
 dict5   = 0.0
 IF (kmb(1) == 0) GO TO 170
 
!     CONDUCTIVITY
 
 rhocp = est(4)
 vdot  = est(5)
 
!     STORE CONDUCTIVITY BY COLUMNS
 
 k(1) = DBLE(rhocp*vdot)
 
 k(2) = -k(1)
 k(3) = 0.0D0
 k(4) = 0.0D0
 
!     OUTPUT VIA EMGOUT THE FULL MATRIX IN GLOBAL, UNSYMETRIC
 
 ifil = 1
 isze = 4
 ASSIGN 170 TO irtn
 IF (iprec == 2) GO TO 160
 150 rk(1) = SNGL(k(1))
 rk(2) = SNGL(k(2))
 rk(3) = SNGL(k(3))
 rk(4) = SNGL(k(4))
 160 CALL emgout (rk(1),k(1),isze,1,dict,ifil,iprec)
 GO TO irtn, (170,240)
 
!     CAPACITY MATRIX
 
 170 IF (kmb(3) == 0) GO TO 240
 rhocp = est( 4)
 vdot  = est( 5)
 id1   = est( 6)
 IF (est(7) == 0.0) THEN
   GO TO   180
 ELSE
   GO TO   190
 END IF
 180 id2 = id1
 GO TO 200
 190 id2   = est( 7)
 200 xa    = est( 9)
 ya    = est(10)
 za    = est(11)
 xb    = est(13)
 yb    = est(14)
 zb    = est(15)
 length = DBLE((xb-xa))**2 + DBLE((yb-ya))**2 + DBLE((zb-za))**2
 IF (length > 0.0D0) GO TO 220
 length = DSQRT(length)
 WRITE  (ioutpt,210) uim,iest(1)
 210 FORMAT (a29,' FROM ELEMENT FTUBE -', /5X,'ELEMENT WITH ID =',i9,  &
     ' HAS A ZERO LENGTH.')
 error = .true.
 
!     FILL AND OUTPUT CAPACITY MATRIX BY COLUMNS IN GLOBAL, SYMMETRIC.
 
 220 k(1) = (DBLE(rhocp*pi*(id1+id2)))**2*length/32.0D0
 k(2) = 0.0D0
 k(3) = 0.0D0
 k(4) = k(1)
 dict(2) = 2
 ifil = 3
 isze = 2
 ASSIGN 240 TO irtn
 IF (iprec-1 < 0) THEN
   GO TO   240
 ELSE IF (iprec-1 == 0) THEN
   GO TO   150
 ELSE
   GO TO   160
 END IF
 
 240 RETURN
END SUBROUTINE ftube
