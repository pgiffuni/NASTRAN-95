SUBROUTINE shstns (numpx,elid,igrid,z12,epslni,bendng,idr)
     
!     TO CALCULATE SHELL ELEMENT STRAINS FOR A 2-D FORMULATION BASE.
!     COMPOSITE LAYER STRAINS ARE NOT CALCULATED IN THIS ROUTINE.
 
 
!     INPUT :
!           NUMPX  - NUMBER OF EVALUATION POINTS
!           ELID   - ELEMENT ID
!           IGRID  - ARRAY IF EXTERNAL GRID IDS
!           Z12    - EVALUATION POINT FIBER DISTANCES
!           EPSLNI - CORRECTED STRAINS AT EVALUATION POINTS
!           BENDNG - INDICATES THE PRESENCE OF BENDING BEHAVIOR
!           IDR    - REORDERING ARRAY BASED ON EXTERNAL GRID POINT ID'S
!          /OUTREQ/- OUTPUT REQUEST LOGICAL FLAGS
 
!     OUTPUT:
!           STRAINS ARE PLACED AT THE PROPER LOCATION IN /SDR2X7/.
 
 
!     THE STRAIN OUTPUT DATA BLOCK, UAI CODE
 
!     ADDRESS    DESCRIPTIONS
 
!        1       ELID
!     --------------------------------------------------------------
!        2       GRID POINT NUMBER OR 'CNTR'
!      3 - 10    STRAINS FOR LOWER POINTS OR MEMBRANE STRAINS
!     11 - 18    STRAINS FOR UPPER POINTS OR BENDING CURVATURES
!     ---------- ABOVE DATA REPEATED 3 TIMES
!                FOR GRID POINTS
 
 
!     THE STRAIN OUTPUT DATA BLOCK, AT ELEMENT CENTER ONLY, COSMIC
 
!     ADDRESS    DESCRIPTIONS
 
!        1       ELID
!     --------------------------------------------------------------
!        2       LOWER FIBER DISTANCE
!      3 -  9    STRAINS FOR LOWER POINTS OR MEMBRANE STRAINS
!       10       UPPER FIBER DISTANCE
!     11 - 17    STRAINS FOR UPPER POINTS OR BENDING CURVATURES
!     ---------- ABOVE DATA REPEATED 3 TIMES
!                FOR GRID POINTS
 
 
 
 INTEGER, INTENT(IN)                      :: numpx
 INTEGER, INTENT(IN)                      :: elid
 INTEGER, INTENT(IN)                      :: igrid(1)
 REAL, INTENT(IN)                         :: z12(2,1)
 REAL, INTENT(IN)                         :: epslni(6,1)
 LOGICAL, INTENT(IN OUT)                  :: bendng
 INTEGER, INTENT(IN OUT)                  :: idr(1)
 LOGICAL :: stsreq,stnreq,forreq,strcur,  &
     grids,vonms,layer,gridss,vonmss,layers,cosmic
 INTEGER :: nstrin(1)
!WKBI NCL93012 3/94
 INTEGER :: nstres(1)
 REAL :: epsil(3),epss,fiber,epsilp(4)
 COMMON /sdr2x7/ dum71(100),stres(100),forsul(200),strin(100)
 COMMON /outreq/ stsreq,stnreq,forreq,strcur,grids,vonms,layer  &
     ,               gridss,vonmss,layers
 EQUIVALENCE     (nstrin(1),strin(1))
!WKBI NCL93012 3/94
 EQUIVALENCE     (nstres(1), stres(1))
!WKBNB 7/94 SPR94004
 LOGICAL :: ostrai
 COMMON / BLANK/ app(2), sort2, idum(2), comps, skp(4), ostrai  &
     ,               sk2(39), midve
!WKBNE 7/94 SPR94004
 DATA    cosmic, epss / .true., 1.0E-17 /
 
 
!     ELEMENT ENTER COMPUATION ONLY FOR COSMIC
!     I.E. CALLER SHOULD PASS 1 IN NUMPX FOR COSMIC, 4 FOR UAI
 
 nump = numpx
 IF (cosmic) nump = 1
 
 nstrin(1) = elid
 
!     START THE LOOP ON EVALUATION POINTS
 
 nump1 = nump - 1
 DO  inplan = 1,nump
   
   istrin = 1
   IF (cosmic) GO TO 140
   
   istrin = (inplan-1)*17 + 2
   nstrin(istrin) = inplan - 1
   IF (.NOT.gridss .OR. inplan <= 1) GO TO 130
   DO  inptmp = 1,nump1
     IF (idr(inptmp) == igrid(inplan)) GO TO 120
   END DO
   CALL errtrc ('SHSTNS  ',100)
   120 istrin = inptmp*17 + 2
   nstrin(istrin) = igrid(inplan)
   130 IF (inplan == 1) nstrin(istrin) = igrid(inplan)
   
!     START THE LOOP ON FIBERS
   
   140 DO  iz = 1,2
     IF (.NOT.strcur) GO TO 190
     
!     IF STRAIN/CURVATURE IS REQUESTED, SIMPLY OUTPUT THE AVAILABLE
!     STRAINS.
     
     strin(istrin+1) = 0.0
     DO  i = 1,3
       epsil(i) = 0.0
     END DO
!WKBI 7/94 SPR94004
     IF ( ostrai .AND. iz == 2 ) GO TO 171
     IF (iz /= 1) GO TO 170
     DO  i = 1,3
       epsil(i) = epslni(i,inplan)
     END DO
!WKBI 7/94 SPR94004
     IF ( ostrai .AND. iz == 1 ) GO TO 220
     170 IF (.NOT.bendng .OR. iz /= 2) GO TO 190
!WKBI 7/94 SPR94004
     171 CONTINUE
     DO  i = 1,3
       epsil(i) = epslni(i+3,inplan)
     END DO
     GO TO 220
     
!     IF FIBER STRAINS ARE REQUESTED, EVALUATE STRAINS AT THIS FIBER
!     DISTANCE
     
     190 fiber = z12(iz,inplan)
     strin(istrin+1) = fiber
     DO  i = 1,3
       epsil(i) = epslni(i,inplan) - epslni(i+3,inplan)*fiber
     END DO
     
!     CLEANUP AND SHIP CALCULATED STRAINS
     
     220 DO  its = 1,3
       IF (ABS(epsil(its)) <= epss) epsil(its) = 0.0
!WKBR NCL93012 3/94      STRIN(ISTRIN+1+ITS) = EPSIL(ITS)
       stres(istrin+1+its) = epsil(its)
     END DO
     
!     CALCULATE PRINCIPAL STRAINS
     
     CALL shpsts (epsil(1),vonmss,epsilp)
!WKBDB NCL93012 3/94
!      STRIN(ISTRIN+5) = EPSILP(1)
!      STRIN(ISTRIN+6) = EPSILP(2)
!      STRIN(ISTRIN+7) = EPSILP(3)
!      STRIN(ISTRIN+8) = EPSILP(4)
!WKBDE NCL93012 3/94
!WKBNB NCL93012 3/94
     nstres( istrin+1 ) = 0
     IF ( iz == 2 ) nstres( istrin+1) = -1
     stres(istrin+5) = epsilp(1)
     stres(istrin+6) = epsilp(2)
     stres(istrin+7) = epsilp(3)
     stres(istrin+8) = epsilp(4) * 2.
!WKBNE NCL93012 3/94
     
     istrin = istrin + 8
   END DO
 END DO
 
 RETURN
END SUBROUTINE shstns
