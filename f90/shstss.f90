SUBROUTINE shstss (numpx,elid,igrid,thikns,z12,g,epscsi,stemp,  &
        tbar,g2alfb,bendng,idr)
     
!     TO CALCULATE SHELL ELEMENT STRESSES FOR A 2-D FORMULATION BASE.
!     COMPOSITE LAYER STRESSES ARE NOT CALCULATED IN THIS ROUTINE.
 
 
!     INPUT :
!           NUMPX  - NUMBER OF EVALUATION POINTS
!           ELID   - ELEMENT ID
!           IGRID  - ARRAY IF EXTERNAL GRID IDS
!           THIKNS - EVALUATION POINT THICKNESSES
!           Z12    - EVALUATION POINT FIBER DISTANCES
!           G      - 6X6 STRESS-STRAIN MATRIX
!           EPSCSI - CORRECTED STRAINS AT EVALUATION POINTS
!           STEMP  - TEMPERATURE DATA FOR STRESS RECOVERY
!           TBAR   - AVERAGE ELEMENT TEMPERATURE
!           G2ALFB - MATRIX USED IN RECORRECTING OF STRESSES
!           BENDNG - INDICATES THE PRESENCE OF BENDING BEHAVIOR
!           IDR    - REORDERING ARRAY BASED ON EXTERNAL GRID POINT ID'S
!          /TMPDAT/- TEMPERATURE-RELATED LOGICAL FLAGS
!          /OUTREQ/- OUTPUT REQUEST LOGICAL FLAGS
 
!     OUTPUT:
!           STRESSES ARE PLACED AT THE PROPER LOCATION IN /SDR2X7/.
 
 
!     THE STRESS OUTPUT DATA BLOCK (UAI CODE)
 
!     ADDRESS    DESCRIPTIONS
 
!        1       ELID
!     -------------------------------------------------------
!        2       'CNTR'
!        3       LOWER FIBER DISTANCE
!      4 - 10    STRESSES FOR LOWER POINTS AT ELEMENT CENTER POINT
!       11       UPER  FIBER DISTANCE
!     12 - 18    STRESSES FOR UPPER POINTS AT ELEMENT CENTER POINT
!       19       FIRST GRID POINT NUMBER
!     20 - 35    REPEAT  3 TO 18 ABOVE FOR FIRST  GRID POINT
!     36 - 52    REPAET 19 TO 36 ABOVE FOR SECOND GRID POINT
!     53 - 69    REPAET 19 TO 36 ABOVE FOR THIRD  GRID POINT
 
 
!     THE STRESS OUTPUT DATA BLOCK AT ELEMENT CENTER ONLY, COSMIC
 
!     ADDRESS    DESCRIPTIONS
 
!        1       ELID
!     -------------------------------------------------------
!        2       LOWER FIBER DISTANCE
!      3 -  9    STRESSES FOR LOWER POINTS AT ELEMENT CENTER POINT
!       10       UPER  FIBER DISTANCE
!     11 - 17    STRESSES FOR UPPER POINTS AT ELEMENT CENTER POINT
 
 
 
 INTEGER, INTENT(IN)                      :: numpx
 INTEGER, INTENT(IN)                      :: elid
 INTEGER, INTENT(IN)                      :: igrid(1)
 REAL, INTENT(IN)                         :: thikns(1)
 REAL, INTENT(IN)                         :: z12(2,1)
 REAL, INTENT(IN)                         :: g(6,6)
 REAL, INTENT(IN)                         :: epscsi(6,1)
 REAL, INTENT(IN)                         :: stemp(2)
 REAL, INTENT(IN)                         :: tbar
 REAL, INTENT(IN)                         :: g2alfb(3,1)
 LOGICAL, INTENT(IN OUT)                  :: bendng
 INTEGER, INTENT(IN OUT)                  :: idr(1)
 LOGICAL :: grids, vonms, layer, strcur, stsreq,stnreq,  &
     gridss,vonmss,layers,forreq,temper,tempp1,tempp2, cosmic
 INTEGER :: nstres(1)
 
 REAL :: s1mat(3,3),s2mat(3,3),sigma(3),epss,sigmap(4),  &
     thick,t3ov12,fiber,const, tprime,tsubi
 COMMON /sdr2x7/ dum71(100),stres(100),forsul(200),strin(100)
 COMMON /tmpdat/ temper,tempp1,tempp2
 COMMON /outreq/ stsreq,stnreq,forreq,strcur,grids,vonms,layer  &
     ,               gridss,vonmss,layers
 EQUIVALENCE     (nstres(1),stres(1))
 DATA    cosmic, epss  / .true., 1.0E-11 /
 
 
!     ELEMENT CENTER POINT COMPUTAION ONLY FOR COSMIC,
!     I.E. THE CALLER SHOULD PASS 1 IN NUMPX FOR COSMIC, 4 FOR UAI
 
 nump = numpx
 IF (cosmic) nump = 1
 
 nstres(1) = elid
 
!     START THE LOOP ON EVALUATION POINTS
 
 nump1 = nump - 1
 DO  inplan = 1,nump
   thick  = thikns(inplan)
   t3ov12 = thick*thick*thick/12.0
   
   istres = 1
   IF (cosmic) GO TO 140
   
   istres = (inplan-1)*17 + 2
   nstres(istres) = inplan - 1
   IF (.NOT.grids .OR. inplan <= 1) GO TO 130
   DO  inptmp = 1,nump1
     IF (idr(inptmp) == igrid(inplan)) GO TO 120
   END DO
   CALL errtrc ('SHSTSS  ',100)
   120 istres = inptmp*17 + 2
   nstres(istres) = igrid(inplan)
   130 IF (inplan == 1) nstres(istres) = igrid(inplan)
   
   
!     START THE LOOP ON FIBERS
   
   140 DO  iz = 1,2
     fiber = z12(iz,inplan)
     stres(istres+1) = fiber
     const = 12.0*fiber/thick
     
!     CREATE [S1] AND [S2]
     
     DO  i = 1,3
       DO  j = 1,3
         s1mat(i,j) = g(i  ,j  ) - const*g(i,j+3)
         s2mat(i,j) = g(i+3,j+3) - const*g(i,j+3)
       END DO
     END DO
     
!     EVALUATE STRESSES AT THIS FIBER DISTANCE
     
     DO  i = 1,3
       sigma(i) = 0.0
       DO  j = 1,3
         sigma(i) = sigma(i) + s1mat(i,j)*epscsi(j  ,inplan)  &
             - fiber    * s2mat(i,j)*epscsi(j+3,inplan)
       END DO
     END DO
     
!     IF TEMPERATURES ARE PRESENT, RECORRECT STRESSES FOR THERMAL
!     STRESSES RESULTING FROM TEMPERATURE VALUES AT FIBER DISTANCES.
     
     IF (.NOT.temper .OR. .NOT.bendng) GO TO 250
     IF (.NOT.tempp1) GO TO 180
     tprime = stemp(2   )
     tsubi  = stemp(2+iz)
     IF (ABS(tsubi) < epss) GO TO 250
     tsubi  = tsubi - tprime*fiber
     GO TO 220
     
     180 IF (.NOT.tempp2) GO TO 250
     tsubi = stemp(4+iz)
     IF (ABS(tsubi) < epss) GO TO 250
     DO  ist = 1,3
       sigma(ist) = sigma(ist) - stemp(ist+1)*fiber/t3ov12
     END DO
     
     220 tsubi = tsubi - tbar
     DO  its = 1,3
       sigma(its) = sigma(its) - tsubi*g2alfb(its,inplan)
     END DO
     
!     CLEANUP AND SHIP CORRECTED STRESSES
     
     250 DO  its = 1,3
       IF (ABS(sigma(its)) <= epss) sigma(its) = 0.0
       stres(istres+1+its) = sigma(its)
     END DO
     
!     CALCULATE PRINCIPAL STRESSES
     
     CALL shpsts (sigma,vonms,sigmap)
     stres(istres+5) = sigmap(1)
     stres(istres+6) = sigmap(2)
     stres(istres+7) = sigmap(3)
     stres(istres+8) = sigmap(4)
     
     istres = istres + 8
   END DO
 END DO
 
 RETURN
END SUBROUTINE shstss
