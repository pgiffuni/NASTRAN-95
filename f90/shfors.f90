SUBROUTINE shfors (numpx,elid,igrid,thikns,g,epscsi,qveci,idr)
     
!     TO CALCULATE SHELL ELEMENT FORCES FOR A 2-DL FORMULATION BASE.
 
 
!     INPUT :
!           NUMPX  - NUMBER OF EVALUATION POINTS
!           ELID   - ELEMENT ID
!           IGRID  - ARRAY IF EXTERNAL GRID IDS
!           THIKNS - EVALUATION POINT THICKNESSES
!           G      - 6X6 STRESS-STRAIN MATRIX
!           EPSCSI - CORRECTED STRAINS AT EVALUATION POINTS
!           QVECI  - CALCULATED SHEAR FORCES READY FOR OUTPUT
!           IDR    - REORDERING ARRAY BASED ON EXTERNAL GRID POINT ID'S
!          /OUTREQ/- OUTPUT REQUEST LOGICAL FLAGS
 
!     OUTPUT:
!            FORCES ARE PLACED AT THE PROPER LOCATION IN /SDR2X7/.
 
 
!     THE FORCE RESULTANT OUTPUT DATA BLOCK, UAI CODE
 
!     ADDRESS    DESCRIPTIONS
 
!        1       ELID
!     ------------------------------------------------
!        2       GRID POINT NUMBER OR 'CNTR'
!      3 - 10    FORCES AT ELEMENT CENTER POINT
!     ---------- ABOVE DATA REPEATED 3 TIMES
!                FOR GRID POINTS
 
 
!     THE FORCE RESULTANT OUTPUT DATA BLOCK AT ELEMETN CENTER, COSMIC
 
!     ADDRESS    DESCRIPTIONS
 
!        1       ELID
!     ------------------------------------------------
!      2 - 9     FORCES AT ELEMENT CENTER POINT
 
 
 
 INTEGER, INTENT(IN)                      :: numpx
 INTEGER, INTENT(IN)                      :: elid
 INTEGER, INTENT(IN)                      :: igrid(1)
 REAL, INTENT(IN)                         :: thikns(1)
 REAL, INTENT(IN)                         :: g(6,6)
 REAL, INTENT(IN OUT)                     :: epscsi(6,1)
 REAL, INTENT(IN)                         :: qveci(2,1)
 INTEGER, INTENT(IN OUT)                  :: idr(1)
 LOGICAL :: grids,vonms,layer,strcur,stsreq,stnreq,forreq  &
     ,               gridss,vonmss,layers,cosmic
 INTEGER :: nfors(1)
 REAL :: thick,thick2,t3ov12,dforce(8),GT(6,6)
 COMMON /sdr2x7/ dum71(100),stres(100),forsul(200),strin(100)
 COMMON /outreq/ stsreq,stnreq,forreq,strcur,grids,vonms,layer  &
     ,               gridss,vonmss,layers
 EQUIVALENCE     (nfors(1),forsul(1))
 DATA    cosmic/ .true. /
 
 
!     ELEMENT CENTER POINT COMPUTAION ONLY FOR COSMIC
!     IE. CALLER SHOULD PASS 1 IN NUMPX FOR COSMIC, 4 FOR UAI
 
 nump = numpx
 IF (cosmic) nump = 1
 
 nfors(1)  = elid
 
!     START THE LOOP ON EVALUATION POINTS
 
 nump1 = nump - 1
 DO  inplan = 1,nump
   thick  = thikns(inplan)
   thick2 = thick*thick
   t3ov12 = thick2*thick/12.0
   
   iforce = 1
   IF (cosmic) GO TO 250
   
   iforce = (inplan-1)*9 + 2
   IF (.NOT.(grids .AND. gridss) .OR. inplan <= 1) GO TO 230
   DO  inptmp = 1,nump1
     IF (idr(inptmp) == igrid(inplan)) EXIT
   END DO
   220 CONTINUE
   iforce = inptmp*9 + 2
   nfors(iforce) = igrid(inplan)
   GO TO 240
   230 nfors(iforce) = inplan - 1
   240 IF (inplan == 1) nfors(iforce) = igrid(inplan)
   
!     MODIFY [G], THEN CALCULATE FORCES AND MOMENTS
   
   250 DO  ig = 1,3
     DO  jg = 1,3
       GT(ig  ,jg  ) = thick *g(ig  ,jg  )
       GT(ig+3,jg  ) = thick2*g(ig+3,jg  )
       GT(ig  ,jg+3) = thick2*g(ig  ,jg+3)
       GT(ig+3,jg+3) = t3ov12*g(ig+3,jg+3)
     END DO
   END DO
   CALL gmmats (GT,6,6,0, epscsi(1,inplan),6,1,0, dforce(1))
   
!     OUTPUT QX AND QY (WE HAVE CALCULATED QY AND QX)
   
   dforce(7) = qveci(2,inplan)
   dforce(8) = qveci(1,inplan)
   
!     SHIP OUT
   
   DO  ifor = 1,8
     forsul(iforce+ifor) = dforce(ifor)
   END DO
 END DO
 
 RETURN
END SUBROUTINE shfors
