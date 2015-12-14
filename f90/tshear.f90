SUBROUTINE tshear
     
!     ELEMENT TEMPERATURE AND DEFORMATION LOADING FOR THE SHEAR PANEL.
 
!     FORMULATION IS THAT OF A PSEUDO-ROD ON EACH EDGE OF THE SHEAR
!     PANEL.
 
!     ECPT( 1)         - ELEMENT ID
!     ECPT( 2 THRU 5)  - 4 GRID SILS
!     ECPT( 6)         - MATERIAL ID
!     ECPT( 7)         - THICKNESS
!     ECPT( 8)         - NON-STRUCTURAL MASS
!     ECPT( 9 THRU 24) - 4 POINTS (CSID,X,Y,Z) REPEATS
!     ECPT(25)         - ELEMENT TEMPERATURE
!     ECPT(26)         - F1 EFFECTIVENESS FACTOR DIRECTION 1, (NOT USED)
!     ECPT(27)         - F2 EFFECTIVENESS FACTOR DIRECTION 2, (NOT USED)
 
 INTEGER :: ncsid(4,4)
 COMMON /ssgett/ eltype   ,oldel    ,eorflg   ,endid   ,bufflg   ,  &
     itemp    ,ideft    ,idefm
 COMMON /trimex/ ecpt(1)  ,isils(4) ,mid      ,thick   ,fmu      ,  &
     csid(4,4),eltemp   ,f12(2)
 COMMON /matin / matid    ,inflag   ,temp     ,stress  ,sinth    , costh
 COMMON /matout/ e1       ,g        ,nu       ,rho     ,alpha    ,  &
     to1      ,GE       ,sigmat   ,sigmac  ,sigmas
 COMMON /ssgwrk/ vec(3,4) ,xl(4)    ,diag1(3) ,diag2(3),ti(16)   ,  &
     pa       ,tsq      ,veca(3)  ,vecb(3) ,area     ,  &
     vmag     ,i        ,j        ,ia      ,ib       , i12      ,tbar     ,in
 COMMON /zzzzzz/ pg(1)
 EQUIVALENCE     (ncsid,csid)
 
 f12(1) = 1.00
 f12(2) = 1.00
 
 IF (f12(1) == 0.0 .AND. f12(2) == 0.0) RETURN
 
!     MATERIAL DATA ACQUISITION
 
 inflag = 1
 matid  = mid
 temp   = eltemp
 CALL mat (ecpt(1))
 
!     GRID POINT TEMPERATURES
 
 IF (itemp == 0) THEN
   GO TO   800
 END IF
 100 CALL ssgetd (ecpt(1),ti,4)
 
!     ELEMENT DEFORMATION (NOT USED)
 
!     4 NORMALIZED EDGE VECTORS AND LENGTHS
 
 DO  i = 1,4
   igrid2 = i + 1
   IF (i == 4) igrid2 = 1
   
   DO  j = 1,3
     vec(j,i) = csid(j+1,i) - csid(j+1,igrid2)
   END DO
   
   CALL norm (vec(1,i),xl(i))
 END DO
 
 IF (f12(1) > 1.01 .AND. f12(2) > 1.01) GO TO 500
 
!     PROJECTED AREA IS NEEDED. FIRST OBTAIN THE DIAGONAL VECTORS.
 
 DO  i = 1,3
   diag1(i) = csid(i+1,3) - csid(i+1,1)
   diag2(i) = csid(i+1,4) - csid(i+1,2)
 END DO
 
!     NORMAL VECTOR (DIAG1 X DIAG2)
 
 CALL saxb (diag1,diag2,diag2)
 CALL norm (diag2,vkl)
 pa = 0.5*vkl
 
!     LOOP THROUGH LOADS ON 4 EDGES.
 
 500 tsq = thick*thick
 DO  i = 1,4
   i12 = MOD(i,2)
   IF (i12 == 0) i12 = 2
   ia = i
   ib = ia + 1
   IF (i == 4) ib = 1
   
!     TEMPERATURE
   
   tbar = (ti(ia+1) + ti(ib+1))/2.0 - to1
   
!     EXTENSIONAL AREA
   
   IF (f12(i12) <= 1.01) GO TO 550
   area = 0.50*f12(i12)*tsq
   GO TO 560
   550 area = f12(i12)*pa*thick/(xl(i12) + xl(i12+2))
   
   560 vmag = e1*area*alpha*tbar
   DO  j = 1,3
     veca(j) = vmag*vec(j,i)
     vecb(j) =-veca(j)
   END DO
   
   IF (ncsid(1,ib) == 0) THEN
     GO TO   590
   END IF
   580 CALL basglb (vecb(1),vecb(1),isils(ib),csid(1,ib))
   590 in = isils(ib) - 1
   DO  j = 1,3
     in = in + 1
     pg(in) = pg(in) + vecb(j)
   END DO
   
   IF (ncsid(1,ia) == 0) THEN
     GO TO   630
   END IF
   620 CALL basglb (veca(1),veca(1),isils(ia),csid(1,ia))
   630 in = isils(ia) - 1
   DO  j = 1,3
     in = in + 1
     pg(in) = pg(in) + veca(j)
   END DO
   
 END DO
 800 RETURN
END SUBROUTINE tshear
