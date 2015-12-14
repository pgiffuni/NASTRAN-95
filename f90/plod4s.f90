SUBROUTINE plod4s
     
!     ROUTINE TO PROCESS PLOAD4 BULK DATA TO CREATE LOADS ON
!     QUAD4 ELEMENTS
 
!     SINGLE PRECISION VERSION.
 
!     GRID POINT NUMBERING IS COUNTER-CLOCKWISE
!     GRIDS 1,2,3, AND 4 ARE AT THE CORNERS
 
 INTEGER :: nest(125),iorder(4),sil(4),ksil(4),nout,cid,  &
     nsurf,swp,sysbuf,iz(1),itgrid(4,4),ibgpdt(4,4), islt(11),nogo
 REAL :: pe(3,4),bgpdt(4,4),ppp(4),nv(3),nvx(3),locate(3)
 REAL :: dpe(3,4),weight,xsi,eta,eps,area,shp(4),  &
     dshp(8),tmpshp(4),dshptp(8),vi(3),vj(3),v3t(3), gauss(3),wtgaus(3)
 COMMON /loadx / idum1(4),cstm,idum2(13),icm
 COMMON /pindex/ best(45),slt(11)
 COMMON /zzzzzz/ z(1)
 COMMON /system/ sysbuf,nout,nogo
 EQUIVALENCE     (iz(1)  ,z(1)   ), (slt(1),islt(1) )
 EQUIVALENCE     (nest(1),best(1)), (numint,nest(25))
 EQUIVALENCE     (sil(1) ,nest(2)), (bgpdt(1,1),best(29))
 EQUIVALENCE     (ibgpdt(1,1),bgpdt(1,1))
 EQUIVALENCE     (pe(1,1),dpe(1,1))
 
 DATA    ndof  / 3 /
 
!                 EST LISTING
!     ----------------------------------------------------------
!      1          EID
!      2 THRU 5   SILS, GRIDS 1 THRU 4
!      6 THRU 9   T (MEMBRANE), GRIDS 1 THRU 4
!     10          THETA (MATERIAL)
!     11          TYPE FLAG FOR WORD 10
!     12          ZOFF  (OFFSET)
!     13          MATERIAL ID FOR MEMBRANE
!     14          T (MEMBRANE)
!     15          MATERIAL ID FOR BENDING
!     16          I FACTOR (BENDING)
!     17          MATERIAL ID FOR TRANSVERSE SHEAR
!     18          FACTOR FOR T(S)
!     19          NSM (NON-STRUCTURAL MASS)
!     20 THRU 21  Z1, Z2  (STRESS FIBRE DISTANCES)
!     22          MATERIAL ID FOR MEMBRANE-BENDING COUPLING
!     23          THETA (MATERIAL) FROM PSHELL CARD
!     24          TYPE FLAG FOR WORD 23
!     25          INTEGRATION ORDER
!     26          THETA (STRESS)
!     27          TYPE FLAG FOR WORD 26
!     28          ZOFF1 (OFFSET)  OVERRIDDEN BY EST(12)
!     29 THRU 44  CID,X,Y,Z - GRIDS 1 THRU 4
!     45          ELEMENT TEMPERATURE
 
 
!          DATA FROM THE PLOAD4 CARD DESCRIBED HERE
!     ----------------------------------------------------------
!     EID  ELEMENT ID
!     P1,P2,P3,P4 CORNER GRID POINT PRESSURES PER UNIT SURFACE AREA
!     G1,G3  DEFINES QUADRILATERAL SURFACE OF HEXA, QUAD8, AND
!                 PENTA SURFACES ON WHICH PRESSURE LOADS EXIST
!                 OTHERWISE SURFACE IS TRIANGULAR IF G3 IS ZERO OR
!                 BLANK, SURFACE IS TRIANGULAR
!     CID  COORDINATE SYSTEM FOR DEFINITION OF PRESSURE VECTOR
!     N1,N2,N3 COMPONENTS OF PRESSURE DIRECTION VECTOR IF CID
!                 'BLANK' OR ZERO, THE PRESSURE ACTS NORMAL TO THE
!                 SURFACE OF THE ELEMENT
 
!     EQUIVALENT NUMERICAL INTEGRATION POINT LOADS PP(III) ARE
!     OBTAINED VIA BI-LINEAR INTERPOLATION
!*****
!     GENERAL INITIALIZATION.
!*****
!     BEST(45) IS THE DATA FOR EST WHICH IS READ IN EXTERN AND IS
!     READY TO BE USED.
 
!     READ FROM PLOAD4 CARDS
!     P1 = PPP(1)
!     P2 = PPP(2)
!     P3 = PPP(3)
!     P4 = PPP(4)
!     CID,N1,N2,N3
!*****
 
 
!     X WILL BE THE LENGTH OF THE PRESSURE VECTOR FOR NORMALIZATION.
!     NV(I) WILL BE THE NORMALIZED PRESSURE VECTOR
 
 x = 0.0
 DO  i = 1,4
   ppp(i) = slt(i+1)
 END DO
 DO  i = 1,3
   nv(i)  = slt(i+8)
   x = x + nv(i)**2
 END DO
 cid = islt(8)
 
 IF (x == 0.0) GO TO 40
 x = SQRT( x )
 DO  i = 1,3
   nv(i) = nv(i) / x
 END DO
 
 40 ncrd = 3
!*****
!     PERFORM TEST FOR PRESENCE OF CONSTANT PRESSURE SET SWP
!*****
 swp = 1
 IF (ppp(2) == 0. .AND. ppp(3) == 0. .AND. ppp(4) == 0.) swp = 0
 nsurf = 4
 
!     THE ARRAY IORDER STORES THE ELEMENT NODE ID IN
!     INCREASING SIL ORDER.
 
!     IORDER(1) = NODE WITH LOWEST  SIL NUMBER
!     IORDER(4) = NODE WITH HIGHEST SIL NUMBER
 
!     ELEMENT NODE NUMBER IS THE INTEGER FROM THE NODE LIST G1,G2,G3,G4.
!     THAT IS, THE 'I' PART OF THE 'GI' AS THEY ARE LISTED ON THE
!     CONNECTIVITY BULK DATA CARD DESCRIPTION.
 
 ksild = 99999995
 DO  i=1,4
   iorder(i) = 0
   ksil(i) = sil(i)
 END DO
 DO  i=1,4
   itemp = 1
   isil  = ksil(1)
   DO  j=2,4
     IF (isil <= ksil(j)) CYCLE
     itemp = j
     isil  = ksil(j)
   END DO
   iorder(i) = itemp
   ksil(itemp) = 99999999
 END DO
!*****
!     ADJUST EST DATA
 
!     USE THE POINTERS IN IORDER TO COMPLETELY REORDER THE
!     GEOMETRY DATA INTO INCREASING SIL ORDER.
!*****
 DO  i=1,4
   ksil(i) = sil(i)
   DO  j=1,4
     itgrid(j,i) = ibgpdt(j,i)
   END DO
 END DO
 DO  i=1,4
   ipoint = iorder(i)
   sil(i) = ksil(ipoint)
   DO  j=1,4
     ibgpdt(j,i) = itgrid(j,ipoint)
   END DO
 END DO
 
 nvct = ncrd*4
 eps  = 0.001
!*****
!     SET VALUES FOR NUMERICAL INTEGRATION POINTS AND WEIGHT FACTORS
 
!     DEFAULT INTEGRATION ORDER IS 2X2
!*****
 numint = 2
 gauss(1) = -0.57735026918962
 gauss(2) = +0.57735026918962
 wtgaus(1) = 1.0
 wtgaus(2) = 1.0
!*****
!     ZERO OUT THE LOAD ROW SET
!*****
 DO  i=1,ndof
   DO  j=1,4
     dpe(i,j) = 0.0
   END DO
 END DO
!*****
!     SET UP THE LOOPS FOR NUMERICAL INTEGRATION
!*****
 DO  ieta=1,numint
   eta = gauss(ieta)
   DO  ixsi=1,numint
     xsi = gauss(ixsi)
     weight = wtgaus(ixsi)*wtgaus(ieta)
     p = 0.0
!*****
!     P1,P2,P3,P4 ARE THE GRID POINT PRESSURE LOADS PER UNIT
!     AREA FROM THE PLOAD4 CARD.  THESE WILL BE USED WITH A
!     BILINEAR SHAPE FUNCTION ROUTINE TO CALCULATE THE NODAL
!     LOADS.
     
!     BILINEAR CASE WHERE THE VALUES OF XSI,ETA ARE INPUT IN
!     EXPLICIT FORM DEPENDING UPON WHICH NUMERICAL INTEGRATION
!     SCHEME IS BEING USED.
     
     
!     NSURF IS AN INTEGER WHICH KEEPS TRACK OF THE SURFACE TYPE
!              NSURF = 3 . . .  TRIANGULAR SURFACE
!              NSURF = 4 . . .  QUADRILATERAL SURFACE
     
!*****
!     CALL SHAPE FCN. ROUTINE FOR THE BILINEAR QUAD4.  INPUT IS
!     XSI,ETA,III AND EVALUATION OF SHAPE FCN. AT INTEG.PTS
!     WILL BE PERFORMED.
!*****
     CALL q4shps (xsi,eta,shp,dshp)
     
     IF (swp == 0) p=ppp(1)
     IF (swp == 0) GO TO 200
     
     DO  iii=1,nsurf
       p = p + shp(iii)*ppp(iii)
     END DO
     200 CONTINUE
     
!     SORT THE SHAPE FUNCTIONS AND THEIR DERIVATIVES INTO SIL ORDER.
     
     DO  i=1,4
       tmpshp(i) = shp(i)
       dshptp(i) = dshp(i)
       dshptp(i+4) = dshp(i+4)
     END DO
     DO  i=1,4
       kk = iorder(i)
       shp (i  ) = tmpshp(kk)
       dshp(i  ) = dshptp(kk)
       dshp(i+4) = dshptp(kk+4)
     END DO
!*****
!     COMPUTE THE UNIT NORMALS V3T AT EACH GRID POINT.  THESE WILL
!     BE USED TO GET COMPONENTS OF PRESSURE VECTOR ACTING NORMAL TO
!     THE SURFACE.  AREA CALCULATION CHECKS THE GEOMETRY OF THE
!     ELEMENT.
!*****
     DO  i=1,3
       vi(i) = 0.0
       vj(i) = 0.0
       DO  j=1,4
         ii=i+1
         vi(i) = vi(i) + bgpdt(ii,j) * dshp(j)
         vj(i) = vj(i) + bgpdt(ii,j) * dshp(j+4)
       END DO
     END DO
     
!     CHECK FOR USER INPUT VECTOR TO ROTATE LOADS
     
     CALL saxb (vi,vj,v3t)
     area = SQRT(v3t(1)**2+v3t(2)**2+v3t(3)**2)
     IF (area > 0.0) GO TO 300
     
     WRITE (nout,240) nest(1)
     240 FORMAT ('0*** SYSTEM FATAL ERROR.  BAD GEOMETRY DETECTED FOR ',  &
         'QUAD4 ELEMENT ',i8,' WHILE PROCESSING PLOAD4 DATA.')
     nogo = 1
     RETURN
     
     300 CONTINUE
     IF (x == 0.0) GO TO 330
     
!     CHECK FOR NON-ZERO CID AND NEED TO ROTATE USER'S VECTOR
     
     IF (cid == 0) GO TO 320
     
!     COMPUTE THE LOCATION OF THE INTEGRATION POINT SO THAT WE CAN
!     ROTATE THE USER VECTOR PER CID. THIS LOCATION REQUIRED ONLY IF
!     CID IS CYLINDRICAL OR SPHERICAL.
     
     locate(1) = 0.
     locate(2) = 0.
     locate(3) = 0.
     DO  j=1,4
       locate(1) = locate(1) + bgpdt(2,j) * shp(j)
       locate(2) = locate(2) + bgpdt(3,j) * shp(j)
       locate(3) = locate(3) + bgpdt(4,j) * shp(j)
     END DO
     CALL glbbas (nv(1),nvx(1),locate(1),cid)
     
!     NOW ROTATE THE PRESSURE LOAD
     
     v3t(1) = nvx(1) * area
     v3t(2) = nvx(2) * area
     v3t(3) = nvx(3) * area
     GO TO 330
     320 CONTINUE
     
!     NOW ROTATE THE PRESSURE LOAD
     
     v3t(1) = nv(1) * area
     v3t(2) = nv(2) * area
     v3t(3) = nv(3) * area
     330 CONTINUE
     
!*****
!     COMPUTE THE CONTRIBUTION TO THE LOAD MATRIX FROM THIS
!     INTEGRATION POINT AS NT * P * V3T
!*****
     DO  i=1,4
       DO  j=1,ndof
         dpe(j,i) = dpe(j,i) + weight * p * shp(i) * v3t(j)
       END DO
     END DO
   END DO
 END DO
!*****
!     END OF NUMERICAL INTEGRATION LOOPS
 
!     MOVE DATA FROM REAL ARRAY DPE TO SINGLE PRECISION PE
!     (NO MOVE, SINCE DPE IS EQUIVALENT TO PE)
!*****
!     DO 400 J=1,4
!     PE(1,J) = DPE(1,J)
!     PE(2,J) = DPE(2,J)
!     PE(3,J) = DPE(3,J)
! 400 CONTINUE
!*****
!     ADD ELEMENT LOAD TO OVERALL LOAD.
!*****
 jb = 25
 DO  j=1,4
   jb = jb + 4
   IF (nest(jb) == 0) GO TO 410
   CALL basglb (pe(1,j),pe(1,j),best(jb+1),nest(jb))
   410 CONTINUE
   jp = sil(j) - 1
   DO  i=1,3
     z(jp+i) = z(jp+i) + pe(i,j)
   END DO
 END DO
 RETURN
END SUBROUTINE plod4s
