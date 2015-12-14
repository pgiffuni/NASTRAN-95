SUBROUTINE t3pl4s
     
!     SINGLE PRECISION ROUTINE TO PROCESS PLOAD4 PRESSURE DATA AND
!     GENERATE EQUIVALENT NODAL LOADS FOR A TRIA3 ELEMENT.
 
!     WAS NAMED T3PRSS (LOADVC,RPDATA,IPDATA) IN UAI
 
!                 EST  LISTING
 
!        WORD     TYP       DESCRIPTION
!     ----------------------------------------------------------------
!     ECT:
!         1        I   ELEMENT ID, EID
!         2-4      I   SIL LIST, GRIDS 1,2,3
!         5-7      R   MEMBRANE THICKNESSES T, AT GRIDS 1,2,3
!         8        R   MATERIAL PROPERTY ORIENTAION ANGLE, THETA
!               OR I   COORD. SYSTEM ID (SEE TM ON CTRIA3 CARD)
!         9        I   TYPE FLAG FOR WORD 8
!        10        R   GRID OFFSET, ZOFF
!    EPT:
!        11        I   MATERIAL ID FOR MEMBRANE, MID1
!        12        R   ELEMENT THICKNESS,T (MEMBRANE, UNIFORMED)
!        13        I   MATERIAL ID FOR BENDING, MID2
!        14        R   MOMENT OF INERTIA FACTOR, I (BENDING)
!        15        I   MATERIAL ID FOR TRANSVERSE SHEAR, MID3
!        16        R   TRANSV. SHEAR CORRECTION FACTOR, TS/T
!        17        R   NON-STRUCTURAL MASS, NSM
!        18-19     R   STRESS FIBER DISTANCES, Z1,Z2
!        20        I   MATERIAL ID FOR MEMBRANE-BENDING COUPLING, MID4
!        21        R   MATERIAL ANGLE OF ROTATION, THETA
!               OR I   COORD. SYSTEM ID (SEE MCSID ON PSHELL CARD)
!                      (DEFAULT FOR WORD 8)
!        22        I   TYPE FLAG FOR WORD 21 (DEFAULT FOR WORD 9)
!        23        I   INTEGRATION ORDER FLAG
!        24        R   STRESS ANGLE OF RATATION, THETA
!               OR I   COORD. SYSTEM ID (SEE SCSID ON PSHELL CARD)
!        25        I   TYPE FLAG FOR WORD 24
!        26        R   OFFSET, ZOFF1 (DEFAULT FOR WORD 10)
!    BGPDT:
!        27-38   I/R   CID,X,Y,Z  FOR GRIDS 1,2,3
!    ETT:
!        39        I   ELEMENT TEMPERATURE
 
 
!    DATA IN THE PLOAD4 ENTRY, 11 WORDS IN ISLT ARRAY
 
!       EID - ELEMENT ID, IPDATA(0)=ISLT(1)
!       PPP - CORNER GRID POINT PRESSURES PER UNIT SURFACE AREA,
!             RPDATA (1-4)
!       DUM - DUMMY DATA WORDS, IPDATA (5-6)
!       CID - COORDINATE SYSTEM FOR DEFINITION OF PRESSURE VECTOR,
!             IPDATA(7)
!       NV  - PRESSURE DIRECTION VECTOR, RPDATA(8-10)
!           - IF CID IS BLANK OR ZERO, THE PRESSURE ACTS NORMAL TO THE
!             SURFACE OF THE ELEMENT.
 
!     EQUIVALENT NUMERICAL INTEGRATION POINT LOADS PP(III) ARE OBTAINED
!     VIA BI-LINEAR INTERPOLATION
 
 
 LOGICAL :: constp,sheart,normal
 INTEGER :: ipdata(7),islt(1),igpdt(4,3),sil(3),iorder(3),  &
     elid,cid,sysbuf,nout,nogo
 REAL :: gpth(3),bgpdt(4,3),nv(3),nvx(3),locate(3), pe(3,3),rpdata(1),loadvc
 REAL :: dpe(3,3),shp(3),weight,detjac,x,unv(3),v3t(3),  &
     p,ppp(3),bterms(6),bmatrx(162),egpdt(4,3),  &
     cente(3),gpnorm(4,3),epnorm(4,3),teb(9),tub(9),  &
     dgpth(3),th,avgthk,aic(1),edglen(3),lx,ly
 COMMON /system/  sysbuf,nout,nogo
 COMMON /zzzzzz/  loadvc(1)
 COMMON /pindex/  est(45),slt(11)
 EQUIVALENCE      (slt( 1),islt(1)),(est( 1),elid),(est(2),sil(1)),  &
     (est( 5),gpth(1)),(est(12),elth), (est(27),bgpdt(1,1),igpdt(1,1)),  &
     (slt( 2),ipdata(1) ,rpdata(1))
 
 
!     INITIALIZE
 
 weight = 1.0/6.0
 sheart = .false.
 nnode  = 3
 ndof   = 3
 DO  i = 1,ndof
   DO  j = 1,nnode
     dpe(i,j) = 0.0
   END DO
 END DO
 
!     GET THE PRESSURE INFORMATION
 
!     EST (45 WORDS) AND SLT (11 WORDS) ARE THE DATA FOR EST AND SLT
!     WHICH ARE READ IN BY EXTERN AND ARE READY TO BE USED
 
 
!     IF ISLT(1).GT.0, GET THE PLOAD4 DATA FROM THE PROCESSED PLOAD2
!                      INFORMATION IN ARRAY SLT.
!                      (NOT AVAILABLE IN COSMIC.NASTRAN)
!     IF ISLT(1).LT.0, GET THE PLOAD4 DATA FROM THE ORIGINAL PLOAD4
!                      INFORMATION IN ARRAY PDATA.
 
 IF (islt(1) < 0) GO TO 20
 normal = .true.
 constp = .true.
 p = slt(2)
 GO TO 60
 
 20 DO  i = 1,nnode
   ppp(i) = rpdata(i)
 END DO
 constp = ppp(2) == 0.0 .AND. ppp(3) == 0.0
 IF (constp) p = ppp(1)
 cid = ipdata(7)
 
!     GET THE DIRECTION VECTOR AND NORMALIZE IT
 
 x = 0.0
 DO  i = 1,nnode
   unv(i) = rpdata(i+7)
   x = x + unv(i)*unv(i)
 END DO
 
 normal = .true.
 IF (x <= 0.0) GO TO 60
 normal = .false.
 x = SQRT(x)
 DO  i = 1,nnode
   nv(i) = unv(i)/x
 END DO
 
!     SET UP THE ELEMENT FORMULATION
 
 60 CALL t3sets (ierr,sil,igpdt,elth,gpth,dgpth,egpdt,gpnorm,epnorm,  &
     iorder,teb,tub,cente,avgthk,lx,ly,edglen,elid)
 IF (ierr /= 0) GO TO 200
 
!     START THE LOOP ON INTEGRATION POINTS
 
 DO  ipt = 5,7
   
   CALL t3bmgs (ierr,sheart,ipt,iorder,egpdt,dgpth,aic,th,detjac,shp,  &
       bterms,bmatrx)
   IF (ierr /= 0) GO TO 200
   
!     CALCULATE THE PRESSURE AT THIS POINT
   
   IF (constp) GO TO 80
   p = 0.0
   DO  i = 1,nnode
     p = p + shp(i)*ppp(i)
   END DO
   
!     SET THE DIRECTION OF PRESSURE AT THIS POINT.
!     THE RESULTING VECTOR MUST BE IN THE BASIC COORD. SYSTEM
   
   80 IF (.NOT.normal) GO TO 90
   v3t(1) = teb(7)*detjac
   v3t(2) = teb(8)*detjac
   v3t(3) = teb(9)*detjac
   GO TO 120
   
   90 IF (cid /= 0) GO TO 100
   v3t(1) = nv(1)*detjac
   v3t(2) = nv(2)*detjac
   v3t(3) = nv(3)*detjac
   GO TO 120
   
!     FOR NON-ZERO CID, COMPUTE THE LOCATION OF THE INTEGRATION POINT SO
!     THAT WE CAN ROTATE THE USER VECTOR PER CID.  THIS LOCATION IS
!     REQUIRED ONLY IF CID IS CYLINDRICAL OR SPHERICAL.
   
   100 locate(1) = 0.0
   locate(2) = 0.0
   locate(3) = 0.0
   DO  j  = 1,nnode
     locate(1) = locate(1) + bgpdt(2,j)*shp(j)
     locate(2) = locate(2) + bgpdt(3,j)*shp(j)
     locate(3) = locate(3) + bgpdt(4,j)*shp(j)
   END DO
   
!     NOW ROTATE THE VECTOR
   
   CALL glbbas (nv(1),nvx(1),locate(1),cid)
   v3t(1) = nvx(1)*detjac
   v3t(2) = nvx(2)*detjac
   v3t(3) = nvx(3)*detjac
   
!     COMPUTE THE CONTRIBUTION TO THE LOAD MATRIX FROM THIS INTEGRATION
!     POINT AS NT*P*V3T
   
   120 DO  i = 1,nnode
     DO  j = 1,ndof
       dpe(j,i) = dpe(j,i) + weight*p*shp(i)*v3t(j)
       pe(j,i)  = dpe(j,i)
     END DO
   END DO
   
 END DO
 
!     END OF NUMERICAL INTEGRATION LOOP
!     ADD ELEMENT LOAD TO OVERALL LOAD.
 
 DO  j = 1,nnode
   IF (igpdt(1,j) /= 0) CALL basglb (pe(1,j),pe(1,j),bgpdt(2,j), igpdt(1,j))
   
   jp = sil(j) - 1
   DO  i = 1,ndof
     loadvc(jp+i) = loadvc(jp+i) + pe(i,j)
   END DO
 END DO
 GO TO 250
 
!     FATAL ERROR
 
 200 islt(1) = IABS(islt(1))
 CALL mesage (30,224,islt(1))
 nogo = 1
 
 250 RETURN
END SUBROUTINE t3pl4s
