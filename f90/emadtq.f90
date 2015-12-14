SUBROUTINE emadtq(narg,mass)
!  THE EMG MASS DOUBLE PRECISION ROUTINE FOR TRI S, QUAD S, TWIST +
!     SHEAR ELEMENTS
 
!  THIS SUBROUTINE CALCULATES THE MASS MATRIX FOR THE ELEMENTS LISTED
!   BELOW
 
!   NOTE THAT THE OUTPUT MASS MATRIX IS NOT ORDERED BY INCREASING SIL
!     DOUBLE PRECISION VERSION
 
!     ******************************************************************
!     E C P T    L I S T I N G S
!     **************************
!         MTWIST              MQDMEM                        MTRMEM
!         MSHEAR    MQUAD1    MQUAD2    MTRIA1    MTRBSC    MTRIA2
! **********************************************************************
! ECPT( 1)ELEM. ID  ELEM. ID  ELEM. ID  ELEM. ID  ELEM. ID  ELEM. ID
! ECPT( 2)GR.PT. A  GR.PT. A  GR.PT. A  GR.PT. A  GR.PT. A  GR.PT. A
! ECPT( 3)GR.PT. B  GR.PT. B  GR.PT. B  GR.PT. B  GR.PT. B  GR.PT. B
! ECPT( 4)GR.PT. C  GR.PT. C  GR.PT. C  GR.PT. C  GR.PT. C  GR.PT. C
! ECPT( 5)GR.PT. D  GR.PT. D  GR.PT. D  THETA     THETA     THETA
! ECPT( 6)MAT ID    THETA     THETA     MAT ID 1  MAT ID 1  MAT ID
! ECPT( 7)T         MAT ID 1  MAT ID    T1        I         T
! ECPT( 8)N S MASS  T1        T         MAT ID 2  MAT ID 2  NS MASS
! ECPT( 9)CSID 1    MAT ID 2  N S MASS  I         T2        CSID 1
! ECPT(10)X1        I         CSID 1    MAT ID 3  N S MASS  X1
! ECPT(11)Y1        MAT ID 3  X1        T2        Z1        Y1
! ECPT(12)Z1        T2        Y1        N S MASS  Z2        Z1
! ECPT(13)CSID 2    N S MASS  Z1        Z1        CSID 1    CSID 2
! ECPT(14)X2        Z1        CSID 2    Z2        X1        X2
! ECPT(15)Y2        Z2        X2        CSID 1    Y1        Y2
! ECPT(16)Z2        CSID 1    Y2        X1        Z1        Z2
! ECPT(17)CSID 3    X1        Z2        Y1        CSID 2    CSID 3
! ECPT(18)X3        Y1        CSID 3    Z1        X2        X3
! ECPT(19)Y3        Z1        X3        CSID 2    Y2        Y3
! ECPT(20)Z3        CSID 2    Y3        X2        Z2        Z3
! ECPT(21)CSID 4    X2        Z3        Y2        CSID 3    TEMP
! ECPT(22)X4        Y2        CSID 4    Z2        X3
! ECPT(23)Y4        Z2        X4        CSID 3    Y3
! ECPT(24)Z4        CSID 3    Y4        X3        Z3
! ECPT(25)TEMP      X3        Z4        Y3        TEMP
! ECPT(26)          Y3        TEMP      Z3
! ECPT(27)          Z3                  TEMP
! ECPT(28)          CSID 4
! ECPT(29)          X4
! ECPT(30)          Y4
! ECPT(31)          Z4
! ECPT(32)          TEMP
! **********************************************************************
 
 
 INTEGER, INTENT(IN)                      :: narg
 DOUBLE PRECISION, INTENT(OUT)            :: mass(100)
 DOUBLE PRECISION :: v1(3),v2(3),v1xv2(3), fmu, t, area , term, rhod
 DIMENSION necpt (7)
 INTEGER :: heat
 
 COMMON /hmtout/ cp
 COMMON /matin/ matid,inflag,eltemp
 COMMON /matout/ rho
 COMMON/system/ ksystm(55),heat
 
 COMMON / emgest/ ecpt(100)
 
 EQUIVALENCE ( necpt(1) , ecpt(1) )
 EQUIVALENCE (iflag , ecpt(8) )
 DATA pi23/2.0943952/
 
!     THIS ROUTINE COMPUTES A MASS MATRIX OF THE FOLLOWING FORM.
 
!   MASS MATRIX = (T1,T1,T1,T2,T2,T2,T3,T3,T3,IF REQ-D T4,T4,T4)  )
 
!                   *******************
!                   NTYPE = 1  -MQDMEM-
!                   NTYPE = 1  -MQUAD2-
!                   NTYPE = 2  -MQUAD1-
!                   NTYPE = 3  -MTRBSC-
!                   NTYPE = 3  -MTRPLT-
!                   NTYPE = 4  -MTRMEM-
!                   NTYPE = 4  -MTRIA2-
!                   NTYPE = 5  -MTRIA1-
!                   NTYPE = 6  -MSHEAR-
!                   NTYPE = 6  -MTWIST-
!                   NTYPE = 7  -MQDPLT-
!                   *******************
 
 ntype = narg
 ndof = 3
 
!            -MQDMEM-      -MTRPLT-MTRMEM-      -MTWIST-
!            -MQUAD2-MQUAD1-MTRBSC-MTRIA2-MTRIA1-MSHEAR-MQDPLT-
 SELECT CASE ( ntype )
   CASE (    1)
     GO TO 10
   CASE (    2)
     GO TO 20
   CASE (    3)
     GO TO 30
   CASE (    4)
     GO TO 40
   CASE (    5)
     GO TO 50
   CASE (    6)
     GO TO 60
   CASE (    7)
     GO TO 70
 END SELECT
 
 10 ncsid = 10
 ngrids = 4
 matid = necpt(7)
 t     = ecpt(8)
 fmu   = ecpt(9)
 GO TO 80
 
 20 ncsid = 16
 ngrids = 4
 matid = necpt(7)
 t     = ecpt(8)
 fmu   = ecpt(13)
 GO TO 80
 
 30 ncsid =  13
 ngrids = 3
 matid = necpt( 6)
 t     =  0.0E0
 fmu   =  ecpt(10)
 GO TO 80
 
 40 ncsid = 9
 ngrids = 3
 matid = necpt(6)
 t     = ecpt(7)
 fmu   = ecpt(8)
 GO TO 80
 
 50 ncsid = 15
 ngrids = 3
 matid = necpt( 6)
 t     =  ecpt( 7)
 fmu   =  ecpt(12)
 GO TO 80
 60 ncsid = 9
 ngrids = 4
 matid = necpt(6)
 t     =  ecpt(7)
 fmu   =  ecpt(8)
 GO TO 80
 70 ncsid = 14
 ngrids = 4
 matid = necpt(7)
 t     = 0.0E0
 fmu   = ecpt(11)
 
!  30 COMPUTE PIVOT TRIANGLE AREA
 
!     FIRST SET UP THE POINTERS TO THE CSID OF THE 3 POINTS FROM THE
!     BASE CSID
 
 80    DO  npvt = 1,ngrids
   npt1 = 0
   npt2 = 4
   npt3 = 8
   IF (ngrids == 3 )  GO TO 140
   ichek = 1
!     SELECT 3 POINTS FOR THE PIVOT TRIANGLE OF A QUADRILATERAL
   SELECT CASE ( npvt )
     CASE (    1)
       GO TO 110
     CASE (    2)
       GO TO 140
     CASE (    3)
       GO TO 130
     CASE (    4)
       GO TO 120
   END SELECT
   110   npt3 = 12
   GO TO 140
   120   npt2  =  12
   GO TO 140
   130   npt1 = 12
   
   140 DO  i=1,3
     isub1 = ncsid + npt1 + i
     isub2 = ncsid + npt2 + i
     isub3 = ncsid + npt3 + i
     v1(i) = ecpt(isub3) - ecpt(isub1)
     v2(i) = ecpt(isub3) - ecpt(isub2)
   END DO
   
!     COMPUTE AREA OF QUAD OR TRI USING V1 AND V2
   area= 0.d0
   
   160 v1xv2(1) = v1(2) * v2(3)  -  v1(3) * v2(2)
   v1xv2(2) = v1(3) * v2(1)  -  v1(1) * v2(3)
   v1xv2(3) = v1(1) * v2(2)  -  v1(2) * v2(1)
   
   area = area + DSQRT(v1xv2(1)**2 + v1xv2(2)**2 + v1xv2(3)**2)/2.d0
   
   IF (ngrids == 3)  GO TO 190
   IF( ichek   == 0) THEN
     GO TO   190
   END IF
   
!     COMPUTE AREA OF WHOLE QUAD, FIRST SET UP V1 + V2 THEN TRA TO 600.
   
   170 IF ( narg /= 1 .OR. iflag /= 1 ) GO TO 175
   isub1 = ncsid + npt1 + 1
   isub2 = ncsid + npt2 + 1
   isub3 = ncsid + npt3 + 1
   t = pi23 * ( ecpt(isub1) + ecpt(isub2) + ecpt(isub3) )
   175 npt1 = ncsid
   npt2 = ncsid + 4
   npt3 = ncsid + 8
   npt4 = ncsid +12
   DO  i=1,3
     npt1 = npt1 + 1
     npt2 = npt2 + 1
     npt3 = npt3 + 1
     npt4 = npt4 + 1
     v1(i) = ecpt(npt1) - ecpt(npt3)
     v2(i) = ecpt(npt2) - ecpt(npt4)
   END DO
   ichek = 0
   
   GO TO 160
!     ******************************************************************
!     FINAL COMPUTATION OF TERM AND SHIP OUT OF MATRIX.
   
   190   CONTINUE
   IF( t  == 0.0) THEN
     GO TO   220
   END IF
!     RHO NOT NEEDED IF T = 0
   
   210 inflag = 4
   IF (heat == 1)  GO TO 240
   CALL mat( ecpt(1) )
   rhod = rho
   
   
   220 term =  (fmu + rhod*t)* area/3.d0
   IF (ngrids == 4)  term = term/2.
   i1  = (npvt-1)*3  + 1
   i2  =  i1 + 2
   DO   i = i1,i2
     mass(i) = term
   END DO
   CYCLE
   
!      HEAT FORMULATION
   
   240   CALL hmat(ecpt)
   cpd = cp
   mass(npvt) = (cpd*t)*area/3.d0
   IF (ngrids == 4)  mass(npvt) = mass(npvt) / 2.
   
 END DO
 RETURN
END SUBROUTINE emadtq
