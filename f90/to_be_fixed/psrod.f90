SUBROUTINE psrod
!*****
! THIS ROUTINE COMPUTES STRESSES AND FORCES FOR THE ROD ELEMENT FOR THE
! PLA3 FUNCTIONAL MODULE.
!*****
 
!                        E C P T  F O R  T H E  R O D
 
!                                                                CARD
!                                                 TYPE   TABLE   TYPE
! ECPT( 1)ELEMENT ID.                               I     ECT    CROD
! ECPT( 2)SCALAR INDEX NUMBER FOR GRID POINT A      I     ECT    CROD
! ECPT( 3)SCALAR INDEX NUMBER FOR GRID POINT B      I     ECT    CROD
! ECPT( 4)MATERIAL ID.                              I     EPT    PROD
! ECPT( 5)AREA  (A)                                 R     EPT    PROD
! ECPT( 6)POLAR MOMENT OF INERTIA (J)               R     EPT    PROD
! ECPT( 7) TORSIONAL STRESS COEFF (C)                R    EPT    PROD
! ECPT( 8) NON-STRUCTRAL MASS (MU)                   R    EPT    PROD
! ECPT( 9) COOR. SYS. ID. NO. FOR GRID POINT A       I   BGPDT   GRID
! ECPT(10) X-COORDINATE OF GRID PT. A (IN BASIC COOR)R   BGPDT
! ECPT(11) Y-COORDINATE OF GRID PT. A (IN BASIC COOR)R   BGPDT
! ECPT(12) Z-COORDINATE OF GRID PT. A (IN BASIC COOR)R   BGPDT
! ECPT(13) COOR. SYS. ID. NO. FOR GRID POINT B       I   BGPDT
! ECPT(14) X-COORDINATE OF GRID PT. B (IN BASIC COOR)R   BGPDT
! ECPT(15) Y-COORDINATE OF GRID PT. B (IN BASIC COOR)R   BGPDT
! ECPT(16) Z-COORDINATE OF GRID PT. B (IN BASIC COOR)R   BGPDT
! ECPT(17) ELEMENT TEMPERATURE
! ECPT(18) PREVIOUS STRAIN VALUE, ONCE REMOVED (EPS STAR SUB 0)
! ECPT(19) PREVIOUS STRAIN VALUE  (EPS STAR)
! ECPT(20) PREVIOUSLY COMPUTED VALUE OF MODULUS OF ELASTICITY (ESTAR)
! ECPT(21) PREVIOUSLY COMPUTED TORSIONAL MOMENT (TSTAR)
! ECPT(22) INCREMENTAL DISPLACEMENT VECTOR FOR GRID POINT A
! ECPT(23)                       ...
! ECPT(24)                       ...
! ECPT(25)                       ...
! ECPT(26)                       ...
! ECPT(27)                       ...
! ECPT(28) INCREMENTAL DISPLACEMENT VECTOR FOR GRID POINT B
! ECPT(29)                       ...
! ECPT(30)                       ...
! ECPT(31)                       ...
! ECPT(32)                       ...
! ECPT(33)                       ...
 
 
 
 DIMENSION ecpt(100)          ,iecpt(100)  &
     ,                  xn(3)              ,ua(9)  &
     ,                  ub(9)              ,diff(3)  &
     ,                  ta(9)              ,tb(9)
 
! EST (ECPT) COMMON BLOCK
 
 COMMON   /pla32e/ ecpt
 
! SCRATCH BLOCK FOR VARIABLES LOCAL TO PLA3 ELEMENT ROUTINES.
 
 COMMON   /pla32s/ xl                 ,xn  &
     ,                  ua                 ,ub  &
     ,                  ta                 ,tb ,                  diff
 
! PLA32 COMMUNICATION BLOCK
 
 COMMON   /pla32c/  gamma              ,gammas
 
! OUTPUT BLOCK FOR ELEMENT STRESSES
 
 COMMON   /sout  / jselid             ,sigma  &
     ,                  smsig              ,tau ,                  smtau
 
! INPUT AND OUTPUT BLOCKS FOR SUBROUTINE MAT
 
 COMMON   /matin/ matidc             ,matflg  &
     ,                  temdum             ,plaarg ,                  dum2(2)
 
 
 
 COMMON   /matout/ e sub 0            ,g sub 0  &
     ,                  matdum(5)          ,sigmat  &
     ,                  sigmac             ,sigmas
 
 
 
 EQUIVALENCE (iecpt(1),ecpt(1)) ,(e sub 0,plaans)  &
     ,                  (smsig,mssig)      ,(smtau,mstau)
 
! CALL MAT ROUTINE TO GET MATERIAL PROPERTIES AND STORE IN LOCAL NAMES.
 
 matidc = iecpt(4)
 matflg = 1
 CALL mat (iecpt(1))
 e sub 0 l = e sub 0
 g sub 0 l = g sub 0
 
! SET UP VECTOR ALONG THE ROD, COMPUTE LENGTH AND NORMALIZE
 
 xn(1) = ecpt(10) - ecpt(14)
 xn(2) = ecpt(11) - ecpt(15)
 xn(3) = ecpt(12) - ecpt(16)
 xl    =  xn(1)**2 + xn(2)**2 + xn(3)**2
 xl    = SQRT (xl)
 xn(1) = xn(1) / xl
 xn(2) = xn(2) / xl
 xn(3) = xn(3) / xl
 
! STORE DISPLACEMENT VECTORS IN LOCAL VARIABLES
 
 DO  i = 1,6
   ua(i) = ecpt(i+21)
   ub(i) = ecpt(i+27)
 END DO
 
! TRANSFORM DISPLACEMENT VECTOR TRANSLATIONAL COMPONENTS IF NECESSARY
 
 ibasea = 0
 IF (iecpt(9) == 0) GO TO 20
 ibasea = 6
 CALL transs (iecpt(9),ta)
 CALL gmmats (ta,3,3,0, ua(1),3,1,0, ua(7))
 20 ibaseb = 0
 IF (iecpt(13) == 0) GO TO 30
 ibaseb = 6
 CALL transs (iecpt(13),tb)
 CALL gmmats (tb,3,3,0, ub(1),3,1,0, ub(7))
 
! FORM DIFFERENCE VECTOR, DOT PRODUCT AND INCREMENT OF STRAIN
 
 30 diff(1) = ua(ibasea+1) - ub(ibaseb+1)
 diff(2) = ua(ibasea+2) - ub(ibaseb+2)
 diff(3) = ua(ibasea+3) - ub(ibaseb+3)
 CALL gmmats (xn,3,1,1, diff,3,1,0, term)
 deps1 = term / xl
 epsin2 = ecpt(19)
 epsin1 = ecpt(18)
 deps2 = epsin2 - epsin1
 
! COMPUTE EPS1 AND EPS2 AND FETCH VIA MAT STRESSES SIGMA1 AND SIGMA2
 
 eps1 = epsin2 + deps1
 eps2 = epsin2  +  (deps1+gammas**2*deps2)*(gamma+1.0)/(gammas+1.0)  &
     + gammas*(deps1-gammas*deps2)*(gamma+1.0)**2 / (gammas+1.0)
 matflg = 6
 plaarg = eps1
 CALL mat (iecpt(1))
 sigma1 = plaans
 plaarg = eps2
 CALL mat (iecpt(1))
 sigma2 = plaans
 IF (eps1 == eps2) GO TO 42
 e = (sigma2 - sigma1) / (eps2 - eps1)
 GO TO 44
 42 e = ecpt(20)
 44 g = ecpt(20) * g sub 0 l / e sub 0 l
 
! COMPUTE STRESSES
 
 iselid = iecpt(1)
 sigma  = sigma1
 p      =  ecpt(5) * sigma1
 
! TRANSFORM DISPLACEMENT VECTOR ROTATIONAL DISPLACEMENTS IF NECESSARY.
 
 ibasea = 3
 IF (iecpt(9) == 0) GO TO 60
 CALL gmmats (ta,3,3,0, ua(4),3,1,0, ua(7))
 ibasea = 6
 60 ibaseb = 3
 IF (iecpt(13) == 0) GO TO 70
 ibaseb = 6
 CALL gmmats (tb,3,3,0, ub(4),3,1,0, ub(7))
 70 diff(1) = ua(ibasea+1) - ub(ibaseb+1)
 diff(2) = ua(ibasea+2) - ub(ibaseb+2)
 diff(3) = ua(ibasea+3) - ub(ibaseb+3)
 CALL gmmats (xn,3,1,1, diff,3,1,0, term)
 t = ecpt(6) * g * term / xl  +  ecpt(21)
 IF (ecpt(6) == 0.0) GO TO 80
 tau = ecpt(7) * t / ecpt(6)
 GO TO 90
 80 tau = 0.0
 
! COMPUTE MARGIN OF SAFETY IN EXTENSION
 
 90 IF(sigma <= 0.0)GO TO 101
 IF(sigmat <= 0.0)GO TO 102
 smsig=sigmat/sigma-1.0
 GO TO 180
 101 IF(sigma /= 0.0) GO TO 103
 GO TO 102
 103 sigmac=-ABS(sigmac)
 smsig=sigmac/sigma - 1.0
 GO TO 180
 102 mssig=1
 
!     COMPUTE MARGIN OF SAFETY IN TORSION
 
 180 IF(sigmas <= 0.0) GO TO 190
 IF(tau == 0.0)GO TO 190
 smtau= sigmas/ABS(tau) - 1.0
 GO TO 210
 190 mstau=1
 210 jselid = iecpt(1)
 
! UPDATE EST (ECPT) ENTRY
 
 ecpt(18) = ecpt(19)
 ecpt(19) = eps1
 ecpt(20) = e
 ecpt(21) = t
 RETURN
END SUBROUTINE psrod
