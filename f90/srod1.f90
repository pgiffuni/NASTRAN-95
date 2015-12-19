SUBROUTINE srod1
!*****
! THIS ROUTINE IS PHASE I OF STRESS DATA RECOVERY FOR THE ROD.
!*****
 
 
 
 DIMENSION          iecpt(13)
 
 
 
 
! INPUT AND OUTPUT BLOCK
 
 COMMON   /sdr2x5/ ecpt(17)           ,dummy1(83),  &
     ielid              ,isilno(2)  &
     ,                  sat(3)             ,sbt(3)  &
     ,                  sar(3)             ,sbr(3)  &
     ,                  st                 ,sdelta  &
     ,                  area               ,fjovrc  &
     ,                  t_subc_0           ,sigmat  &
     ,                  sigmac             ,sigmas  &
     ,                  sigvec(77)         ,forvec(25)
 
! SCRATCH BLOCK
 
 COMMON   /sdr2x6/ xn(6)              ,ti(9)  &
     ,                  xl                 ,eoverl ,                  ibase
 
! INPUT AND OUTPUT BLOCKS FOR SUBROUTINE MAT
 
 COMMON   /matin/ matidc             ,matflg  &
     ,                  eltemp             ,stress  &
     ,                  sinth              ,costh
 
 
 
 COMMON   /matout/ e                  ,g  &
     ,                  nu                 ,rho  &
     ,                  alpha              ,t_sub_0  &
     ,                  gsube              ,sigt  &
     ,                  sigc               ,sigs
 
 
 
 EQUIVALENCE        (iecpt(1),ecpt(1))
 
! CALL MAT TO GET MATERIAL PROPERTIES
 
 matidc = iecpt(4)
 matflg = 1
 eltemp = ecpt(17)
 CALL mat (iecpt(1))
 
! SET UP VECTOR ALONG THE ROD, COMPUTE LENGTH AND NORMALIZE
 
 xn(1) = ecpt(10) - ecpt(14)
 xn(2) = ecpt(11) - ecpt(15)
 xn(3) = ecpt(12) - ecpt(16)
 xl =  xn(1)**2  +  xn(2)**2  +  xn(3)**2
 xl =  SQRT(xl)
 xn(1) = xn(1) / xl
 xn(2) = xn(2) / xl
 xn(3) = xn(3) / xl
 eoverl = e / xl
 gcovrl = g * ecpt(6) / xl
 ibase = 0
 
! TRANSFORM XN VECTOR IF POINT A IS NOT IN BASIC COORDINATES.
 
 IF (iecpt(9) == 0) GO TO 10
 ibase = 3
 CALL transs (iecpt(9),ti)
 CALL gmmats (xn(1),3,1,1, ti(1),3,3,0, xn(4) )
 10 sat(1) = xn(ibase +1) * eoverl
 sat(2) = xn(ibase +2) * eoverl
 sat(3) = xn(ibase +3) * eoverl
 sar(1) = xn(ibase +1) * gcovrl
 sar(2) = xn(ibase +2) * gcovrl
 sar(3) = xn(ibase +3) * gcovrl
 
! TRANSFORM XN VECTOR IF POINT B IS NOT IN BASIC COORDINATES.
 
 ibase = 0
 IF (iecpt(13) == 0) GO TO 20
 ibase = 3
 CALL transs (iecpt(13),ti)
 CALL gmmats (xn(1),3,1,1, ti(1),3,3,0, xn(4) )
 20 sbt(1) = - xn(ibase+1) * eoverl
 sbt(2) = - xn(ibase+2) * eoverl
 sbt(3) = - xn(ibase+3) * eoverl
 sbr(1) = - xn(ibase+1) * gcovrl
 sbr(2) = - xn(ibase+2) * gcovrl
 sbr(3) = - xn(ibase+3) * gcovrl
 
! FILL REMAINDER OF OUTPUT BLOCK
 
 st     = - alpha * e
 sdelta = - eoverl
 area   =   ecpt(5)
 IF(ecpt(6) == 0.0) THEN
   GO TO    40
 END IF
 30 fjovrc = ecpt(7) / ecpt(6)
 GO TO 50
 40 fjovrc = 0.0
 50 tsubc0 = tsub0
 sigmat = sigt
 sigmac = sigc
 sigmas = sigs
 ielid = iecpt(1)
 isilno(1) = iecpt(2)
 isilno(2) = iecpt(3)
 RETURN
END SUBROUTINE srod1
