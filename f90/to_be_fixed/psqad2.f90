SUBROUTINE psqad2
!  THIS SUBROUTINE IS THE DRIVER FOR THE QUAD2 CALCULATIONS IN
!  PLA3
 
!     ECPT FOR QUAD2
 
!  1  EL.ID
!  2  GRID A
!  3  GRID B
!  4  GRID C
!  5  GRID D
!  6  THETA
!  7  MAT ID
!  8  T
!  9  MS MASS
! 10  CSID 1
! 11  X1
! 12  Y1
! 13  Z1
! 14  CSID 2
! 15  X2
! 16  Y2
! 17  Z2
! 18  CSID 3
! 19  X3
! 20  Y3
! 21  Z3
! 22  CSID 4
! 23  X4
! 24  Y4
! 25  Z4
! 26  TEMP
! 27  EPS0
! 28  EPSS
! 29  ESTAR
! 30  SIGXS
! 31  SIGYS
! 32  SIGXXS
! 33  MXS
! 34  MYS
! 35  MXYS
! 36  VXS
! 37  VYS
! 38  U(A) (6X1)
! 44  U(B) (6X1)
! 50  U(C) (6X1)
! 56  U(D) (6X1)
 
!     ******************************************************************
 REAL :: nu
 
 DIMENSION necpt(26), necpts(26)
 
 COMMON /pla32e/ ecpt(26),eps0,epss,estar,sigxs,sigys,sigxys,  &
     forvec(5), ui(24),  dummy(39)
 COMMON /pla3es/ ecptsa(100),ph1out(200)
 COMMON /pla3uv/  ivec, z(24)
 
! SCRATCH BLOCK  325 CELLS
 
 COMMON /pla32s/s(3),dum(297),tau0 ,tau1 ,tau2 ,f,sx,sy,deps,depss,  &
     eps1,eps2, dum1,idum2,idum3(3,3) ,              extra(4)
 COMMON /matin/ matid,inflag,eltemp,plaarg,sinth,costh
 COMMON /matout/ g11,g12,g13,g22,g23,g33
 COMMON /pla32c/ gamma, gammas, ipass
 COMMON /plagp/  gp(9), midgp   , elid
 
 EQUIVALENCE (necpt(7),matid1) , (ecpt(1),necpt(1)) , (g11,plaans),  &
     (g13,nu)   , (g11,esub0) ,  (necpts(1),ecptsa(1))  &
     ,   (g12,nirof)
 
! SETUP GP MATRIX FOR PLAMAT
 
 elid = ecpt(1)
 midgp = matid1
 DO  i=1,9
   gp(i)=0.0
 END DO
 tau0  = SQRT(sigxs**2 - sigxs*sigys + sigys**2 + 3.0*sigxys**2)
 IF(estar == 0.0) GO TO 50
 IF(ipass /= 1  ) GO TO 20
 matid = matid1
 costh = 1.0
 sinth = 0.0E0
 inflag= 2
 
 CALL mat(ecpt(1))
 
 gp(1) = g11
 gp(2) = g12
 gp(3) = g13
 gp(4) = g12
 gp(5) = g22
 gp(6) = g23
 gp(7) = g13
 gp(8) = g23
 gp(9) = g33
 GO TO 50
 20 IF(tau0  == 0.0) GO TO 50
 matid = matid1
 inflag = 1
 
 CALL mat(ecpt(1))
 
 f =   9.0*(esub0 - estar) / (4.0 * tau0**2 * estar)
 sx = (2.0*sigxs - sigys)/ 3.0
 sy = (2.0*sigys - sigxs)/ 3.0
 gp(1) = (1.0+sx**2*f) / esub0
 gp(2) = (-nu+sx*sy*f) / esub0
 gp(3) = (2.0*sigxys*sx*f) / esub0
 gp(4) = gp(2)
 gp(5) = (1.0+sy**2*f) / esub0
 gp(6) = (2.0*sigxys*sy*f) / esub0
 gp(7) = gp(3)
 gp(8) = gp(6)
 gp(9) = (2.0*(1.0+nu) + 4.0*f*sigxys**2) / esub0
 
!     NO NEED TO COMPUTE DETERMINANT SINCE IT IS NOT USED SUBSEQUENTLY.
 idum2 = -1
 CALL invers(3,gp,3,0,0,dum1,idum2,idum3)
 
! CHECK SINGULARITY
 
 IF( idum2 == 2) CALL mesage(-30,38,ecpt(1))
 
! CALCULATE PHASE I STRESSES
 
 50 DO  i = 1,32
   ecptsa(i) = ecpt(i)
 END DO
 necpts(2) = 1
 necpts(3) = 7
 necpts(4) = 13
 necpts(5) = 19
 
 CALL pstq1(4)
 
 
! CALCULATE PHASE II STRESSES
 
 ivec = 1
 DO  i = 1,24
   z(i) = ui(i)
 END DO
 DO  i = 1,200
   ecptsa(i) = ph1out(i)
 END DO
 s(1) = sigxs
 s(2) = sigys
 s(3) = sigxys
 i201 = 201
 ecptsa(i201) = ecpt(1)
 DO  i=1,5
   ecptsa(i+201) = forvec(i)
 END DO
 
 CALL pstq2(4)
 
 
!  UPDATE ECPT FOR STRESSES
 
 sigxs = s(1)
 sigys = s(2)
 sigxys = s(3)
 
!     NEW FORCES ARE IN  /PLA3ES/ AT LOCATIONS 202-206
 
 DO  i=1,5
   forvec(i)   = ecptsa(i+201)
 END DO
 tau1  = SQRT(sigxs**2 - sigxs*sigys + sigys**2 + 3.0*sigxys**2)
 matid= matid1
 inflag = 8
 plaarg = tau1
 
 CALL mat(ecpt(1))
 
 
! TEST FOR TAU 1 OUTSIDE THE RANGE OF FUNCTION
 
 IF ( nirof . EQ. 1 ) GO TO 80
 
! RETURNS EPS SUB 1 GIVEN TAU1
 
 eps1 = plaans
 deps = eps1 - epss
 depss= epss - eps0
 eps2 = eps1 + gamma * deps
 inflag=6
 plaarg = eps2
 
 CALL mat(ecpt(1))
 
! RETURNS  TAU2 GIVEN EPS2
 
 tau2  = plaans
 estar = 0.0
 IF( (eps2 - eps1) /= 0.0) estar = (tau2 - tau1) / (eps2-eps1)
 eps0  = epss
 epss  = eps1
 RETURN
 
 80 estar = 0.0
 RETURN
END SUBROUTINE psqad2
