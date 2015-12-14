SUBROUTINE pktri2
!  THIS SUBROUTINE IS THE DRIVER FOR THE  TRIA2 CALCULATIONS IN
!  PLA4
 
!     ECPT FOR  TRIA2
 
!  1  EL.ID
!  2  GRID A
!  3  GRID B
!  4  GRID C
!  5  THETA
!  6  MAT ID
!  7  T
!  8  MS MASS
!  9  CSID 1
! 10  X1
! 11  Y1
! 12  Z1
! 13  CSID 2
! 14  X2
! 15  Y2
! 16  Z2
! 17  CSID 3
! 18  X3
! 19  Y3
! 20  Z3
! 21  TEMP
! 22  EPS0
! 23  EPSS
! 24  ESTAR
! 25  SIGXS
! 26  SIGYS
! 27  SIGXYS
! 28  U(A) (3X1)
! 31  U(B) (3X1)
! 34  U(C) (3X1)
 
!     ******************************************************************
 
 LOGICAL :: istiff
 
 REAL :: nu
 
 DIMENSION necpt(21), necpts(21)
 
 COMMON /pla42e/ ecpt(21),eps0,epss,estar,sigxs,sigys,sigxys,  &
     ui(09),  dummy(64)
 COMMON /pla4es/ ecptsa(100), ph1out(200)
 COMMON /pla4uv/ ivec, z(24)
 
! SCRATCH BLOCK  325 CELLS
 
 COMMON /pla42s/s(3),dum(297),tau0 ,tau1 ,tau2 ,f,sx,sy,deps,depss,  &
     eps1,eps2, dum1,idum2,idum3(3,3) ,              extra(4)
 COMMON /matin/ matid,inflag,eltemp,plaarg,sinth,costh
 COMMON /matout/ g11,g12,g13,g22,g23,g33
 COMMON /pla42c/  npvt, gamma, gammas, ipass  &
     ,                  dumcl(145)         ,nogo
 COMMON /plagp/  gp(9) , midgp  , elid
 
 EQUIVALENCE (necpt(6),matid1) , (ecpt(1),necpt(1)) , (g11,plaans)  &
     ,           (g13,nu)   , (g11,esub0) ,  (necpts(1),ecptsa(1))  &
     ,   (g12,nirof)
 
! SETUP GP MATRIX FOR PLAMAT
 
 istiff = .false.
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
 120 matid = matid1
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
 
 IF (idum2 == 2) GO TO 150
 
 50 IF(istiff) GO TO 130
 istiff = .true.
 
! CALCULATE PHASE I STRESSES
 
 DO  i = 1,32
   ecptsa(i) = ecpt(i)
 END DO
 necpts(2) = 1
 necpts(3) = 4
 necpts(4) = 7
 
 CALL pktq1(2)
 
 
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
 
 CALL pktq2(3)
 
!  UPDATE ECPT FOR STRESSES
 
 sigxs = s(1)
 sigys = s(2)
 sigxys = s(3)
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
 GO TO 100
 80 estar = 0.0
 
!  SETUP STIFFNESS CALCULATIONS FOR GP
 
 100 DO  i = 1,9
   gp(i) = 0.0
 END DO
 tau0  = SQRT(sigxs**2 - sigxs*sigys + sigys**2 + 3.0*sigxys**2)
 IF( estar /= 0.0 .AND. tau0 /= 0.0) GO TO 120
 
!  SETUP CALL TO ELEMENT STIFFNESS ROUTINE IT WILL ALSO INSERT
 
 130 DO  i = 1,32
   ecptsa(i) = ecpt(i)
 END DO
 CALL pktrqd(2)
 RETURN
 150 CALL mesage(30,38,ecpt(1))
 
!  SET FLAG FOR FATAL ERROR WHILE ALLOWING ERROR MESSAGES TO ACCUMULATE
 
 nogo=1
 RETURN
END SUBROUTINE pktri2
