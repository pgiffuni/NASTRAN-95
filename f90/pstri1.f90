SUBROUTINE pstri1
!  THIS ROUTINE CALCULATES GP,SET-S UP THE ECPT AND UPDATES THE ECPT
!  FOR THE TRIA1 ELEMENTS
 
!     ECPT FOR TRIA1
!     EL.ID                                       ECPT( 1)
!     GRID A                                      ECPT( 2)
!     GRID B                                      ECPT( 3)
!     GRID C                                      ECPT( 4)
!     THETA                                       ECPT( 5)
!     MATID1                                      ECPT( 6)
!     T1                                          ECPT( 7)
!     MATID2                                      ECPT( 8)
!     I                                           ECPT( 9)
!     MATID3                                      ECPT(10)
!     T2                                          ECPT(11)
!     NS MASS                                     ECPT(12)
!     Z1                                          ECPT(13)
!     Z2                                          ECPT(14)
!     CSID 1                                      ECPT(15)
!     X1                                          ECPT(16)
!     Y1                                          ECPT(17)
!     Z1                                          ECPT(18)
!     CSID 2                                      ECPT(19)
!     X2                                          ECPT(20)
!     Y2                                          ECPT(21)
!     Z2                                          ECPT(22)
!     CSID3                                       ECPT(23)
!     X3                                          ECPT(24)
!     Y3                                          ECPT(25)
!     Z3                                          ECPT(26)
!     TEMP                                        ECPT(27)
!     EPS SUB 0      (PREVIOUS  STRAIN)           ECPT(28)
!     EPS SUB STAR   (LAST STRAIN)                ECPT(29)
!     MODULUS OF ELASTICITY                       ECPT(30)
!     SIGMA X        STRESS                       ECPT(31)
!     SIGMA Y        STRESS                       ECPT(32)
!     SIGMA XY       STRESS                       ECPT(33)
!     M X STAR       FORCE                        ECPT(34)
!     M Y STAR       FORCE                        ECPT(35)
!     M XX STAR      FORCE                        ECPT(36)
!     V X  STAR      FORCE                        ECPT(37)
!     V Y  STAR      FORCE                        ECPT(38)
!     U A       (6X1 DISPLACEMENT VECTOR)         ECPT(39)
!     U B       (6X1 DISPLACEMENT VECTOR)         ECPT(45)
!     U C       (6X1 DISPLACEMENT VECTOR)         ECPT(51)
 
!     ******************************************************************
 
 REAL :: nu
 
 DIMENSION necpt(27),necpts(27)
 COMMON /pla32e/ ecpt(27), eps0,epss,estar,sigxs,sigys,sigxys,  &
     forvec(5),ui(18),dummy(44)
 COMMON /pla3es/ ecptsa(100), ph1out(200)
 COMMON /pla3uv/ ivec, z(24)
 
! SCRATCH BLOCK  325 CELLS
 
 COMMON /pla32s/s(3),dum(297),tau0 ,tau1 ,tau2 ,f,sx,sy,deps,depss,  &
     eps1,eps2, dum1,idum2,idum3(3,3) ,              extra(4)
 COMMON /matin/ matid,inflag,eltemp,plaarg,sinth,costh
 COMMON /matout/ g11,g12,g13,g22,g23,g33
 COMMON /pla32c/ gamma, gammas, ipass
 COMMON /plagp/  gp(9),  midgp  , elid
 
 EQUIVALENCE (necpt(6),matid1) , (ecpt(1),necpt(1)) , (g11,plaans)  &
     ,           (g13,nu)   , (g11,esub0) ,  (necpts(1),ecptsa(1))  &
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
 
 CALL pstq1 (1)
 
 
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
 
 CALL pstq2 (3)
 
 
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
 
 IF ( nirof == 1 ) GO TO 80
 
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
END SUBROUTINE pstri1
