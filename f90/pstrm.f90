SUBROUTINE pstrm
     
!     THIS SUBROUTINE IS THE DRIVER FOR THE  TRI-MEMBRANE CALCULATIONS
!     IN PLA3
 
!     ECPT LIST
!                                                      IN
!                                                      THIS
!       ECPT       DESCRIPTION                         ROUTINE   TYPE
!     ******************************************************************
!       ECPT( 1) = ELEMENT ID                          NECPT(1)  INTEGER
!       ECPT( 2) = GRID POINT A                        NGRID(1)  INTEGER
!       ECPT( 3) = GRID POINT B                        NGRID(2)  INTEGER
!       ECPT( 4) = GRID POINT C                        NGRID(3)  INTEGER
!       ECPT( 5) = THETA = ANGLE OF MATERIAL           ANGLE     REAL
!       ECPT( 6) = MATERIAL ID                         MATID     INTEGER
!       ECPT( 7) = T                                   T         REAL
!       ECPT( 8) = NON-STRUCTURAL MASS                 FMU       REAL
!       ECPT( 9) = COORD. SYSTEM ID 1                  NECPT(9)  INTEGER
!       ECPT(10) = X1                                  X1        REAL
!       ECPT(11) = Y1                                  Y1        REAL
!       ECPT(12) = Z1                                  Z1        REAL
!       ECPT(13) = COORD. SYSTEM ID 2                  NECPT(13) INTEGER
!       ECPT(14) = X2                                  X2        REAL
!       ECPT(15) = Y2                                  Y2        REAL
!       ECPT(16) = Z2                                  Z2        REAL
!       ECPT(17) = COORD. SYSTEM ID 3                  NECPT(17) INTEGER
!       ECPT(18) = X3                                  X3        REAL
!       ECPT(19) = Y3                                  Y3        REAL
!       ECPT(20) = Z3                                  Z3        REAL
!       ECPT(21) = ELEMENT TEMPERATURE                 ELTEMP    REAL
!       ECPT(22) = STRAIN (MINUS ONE)                  EPS0      REAL
!       ECPT(23) = STRAIN (PRESENT)                    EPSS      REAL
!       ECPT(24) = MODULUS OF ELASTICITY               ESTAR     REAL
!       ECPT(25) = STRESS SUB X                        SIGXS     REAL
!       ECPT(26) = STRESS SUB Y                        SIGYS     REAL
!       ECPT(27) = STRESS SUB XY                       SIGXYS    REAL
!       ECPT(28) = DISPLACEMENT VECTOR   A1            UI(1)     REAL
!       ECPT(29) = DISPLACEMENT VECTOR   A2            UI(2)     REAL
!       ECPT(30) = DISPLACEMENT VECTOR   A3            UI(3)     REAL
!       ECPT(31) = DISPLACEMENT VECTOR   B1            UI(4)     REAL
!       ECPT(32) = DISPLACEMENT VECTOR   B2            UI(5)     REAL
!       ECPT(33) = DISPLACEMENT VECTOR   B3            UI(6)     REAL
!       ECPT(34) = DISPLACEMENT VECTOR   C1            UI(7)     REAL
!       ECPT(35) = DISPLACEMENT VECTOR   C2            UI(8)     REAL
!       ECPT(36) = DISPLACEMENT VECTOR   C3            UI(9)     REAL
 
!     ******************************************************************
 
 REAL :: nu
 DIMENSION       necpt(21),necpts(21)
 COMMON /pla32e/ ecpt(21),eps0,epss,estar,sigxs,sigys,sigxys, ui(9),dummy(64)
 COMMON /pla3es/ ecptsa(100),ph1out(200)
 COMMON /pla3uv/ ivec,z(24)
 
!     SCRATCH BLOCK  325 CELLS
 
 COMMON /pla32s/ s(3),dum(297),tau0,tau1,tau2,f,sx,sy,deps,  &
     depss,eps1,eps2,dum1,idum2,idum3(3,3),extra(4)
 COMMON /matin / matid,inflag,eltemp,plaarg,sinth,costh
 COMMON /matout/ g11,g12,g13,g22,g23,g33
 COMMON /pla32c/ gamma,gammas,ipass
 COMMON /plagp / gp(9),midgp,elid
 EQUIVALENCE     (necpt(6),matid1),(ecpt(1),necpt(1)),(g11,plaans),  &
     (g13,nu),(g11,esub0),(necpts(1),ecptsa(1)), (g12,nirof)
 
!     SETUP GP MATRIX FOR PLAMAT
 
 elid  = ecpt(1)
 midgp = matid1
 DO  i = 1,9
   gp(i) = 0.0
 END DO
 tau0  = SQRT(sigxs**2 - sigxs*sigys + sigys**2 + 3.0*sigxys**2)
 IF (estar == 0.) GO TO 50
 IF (ipass /=  1) GO TO 20
 matid = matid1
 costh = 1.0
 sinth = 0.0
 inflag = 2
 
 CALL mat (ecpt(1))
 
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
 20 IF (tau0 == 0.0) GO TO 50
 matid = matid1
 inflag = 1
 
 CALL mat (ecpt(1))
 
 f     =  9.0*(esub0-estar)/(4.0*tau0**2*estar)
 sx    = (2.0*sigxs-sigys)/3.0
 sy    = (2.0*sigys-sigxs)/3.0
 gp(1) = (1.0+sx**2*f) /esub0
 gp(2) = (-nu+sx*sy*f) /esub0
 gp(3) = (2.0*sigxys*sx*f)/esub0
 gp(4) = gp(2)
 gp(5) = (1.0+sy**2*f)/esub0
 gp(6) = (2.0*sigxys*sy*f)/esub0
 gp(7) = gp(3)
 gp(8) = gp(6)
 gp(9) = (2.0*(1.0+nu) + 4.0*f*sigxys**2)/esub0
 
!     NO NEED TO COMPUTE DETERMINANT SINCE IT IS NOT USED SUBSEQUNTLY.
 
 idum2 = -1
 CALL invers (3,gp,3,0,0,dum1,idum2,idum3)
 
!     CHECK SINGULARITY
 
 IF (idum2 == 2) CALL mesage (-30,38,ecpt(1))
 
!     CALCULATE PHASE I STRESSES
 
 50 DO  i = 1,32
   ecptsa(i) = ecpt(i)
 END DO
 necpts(2) = 1
 necpts(3) = 4
 necpts(4) = 7
 
 CALL pstrm1 (0)
 
!     CALCULATE PHASE II STRESSES
 
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
 
 CALL pstrq2 (1)
 
!     UPDATE ECPT FOR STRESSES
 
 sigxs  = s(1)
 sigys  = s(2)
 sigxys = s(3)
 tau1   = SQRT(sigxs**2 - sigxs*sigys + sigys**2 + 3.0*sigxys**2)
 matid  = matid1
 inflag = 8
 plaarg = tau1
 
 CALL mat (ecpt(1))
 
!     TEST FOR TAU 1 OUTSIDE THE RANGE OF FUNCTION
 
 IF (nirof == 1) GO TO 80
 
!     RETURNS EPS SUB 1 GIVEN TAU1
 
 eps1 = plaans
 deps = eps1 - epss
 depss= epss - eps0
 eps2 = eps1 + gamma * deps
 inflag = 6
 plaarg = eps2
 
 CALL mat (ecpt(1))
 
!     RETURNS  TAU2 GIVEN EPS2
 
 tau2  = plaans
 estar = 0.0
 IF (eps2-eps1 /= 0.0) estar = (tau2-tau1)/(eps2-eps1)
 eps0  = epss
 epss  = eps1
 RETURN
 
 80 estar = 0.0
 RETURN
END SUBROUTINE pstrm
