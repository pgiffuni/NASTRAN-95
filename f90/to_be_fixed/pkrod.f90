SUBROUTINE pkrod
     
!     THIS ROUTINE COMPUTES THE TWO 6 X 6 MATRICES  K(NPVT,NPVT) AND
!     K(NPVT,J) FOR A ROD HAVING END POINTS NUMBERED NPVT AND J.
 
!     ECPT FOR THE ROD
!     ================                                              CARD
!                                                      TYPE  TABLE  TYPE
!     ECPT( 1)ELEMENT ID.                                I    ECT   CROD
!     ECPT( 2)SCALAR INDEX NUMBER FOR GRID POINT A       I    ECT   CROD
!     ECPT( 3)SCALAR INDEX NUMBER FOR GRID POINT B       I    ECT   CROD
!     ECPT( 4)MATERIAL ID.                               I    EPT   PROD
!     ECPT( 5)AREA  (A)                                  R    EPT   PROD
!     ECPT( 6)POLAR MOMENT OF INERTIA (J)                R    EPT   PROD
!     ECPT( 7) TORSIONAL STRESS COEFF (C)                R    EPT   PROD
!     ECPT( 8) NON-STRUCTRAL MASS (MU)                   R    EPT   PROD
!     ECPT( 9) COOR. SYS. ID. NO. FOR GRID POINT A       I   BGPDT  GRID
!     ECPT(10) X-COORDINATE OF GRID PT. A (IN BASIC COOR)R   BGPDT
!     ECPT(11) Y-COORDINATE OF GRID PT. A (IN BASIC COOR)R   BGPDT
!     ECPT(12) Z-COORDINATE OF GRID PT. A (IN BASIC COOR)R   BGPDT
!     ECPT(13) COOR. SYS. ID. NO. FOR GRID POINT B       I   BGPDT
!     ECPT(14) X-COORDINATE OF GRID PT. B (IN BASIC COOR)R   BGPDT
!     ECPT(15) Y-COORDINATE OF GRID PT. B (IN BASIC COOR)R   BGPDT
!     ECPT(16) Z-COORDINATE OF GRID PT. B (IN BASIC COOR)R   BGPDT
!     ECPT(17) ELEMENT TEMPERATURE
!     ECPT(18) PREVIOUS STRAIN VALUE, ONCE REMOVED (EPSIN1)
!     ECPT(19) PREVIOUS STRAIN VALUE (EPSIN2)
!     ECPT(20) PREVIOUSLY COMPUTED VALUE OF MODULUS OF ELASTICITY, ESTAR
!     ECPT(21) DISPLACEMENT COORDINATES FOR GRID POINT A
!     ECPT(22)                   . . .
!     ECPT(23)                   . . .
!     ECPT(24) DISPLACEMENT COORDINATES FOR GRID POINT B
!     ECPT(25)                   . . .
!     ECPT(26)                   . . .
 
 DOUBLE PRECISION :: d(18),x,y,z,xl,xn(3),ua(6),ub(6),ta(9),tb(9),e,g,  &
     diff(3),dpterm,epsin1,epsin2,deps1,deps2,eps1,  &
     eps2,gamma,gammas,sigma1,sigma2,dscl,dscr,ke(36)
 DIMENSION        iecpt(200)
 
!     PLA42 PARAMETERS COMMUNICATION BLOCK
 COMMON /pla42c/  npvt,g NEW,g OLD,dumcl(146),nogo
 
!     ECPT COMMON BLOCK
 COMMON /pla42e/  ecpt(100)
 
!     PLA42 LOCAL VARIABLE (SCRATCH) BLOCK
 COMMON /pla42d/  d,x,y,z,xl,xn,ua,ub,ta,tb,diff,dpterm,epsin1,  &
     epsin2,deps1,deps2,eps1,eps2,gamma,gammas, sigma1,sigma2,dscl,dscr,e,g,ke
 
!     INPUT AND OUTPUT BLOCKS FOR SUBROUTINE MAT
 COMMON /matin /  matidc,matflg,temdum,plaarg,matdum(2)
 COMMON /matout/  e sub 0,g sub 0,dummat(18)
 EQUIVALENCE      (iecpt(1),ecpt(1)) ,(plaans,esub0)
 
!     BEGIN EXECUTION
 
 ind = 0
 IF (iecpt(2) == npvt) GO TO 10
 IF (iecpt(3) /= npvt) CALL mesage (-30,34,iecpt(1))
 ind = 1
 itemp = iecpt(2)
 iecpt(2) = iecpt(3)
 iecpt(3) = itemp
 ka  = 13
 kb  =  9
 idispa = 23
 idispb = 20
 GO TO 20
 10 ka  =  9
 kb  = 13
 idispa = 20
 idispb = 23
 
!     AT THIS POINT KA POINTS TO THE COOR. SYS. ID. OF THE PIVOT GRID
!     POINT. SIMILARLY FOR KB AND THE NON-PIVOT GRID POINT.
!     NOW COMPUTE THE LENGTH OF THE ROD.
 
!     WE STORE THE COORDINATES IN THE D ARRAY SO THAT ALL ARITHMETIC
!     WILL BE DOUBLE PRECISION
 
 20 d(1) = ecpt(ka+1)
 d(2) = ecpt(ka+2)
 d(3) = ecpt(ka+3)
 d(4) = ecpt(kb+1)
 d(5) = ecpt(kb+2)
 d(6) = ecpt(kb+3)
 x    = d(1) - d(4)
 y    = d(2) - d(5)
 z    = d(3) - d(6)
 xl = DSQRT(x**2 + y**2 + z**2)
 IF (xl /= 0.0D0) GO TO 25
 CALL mesage (30,26,iecpt(1))
 
!     SET FLAG FOR FATAL ERROR WHILE ALLOWING ERROR MESSAGES TO
!     ACCUMULATE
 
 nogo = 1
 RETURN
 
!     CALCULATE A NORMALIZED DIRECTION VECTOR IN BASIC COORDINATES.
 
 25 xn(1) = x/xl
 xn(2) = y/xl
 xn(3) = z/xl
 
!     STORE DISPLACEMENT VECTORS IN DOUBLE PRECISION LOCATIONS
 
 ua(1) = ecpt(idispa+1)
 ua(2) = ecpt(idispa+2)
 ua(3) = ecpt(idispa+3)
 ub(1) = ecpt(idispb+1)
 ub(2) = ecpt(idispb+2)
 ub(3) = ecpt(idispb+3)
 
 
!     COMPUTE THE DIFFERENCE VECTOR DIFF =  T  * U   -  T  * U
!                                            A    A      B    B
 
 ibasea = 0
 IF (iecpt(ka) == 0) GO TO 30
 CALL transd (ecpt(ka),ta)
 ibasea = 3
 CALL gmmatd (ta,3,3,0, ua(1),3,1,0, ua(4))
 30 ibaseb = 0
 IF (iecpt(kb) == 0) GO TO 40
 CALL transd (ecpt(kb),tb)
 ibaseb = 3
 CALL gmmatd (tb,3,3,0, ub(1),3,1,0, ub(4))
 40 diff(1) = ua(ibasea+1) - ub(ibaseb+1)
 diff(2) = ua(ibasea+2) - ub(ibaseb+2)
 diff(3) = ua(ibasea+3) - ub(ibaseb+3)
 
!     COMPUTE DOT PRODUCT XN . DIFF
 
 CALL gmmatd (xn,3,1,1, diff,3,1,0, dpterm)
 
!     COMPUTE INCREMENT OF STRAIN
 
 deps1  = dpterm/xl
 epsin1 = ecpt(18)
 epsin2 = ecpt(19)
 deps2  = epsin2 - epsin1
 
!     COMPUTE CURRENT STRAIN AND ESTIMATED NEXT STRAIN
 
 eps1   = epsin2 + deps1
 gamma  = g NEW
 gammas = g OLD
 eps2   = eps1 + gamma*deps1
 
!     CALL MAT ROUTINE TWICE TO GET SIGMA1 AND SIGMA2 AS A FUNCTION OF
!     EPS1 AND EPS2
 
 matidc = iecpt(4)
 matflg = 6
 plaarg = eps1
 CALL mat (iecpt(1))
 sigma1 = plaans
 plaarg = eps2
 CALL mat (iecpt(1))
 sigma2 = plaans
 
!     ON THE FIRST PASS, I.E. WHEN ECPT(19) = 0.0, SIGMA1 = E  * EPS1
!                                                            0
 
 IF (ecpt(19) /= 0.0) GO TO 41
 matflg = 1
 CALL mat (iecpt(1))
 d(2)   = e sub 0
 sigma1 = d(2)*eps1
 
!     FOR STIFFNESS MATRIX GENERATION, COMPUTE THE NEW MATERIAL
!     PROPERTIES
 
 41 IF (eps1 == eps2) GO TO 42
 e = (sigma2-sigma1)/(eps2-eps1)
 GO TO 44
 42 e = ecpt(20)
 
!     CALL MAT ROUTINE TO GET ELASTIC MODULI.  STORE IN D.P. LOCATIONS.
 
 44 matflg = 1
 CALL mat (iecpt(1))
 d(2) = e sub 0
 d(4) = gsub0
 
!     SET UP STIFFNESS MATRIX CONSTANTS IN DSCL AND DSCR
 
 g    = e*d(4)/d(2)
 d(1) = ecpt(5)
 d(3) = ecpt(6)
 dscl = d(1)*e/xl
 dscr = d(3)*g/xl
 
!     SET UP THE -N- MATRIX AND STORE AT D(1)
 
 d(1) = xn(1)*xn(1)
 d(2) = xn(1)*xn(2)
 d(3) = xn(1)*xn(3)
 d(4) = d(2)
 d(5) = xn(2)*xn(2)
 d(6) = xn(2)*xn(3)
 d(7) = d(3)
 d(8) = d(6)
 d(9) = xn(3)*xn(3)
 
!     ZERO OUT THE 6X6 WHICH WILL BE USED FOR STORAGE OF
!     KGG(NPVT,NONPVT), NONPVT = NPVT,J
 
 DO  i = 1,36
   ke(i) = 0.0D0
 END DO
 nonpvt = 2
 k2 = 1
 
!     IF PIVOT GRID POINT IS IN BASIC COORDINATES, GO TO 70
 
 IF (iecpt(ka) == 0) GO TO 70
 CALL gmmatd (ta(1),3,3,1, d(1),3,3,0, d(10))
 CALL gmmatd (d(10),3,3,0, ta(1),3,3,0, d(1))
 
!     AT THIS POINT D(1) CONTAINS THE MATRIX PRODUCT TAT*N*TA
!     AND D(10) CONTAINS THE MATRIX PRODUCT TAT*N.
 
 ASSIGN 100 TO iretrn
 GO TO  80
 70 ASSIGN 90 TO iretrn
 
!     FILL THE KE MATRIX
 
 80 ke( 1) = dscl*d(k2  )
 ke( 2) = dscl*d(k2+1)
 ke( 3) = dscl*d(k2+2)
 ke( 7) = dscl*d(k2+3)
 ke( 8) = dscl*d(k2+4)
 ke( 9) = dscl*d(k2+5)
 ke(13) = dscl*d(k2+6)
 ke(14) = dscl*d(k2+7)
 ke(15) = dscl*d(k2+8)
 ke(22) = dscr*d(k2  )
 ke(23) = dscr*d(k2+1)
 ke(24) = dscr*d(k2+2)
 ke(28) = dscr*d(k2+3)
 ke(29) = dscr*d(k2+4)
 ke(30) = dscr*d(k2+5)
 ke(34) = dscr*d(k2+6)
 ke(35) = dscr*d(k2+7)
 ke(36) = dscr*d(k2+8)
 CALL pla4b (ke,iecpt(nonpvt))
 
!     RETURN FROM FILL CODE W/ IRETRN =  90 IMPLIES G.P. A WAS IN BASIC
!       .     .    .     .      .     = 100 IMPLIES G.P. A WAS NOT BASIC
!       .     .    .     .      .     = 140 IMPLIES THE K(NPVT,NONPVT)
!                                       HAS BEEN COMPUTED AND INSERTED
!                                       AND HENCE WE ARE FINISHED.
 
 GO TO iretrn, (90,100,140)
 90 k1 = 1
 k2 = 10
 GO TO 110
 100 k1 = 10
 k2 = 1
 110 nonpvt = 3
 
!     IF NON-PIVOT GRID POINT IS IN BASIC COORDINATES, GO TO 120
 
 IF (iecpt(kb) == 0) GO TO 120
 
!     RECALL THAT D(K1) CONTAINS TAT*N.
 
 CALL gmmatd (d(k1),3,3,0, tb(1),3,3,0, d(k2))
 
!     AT THIS POINT D(K2) CONTAINS TAT*N*TB.
 
 GO TO 130
 120 k2 = k1
 130 ASSIGN 140 TO iretrn
 
!     SET CONSTANTS NEGATIVE TO PROPERLY COMPUTE K(NPVT,NONPVT)
 
 dscr = -dscr
 dscl = -dscl
 GO TO 80
 
!     A TRANSFER TO STATEMENT NO. 140 IMPLIES KGGNL CALCULATIONS HAVE
!     BEEN COMPLETED.  UPDATE ECPT ARRAY.
 
 140 IF (ind == 0) GO TO 150
 itemp    = iecpt(2)
 iecpt(2) = iecpt(3)
 iecpt(3) = itemp
 150 ecpt(18) = ecpt(19)
 ecpt(19) = eps1
 ecpt(20) = e
 RETURN
END SUBROUTINE pkrod
