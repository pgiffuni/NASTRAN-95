SUBROUTINE shgmgs (*,elid,tem,mid,ts,noalfa,g,rho,gsube,tsub0,  &
        egnor,alpha)
     
!     MATERIAL PROPERTY G-MATRICES GENERATOR FOR SHELL ELEMENTS
 
!     SINGLE PRECISION VERSION
 
!     INPUT :
!           ELID   - ELEMENT ID
!           TEM    - 3X3 TRANSFORMATION BETWEEN ELEMENT AND MATERIAL
!                    COORDINATE SYSTEMS
!           MID    - ARRAY OF LENGTH 4, CONTAINS MATERIAL ID'S
!           TS     - EQUIVALENT SHEAR THICKNESS
!           NOALFA - LOGICAL TO REQUEST OR NOT REQUEST THERMAL EXPANSION
!                    COEFFICIENTS
 
!     OUTPUT:
!           G      - ARRAY OF LENGTH 36 (FOUR 3X3), CONATINS MATERIAL
!                    PROPERTIES IN ELEMENT COORD. SYSTEM
!           RHO    - MASS DENSITY FROM MEMBRANE MATERIAL
!           GSUBE  - DAMPING COEFFICIENT FROM MEMBRANE OR BENDING
!                    MATERIALS
!           TSUB0  - REFERENCE TEMPERATURE
!           EGNOR  - ARRAY OF PSEUDO E'S AND G'S FOR SHEAR FACTOR
!                    CALCULATIONS IN BENDING
!           ALPHA  - ARRAY OF THERMAL EXPANSION COEFFICIENTS
 
!     NOTES:
!           1- THIS ROUTINE BUILDS THE MATERIAL PROPERTY MATRIX USING
!              THE OUTPUT OF SUBROUTINE 'MAT' (/MATOUT/).
!              /MATOUT/ IS IN MAT2 FORMAT IF MAT1 AND/OR MAT2 ARE USED
!              /MATOUT/ IS IN MAT8 FORMAT IF MAT8 CARD IS REQUESTED.
!           2- ISOTROPIC, ORTHOTROPIC, AND ANISOTROPIC PROPERTY TYPES
!              ARE SUPPORTED.
!           3- PROPERTIES FOR MEMBRANE, BENDING, SHEAR FLEXIBILITY, AND
!              MEMBRANE/BENDING COUPLING ARE PROCESSED.
!           4- NON-STANDARD RETURN IS TAKEN WHEN THE MATERIAL FOR SHEAR
!              FLEXIBILITY IS NOT PROPERLY DEFINED.
!           5- SOME OF THE CONTENTS OF /MATIN/ MUST BE DEFINED PRIOR TO
!              A CALL TO THIS ROUTINE.
!           6- CONTENTS OF /TERMS/, MID, AND TS MAY BE CHANGED IN THIS
!              ROUTINE.
 
 
 
 
 , INTENT(IN OUT)                         :: *
 INTEGER, INTENT(IN OUT)                  :: elid
 REAL, INTENT(IN)                         :: tem(9)
 INTEGER, INTENT(IN OUT)                  :: mid(4)
 REAL, INTENT(OUT)                        :: ts
 LOGICAL, INTENT(IN OUT)                  :: noalfa
 REAL, INTENT(OUT)                        :: g(36)
 REAL, INTENT(OUT)                        :: rho
 REAL, INTENT(OUT)                        :: gsube
 REAL, INTENT(OUT)                        :: tsub0
 REAL, INTENT(OUT)                        :: egnor(4)
 REAL, INTENT(OUT)                        :: alpha(6)
 LOGICAL :: membrn,bendng,shrflx,mbcoup,norpth
 INTEGER :: INDEX(3,3), NAME(2)
 REAL :: nu12,nu21,matset
 REAL :: u(9),GT(9), dn12,dn21,ps1, ps2, const, talpha(6),detu,bdum
 COMMON /terms /  membrn,bendng,shrflx,mbcoup,norpth
 COMMON /matin /  matid,inflag,eltemp,dummy,sinmat,cosmat
 COMMON /matout/  g11,g12,g13,g22,g23,g33,rhox,alph1,alph2,alph12,  &
     tref,GE,st,sc,ss,e,dum(8),matset
!                      MAT8 FORMAT...
 EQUIVALENCE      (e1 ,g11),(nu12,g12),(e2,g13),(g2z,g23),  &
     (g1z,g33),(g12x,g22)
!    2,                (GE  ,E  ),(T0,ALPH12)
!             EQUIV    (MATOUT(1),G11))
 DATA    NAME  /  4HSHGM,4HGS    /
 
 
!     INITIALIZE
 
!     SET INFLAG = 12 SO THAT SUBROUTINE MAT WILL SEARCH FOR:
!     ISOTROPIC   MATERIAL PROPERTIES AMONG THE MAT1 CARDS,
!     ORTHOTROPIC MATERIAL PROPERTIES AMONG THE MAT8 CARDS, AND
!     ANISOTROPIC MATERIAL PROPERTIES AMONG THE MAT2 CARDS.
 
 DO  ig = 1,36
   g(ig) = 0.0
 END DO
 DO  ig = 1,4
   egnor(ig) = 0.0
 END DO
 IF (noalfa) GO TO 40
 DO  ig = 1,6
   alpha(ig) = 0.0
 END DO
 
 40 rho   = 0.0
 gsube = 0.0
 tsub0 = 0.0
 inflag= 12
 igobk = 0
 it0   = 0
 
!     BEGIN THE LOOP TO FETCH PROPERTIES FOR EACH MATERIAL ID. FOR SHEAR
!     FLEXIBILITY MATERIAL, DEFAULT TO THE BENDING MATERIAL IF BENDING
!     IS PRESENT.
!     IF SHEAR MATERIAL IS PRESENT, BUT YIELDS ZEROES, GO BACK AND RESET
!     IT TO BENDING MATERIAL.
 
 m = 0
 100 lpoint = m*9
 m = m + 1
 IF (m > 4) GO TO 600
 IF (m == 4 .AND. igobk == 1) GO TO 610
 matid = mid(m)
 IF (matid == 0 .AND. m /= 3) GO TO 100
 IF (matid == 0 .AND. m == 3. AND. .NOT.bendng) GO TO 100
 IF (matid == 0 .AND. m == 3. AND. bendng) matid = mid(2)
 
 IF (m-1 < 0) THEN
   GO TO   130
 ELSE IF (m-1 == 0) THEN
   GO TO   120
 END IF
 110 IF (matid == mid(m-1) .AND. igobk == 0) GO TO 130
 120 CALL mat (elid)
 tmtset = matset
 IF (matset == 8.0) tmtset = 3.0
 mtype = IFIX(tmtset+0.05) - 2
 
!     SET THE MISC ITEMS
 
 130 IF (membrn .AND. m == 1) rho = rhox
 IF (membrn .AND. m /= 1 .OR. .NOT.membrn .AND. m /= 2) GO TO 140
 gsube = GE
 IF (mtype > 0) gsube = e
 140 IF (it0 > 0) GO TO 150
 it0   = 1
 tsub0 = tref
 IF (mtype > 0) tsub0 = alph12
 
!     BRANCH ON MATERIAL TYPE
 
 150 IF (mtype < 0) THEN
   GO TO   200
 ELSE IF (mtype == 0) THEN
   GO TO   210
 ELSE
   GO TO   250
 END IF
!               MAT1,MAT2,MAT8
 
 
!     ISOTROPIC  MATERIALS (MAT1)
!     ---------------------------
 
! 200 IF (M .NE. 3) GO TO 205
 200 IF (m /= 3) GO TO 220
 
!     G(LPOINT+1) = MATOUT(3)   <== G13, SHOULD BE MATOUT(6) <== G33
!     G(LPOINT+4) = G(LPOINT+1)
!     IF (G(LPOINT+1).EQ.0.0 .AND. SHRFLX) GO TO 300
 
 g(lpoint+1) = g33
 g(lpoint+4) = g33
 IF (g33 == 0.0 .AND. shrflx) GO TO 300
 GO TO 400
 
!     ACCORDING TO Q4GMGS, SHOULD TO TO 220 NEXT
 
! 205 G(LPOINT+1) = G22
!     G(LPOINT+2) = G12*G22
!     G(LPOINT+4) = G12*G22
!     G(LPOINT+5) = G22
!     G(LPOINT+9) = G13         <== G13,  SHOULD IT BE G33 ??
!     GO TO 400
 
!     ANISOTROPIC  MATERIALS (MAT2)
!     -----------------------------
 
 210 IF (m == 3) GO TO 230
 220 g(lpoint+1) = g11
 g(lpoint+2) = g12
 g(lpoint+3) = g13
 g(lpoint+4) = g12
 g(lpoint+5) = g22
 g(lpoint+6) = g23
 g(lpoint+7) = g13
 g(lpoint+8) = g23
 g(lpoint+9) = g33
 GO TO 400
 
 230 IF (shrflx) GO TO 240
 IF (g11 == 0.0 .OR. g22 == 0.0) GO TO 400
 dn21  = g12/g11
 dn12  = g12/g22
 const = dn21*dn12
 IF (const < 0.0) GO TO 400
 ps1 = g11*(1.0-const)
 ps2 = g22*(1.0-const)
 IF (const > 0.0) const = SQRT(const)
 const = 2.0*(1.0+const)
 g(lpoint+1) = ps1/const
 g(lpoint+4) = ps2/const
 GO TO 400
 
 240 g(lpoint+1) = g11
 g(lpoint+2) = g12
 g(lpoint+3) = g12
 g(lpoint+4) = g22
 IF (g33 /= 0.0) GO TO 300
 GO TO 400
 
!     ORTHOTROPIC MATERIALS (MAT8)
!     ----------------------------
 
 250 IF (m  ==   3) GO TO 260
 IF (e1 == 0.0) GO TO 400
 nu21 = nu12*e2/e1
 const= 1.0 - nu21*nu12
 IF (const <= 0.0) GO TO 400
 g(lpoint+1) = e1/const
 g(lpoint+2) = nu12*e2/const
 g(lpoint+4) = g(lpoint+2)
 g(lpoint+5) = e2/const
 g(lpoint+9) = g12x
 GO TO 400
 
 260 IF (shrflx) GO TO 270
 IF (e1 == 0.0) GO TO 400
 nu21  = nu12*e2/e1
 const = nu21*nu12
 IF (const <= 0.0) GO TO 400
 const = SQRT(const)
 const = 2.0*(1.0+const)
 g(lpoint+1) = e1/const
 g(lpoint+4) = e2/const
 GO TO 400
 
! 270 G(LPOINT+1) = MATOUT(5)         <== COSMIC (5) & (6) INTERCHANGED
!     G(LPOINT+4) = MATOUT(6)
 270 g(lpoint+1) = g1z
 g(lpoint+4) = g2z
 IF (g1z == 0.0 .AND. g2z == 0.0) GO TO 300
 GO TO 400
 
!     BAD SHEAR MATERIAL
 
 300 IF (.NOT.shrflx .AND. bendng) GO TO 400
 RETURN 1
 
!     TRANSFORM NON-ISOTROPIC MATERIALS
 
 400 IF (mtype < 0) GO TO 430
 IF (m     == 3) GO TO 410
 u(1) = tem(1)*tem(1)
 u(2) = tem(4)*tem(4)
 u(3) = tem(1)*tem(4)
 u(4) = tem(2)*tem(2)
 u(5) = tem(5)*tem(5)
 u(6) = tem(2)*tem(5)
 u(7) = tem(1)*tem(2)*2.0
 u(8) = tem(4)*tem(5)*2.0
 u(9) = tem(1)*tem(5) + tem(2)*tem(4)
 l    = 3
 GO TO 420
 
 410 u(1) = tem(5)*tem(9) + tem(6)*tem(8)
 u(2) = tem(2)*tem(9) + tem(8)*tem(3)
 u(3) = tem(4)*tem(9) + tem(7)*tem(6)
 u(4) = tem(1)*tem(9) + tem(3)*tem(7)
 l    = 2
 
 420 CALL gmmats ( u(1),l,l,1, g(lpoint+1),l,l,0, GT(1))
 CALL gmmats (GT(1),l,l,0, u(1),l,l,0,  g(lpoint+1))
 
!     GET THE THERMAL EXPANSION COEFFICIENTS, IF NEEDED
 
 430 IF (noalfa .OR. m > 2) GO TO 100
 morb = (m-1)*3
 IF (mtype < 0) THEN
   GO TO   500
 ELSE IF (mtype == 0) THEN
   GO TO   510
 ELSE
   GO TO   520
 END IF
!                MAT1,MAT2,MAT8
 
!     MAT1
 
 500 alpha(morb+1) = alph1
 alpha(morb+2) = alph1
 alpha(morb+3) = 0.0
 GO TO 100
 
!     MAT2
 
 510 alpha(morb+1) = alph1
 alpha(morb+2) = alph2
 alpha(morb+3) = alph12
 GO TO 530
 
!     MAT8
 
 520 alpha(morb+1) = alph1
 alpha(morb+2) = alph2
 alpha(morb+3) = 0.0
 
!     TRANSFORM THERMAL EXPANSION COEFFICIENTS AND STORE THEM IN ALPHA.
!     THE ALPHAS NEED TO BE PREMULTIPLIED BY [U] INVERSE.
 
 530 DO  ig = 1,3
   talpha(ig+morb) = alpha(ig+morb)
 END DO
 morb = morb + 1
 CALL invers (3,u,3,bdum,0,detu,isngu,INDEX)
 CALL gmmats (u,3,3,0, talpha(morb),3,1,0, alpha(morb))
 GO TO 100
 
 
!     LOOP IS DONE, CHECK FOR ALL ZEROES FOR SHEAR MATERIAL
 
 600 IF (g(19) /= 0.0 .OR. g(20) /= 0.0 .OR. g(21) /= 0.0 .OR.  &
     g(22) /= 0.0) GO TO 610
 igobk  = 1
 m      = 2
 mid(3) = 0
 shrflx = .false.
 ts     = 0.833333333
!              0.833333333 = 5.0/6.0
 GO TO 100
 
!     SAVE PSEUDO E'S AND G'S FOR SHEAR FACTOR CALCULATIONS
 
 610 IF (.NOT.bendng) GO TO 620
 egnor(1) = g(10)
 egnor(2) = g(14)
 egnor(3) = g(19)
 egnor(4) = g(22)
 
 620 RETURN
END SUBROUTINE shgmgs
