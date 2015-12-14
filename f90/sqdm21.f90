SUBROUTINE sqdm21
     
!     PHASE-I STRESS-DATA-RECOVERY ROUTINE FOR THE -QDMEM2- ELEMENT.
 
!     THIS ROUTINE WILL PREPARE FOR USE BY -SQDM22-, THE PHASE-II
!     ROUTINE, A TABLE CONTAINING THE FOLLOWING.
 
!     TABLE WORDS        DISCRIPTION
!     ------------------------------------------------------
!       1 THRU   1       ELEMENT-ID
!       2 THRU   5       4 SILS
!       6 THRU   6       ELEMENT-THICKNESS
!       7 THRU   7       REFERENCE TEMP -TSUB0-
!       8 THRU 151       16 (3X3) KIJ-G MATRICES
!     152 THRU 187       4 (3X3) STRESS MATRICES
!     188 THRU 199       4 (3X1) TEMP VECTORS
!     200 THRU 202       ST (3X1) STRESS-TEMPERATURE VECTOR
!     203 THRU 206       4 SIDE LENGTHS
 
!     ELEMENT EST ENTRY CONTENTS
!     + + + + + + + + + + + + + + + + + + + + + + + + + +
!     +   1 = ID                                        +
!     +   2 = SIL-PT-A            (ELEMENT CONNECTS     +
!     +   3 = SIL-PT-B             GRID POINTS A,B,     +
!     +   4 = SIL-PT-C             C,D IN THAT ORDER)   +
!     +   5 = SIL-PT-D                                  +
!     +   6 = MATERIAL-ANGLE                            +
!     +   7 = MATERIAL-ID                               +
!     +   8 = THICKNESS OF ELEMENT                      +
!     +   9 = NON-STRUCTURAL-MASS                       +
!     +  10 = COORD-SYS-ID PT-A OR 0                    +
!     +  11 = XA                                        +
!     +  12 = YA                                        +
!     +  13 = ZA                                        +
!     +  14 = COORD-SYS-ID PT-B OR 0                    +
!     +  15 = XB                                        +
!     +  16 = YB                                        +
!     +  17 = ZB                                        +
!     +  18 = COORD-SYS-ID PT-C OR 0                    +
!     +  19 = XC                                        +
!     +  20 = YC                                        +
!     +  21 = ZC                                        +
!     +  22 = COORD-SYS-ID PT-D OR 0                    +
!     +  23 = XD                                        +
!     +  24 = YD                                        +
!     +  25 = ZD                                        +
!     +  26 = AVERAGE OF CONNECTED GRID TEMPERATURES    +
!     + + + + + + + + + + + + + + + + + + + + + + + + + +
 
 LOGICAL :: planar
 INTEGER :: nest(7), map(4,3)
 REAL :: k1sum, k5sum, isinth, icosth, kmat(63), smat(27),  &
     pmat(9), jtemp9, k5mod, ktemp9(9), zmat(9),  &
     itemp9(9), q(3,3,4), imat12, rmat(3,5), eti(36), dvec(3,4), kvec(3)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm, uwm
 
!     FOLLOWING COMMON BLOCK MUST BE DIMENSIONED AT LEAST 350 IN SDR2B
 
 COMMON /sdr2x5/ est(100), id, isils(4), elthik, reftmp,  &
     k1sum(9,16), sg(36),pt(3,4), st(3),rg(4)
 
!     WORKING STORAGE BLOCK (KEEP .LE. 300 WORDS)
 
 COMMON /sdr2x6/ k5sum(9,5), sisum(9,5), pisum(3,5), r(3,4,5),  &
     k5mod(9,5), g(36),  t(9), e(9), imat12(12), jtemp9(9), gsube(9)
 COMMON /matin / matid, inflag, eltemp, stress, sinth, costh
 COMMON /matout/ g11, g12, g13, g22, g23, g33, rho, alps(3), tsub0
 COMMON /system/ ksystm(65)
 COMMON /condas/ pi, twopi, radeg, degra, s4pisq
 EQUIVALENCE     (ksystm(2),ioutpt), (nest(1),est(1))
 DATA    map   / 1, 2, 3, 4, 2, 3, 4, 1,  &
     5, 5, 5, 5  /
 
!     COMPUTE BASIC SIN AND COSINE OF ELEMENT MATERIAL ANGLE.
 
 angl   = est(6)*degra
 isinth = SIN(angl)
 icosth = COS(angl)
 
!     COMPUTE GSUBE MATRIX
 
 inflag = 2
 matid  = nest(7)
 eltemp = est(26)
 sinth  = 0.0
 costh  = 1.0
 CALL mat (nest(1))
 gsube(1) = g11
 gsube(2) = g12
 gsube(3) = g13
 gsube(4) = g12
 gsube(5) = g22
 gsube(6) = g23
 gsube(7) = g13
 gsube(8) = g23
 gsube(9) = g33
 
!     BASIC WHOLE-ELEMENT CALCULATIONS
 
 CALL q2bcs (est,planar,rmat,e,ierror)
 IF (ierror > 0) THEN
   GO TO   400
 END IF
 
!     ZERO SUMMATION ARRAYS
 
 10 DO  i = 1,9
   DO  j = 1,16
     k1sum(i,j) = 0.0
   END DO
   DO  j = 1,5
     k5sum(i,j) = 0.0
     sisum(i,j) = 0.0
   END DO
 END DO
 
 DO  i = 1,5
   pisum(1,i) = 0.0
   pisum(2,i) = 0.0
   pisum(3,i) = 0.0
 END DO
 
!     SUB-TRIANGLES ARE COMPUTED AND RESULTS SUMMED.
 
 DO  i = 1,4
   
!     CALL TRIANGLE CALCULATION ROUTINE TO GET (3X3) SUB-PARTITIONS
   
   ia = map(i,1)
   ib = map(i,2)
   ic = map(i,3)
   
   CALL q2trms (rmat(1,ia),rmat(1,ib),rmat(1,ic),alps ,isinth,icosth,  &
       gsube,est(8),ierror,3,kmat,pmat,smat,zmat)
   IF (ierror > 0) THEN
     GO TO   400
   END IF
   
!     SUM IN KCA,KCB,KCC 3-(3X3)-S STORED FIRST IN KMAT
   
!     ALSO SUM IN KAA,KAB,KBA,KBB = LAST 4-(3X3)-S STORED IN KMAT.
!     THESE GO INTO 4 OF THE 16 POSSIBLE (3X3) SUM MATRICES = ,
   
!     K11,K12,K13,K14,K21,K22,K23,K24,K31,K32,K33,K34,K41,K42,K43,K44
   
!     J1,J2,J3,J4 WILL EACH POINT TO 1 OF THE 16 (3X3)-S.
   
   70 j1 = 5*ia - 4
   j2 = 4*ia - 4 + ib
   j3 = 4*ib - 4 + ia
   j4 = 5*ib - 4
   
   DO  k = 1,9
     k5sum(k,ia) = k5sum(k,ia) + kmat(k   )
     k5sum(k,ib) = k5sum(k,ib) + kmat(k+ 9)
     k5sum(k,ic) = k5sum(k,ic) + kmat(k+18)
     k1sum(k,j1) = k1sum(k,j1) + kmat(k+27)
     k1sum(k,j2) = k1sum(k,j2) + kmat(k+36)
     k1sum(k,j3) = k1sum(k,j3) + kmat(k+45)
     k1sum(k,j4) = k1sum(k,j4) + kmat(k+54)
     sisum(k,ia) = sisum(k,ia) + smat(k   )
     sisum(k,ib) = sisum(k,ib) + smat(k+ 9)
     sisum(k,ic) = sisum(k,ic) + smat(k+18)
   END DO
   
   DO  k = 1,3
     pisum(k,ia) = pisum(k,ia) + pmat(k  )
     pisum(k,ib) = pisum(k,ib) + pmat(k+3)
     pisum(k,ic) = pisum(k,ic) + pmat(k+6)
   END DO
   
 END DO
 
!     FORMATION OF THE FOUR (3X3) G MATRICES.
!                     -1
!     (G ) = -(K5SUM  ) (K  )   NOTE.  IF -PLANAR- THEN MODIFIED
!       I           55    5I           K5SUM MATRICES ARE USED.
 
 IF (planar) GO TO 120
 DO  i = 1,5
   DO  j = 1,9
     k5mod(j,i) = k5sum(j,i)
   END DO
 END DO
 GO TO 140
 
 120 DO  i = 1,5
   k5mod(1,i) = k5sum(1,i)
   k5mod(2,i) = k5sum(2,i)
   k5mod(3,i) = k5sum(3,i)
   k5mod(4,i) = k5sum(4,i)
   k5mod(5,i) = k5sum(5,i)
   k5mod(6,i) = k5sum(6,i)
   k5mod(7,i) = 0.0
   k5mod(8,i) = 0.0
   k5mod(9,i) =-0.25
 END DO
 k5mod(9,5) = 1.0
 
!     INVERT K5MOD   AND NEGATE RESULT.
!                 55
 
 140 CONTINUE
 
!     NO NEED TO COMPUTE DETERMINANT SINCE IT IS NOT USED SUBSEQUENTLY.
 
 ising = -1
 CALL invers (3,k5mod(1,5),3,dummy,0,determ,ising,itemp9)
 IF (ising == 2) GO TO 400
 
 DO  i = 1,9
   k5mod(i,5) = -k5mod(i,5)
 END DO
 
!     FORM G MATRICES
 
 DO  i = 1,4
   CALL gmmats (k5mod(1,5),3,3,0, k5mod(1,i),3,3,0, g(9*i-8))
 END DO
 
!     FORM STIFFNESS MATRIX BY ROW-PARTIONS.
 
 DO  i = 1,4
!                          T
!     IF -PLANAR- FORM (G ) (K  ) FOR USE IN COLUMN-PARTITIONS LOOP.
!                        I    55
   
   IF (.NOT.planar) GO TO 170
   CALL gmmats (g(9*i-8),3,3,1, k5sum(1,5),3,3,0, itemp9)
   
!     COLUMN-PARTITIONS-LOOP
   
   170 DO  j = 1,4
!                                   T
!     FORM (K  ) = (K1SUM  ) + (K  ) (G )
!            IJ          IJ      5I    J
     
     CALL gmmats (k5sum(1,i),3,3,1, g(9*j-8),3,3,0, jtemp9)
     lpart = 4*i - 4 + j
     DO  k = 1,9
       k1sum(k,lpart) = k1sum(k,lpart) + jtemp9(k)
     END DO
     
!     BALANCE OF TERMS IF -PLANAR-
     
!                T            T
!     ADD IN (G ) (K  ) + (G ) (K  )(G )
!              I    5J      I    55   J
     
     IF (.NOT.planar) CYCLE
     CALL gmmats (itemp9,3,3,0, g(9*j-8),3,3,0, jtemp9)
     CALL gmmats (g(9*i-8),3,3,1, k5sum(1,j),3,3,0, ktemp9)
     DO  k = 1,9
       k1sum(k,lpart) = k1sum(k,lpart) + ktemp9(k) + jtemp9(k)
     END DO
   END DO
 END DO
 
!     CALCULATION OF 4 (Q ) MATRICES, EACH 3X3.
!                        I
 
 DO  i = 1,4
   ia = map(i,1)
   ib = map(i,2)
   DO  j = 1,3
     dvec(j,i) = rmat(j,ib) - rmat(j,ia)
   END DO
   fmag  = SQRT(sadotb(dvec(1,i),dvec(1,i)))
   rg(i) = fmag
   IF (fmag > 0.0) THEN
     GO TO   240
   ELSE
     GO TO   400
   END IF
   240 DO  j = 1,3
     dvec(j,i) = dvec(j,i)/fmag
   END DO
 END DO
 
 DO  i = 1,4
   j = i - 1
   IF (j == 0) j = 4
   i1 = map(j,1)
   i2 = map(j,2)
   CALL saxb (dvec(1,i2),dvec(1,i1),kvec)
   
!     NORMALIZE, NEGATE, AND STORE AS DELTA-VEC IN (Q )
!                                                    I
   fmag = SQRT(sadotb(kvec,kvec))
   IF (fmag > 0.0) THEN
     GO TO   270
   ELSE
     GO TO   400
   END IF
   270 q(1,3,i) = -kvec(1)/fmag
   q(2,3,i) = -kvec(2)/fmag
   q(3,3,i) = -kvec(3)/fmag
   
!     STORE D VECTORS AS ALPHA- VECTORS IN (Q )
!                                            I
   q(1,1,i) = -dvec(1,i)
   q(2,1,i) = -dvec(2,i)
   q(3,1,i) = -dvec(3,i)
   
   q(1,2,i) = dvec(1,j)
   q(2,2,i) = dvec(2,j)
   q(3,2,i) = dvec(3,j)
   
!     INVERT 3X3
   
!     AGAIN NO NEED TO COMPUTE DETERMINANT SINCE IT IS NOT USED .
!     SUBSEQUENTLY.
   
   ising = -1
   CALL invers (3,q(1,1,i),3,dummy,0,determ,ising,jtemp9)
   IF (ising == 2) GO TO 400
 END DO
 
!     FORM FINAL OUTPUTS
 
 DO  i = 1,4
   ii = 9*i - 8
   
!     TRANSFORMATION ETI = (E)(T )
!                               I
   
   kk = 4*i
   IF (nest(kk+6) == 0) THEN
     GO TO   300
   END IF
   290 CALL transs (nest(kk+6),t)
   CALL gmmats (e,3,3,0, t,3,3,0, eti(ii))
   GO TO 320
   
   300 kk = ii
   DO  j = 1,9
     eti(kk) = e(j)
     kk = kk + 1
   END DO
   
!       G            E      E
!     (S ) = 0.25( (S ) + (S )(G ) )(E)(T )
!       I            I      5   I        I
   
   320 CALL gmmats (sisum(1,5),3,3,0, g(ii),3,3,0, jtemp9)
   DO  j = 1,9
     jtemp9(j) = 0.25*(jtemp9(j)+sisum(j,i))
   END DO
   CALL gmmats (jtemp9,3,3,0, eti(ii),3,3,0, sg(ii))
   
!       T     -         T -
!     (P ) = (P ) + (G ) (P )
!       I      I      I    5
   
   CALL gmmats (g(ii),3,3,1, pisum(1,5),3,1,0, pt(1,i))
   DO  j = 1,3
     pisum(j,i) = pt(j,i) + pisum(j,i)
   END DO
   CALL gmmats (q(1,1,i),3,3,1, pisum(1,i),3,1,0, pt(1,i))
 END DO
 
!     TRANSFORM STIFFNESS MATRIX TO GLOBAL
 
!        G           E
!     (K  ) = (Q )(K  )(E)(T )
!       IJ      I   IJ      J
 
 jpart = 0
 DO  i = 1,4
   DO  j = 1,4
     jpart = jpart + 1
     CALL gmmats (q(1,1,i),3,3,1, k1sum(1,jpart),3,3,0, jtemp9)
     CALL gmmats (jtemp9,3,3,0, eti(9*j-8),3,3,0, k1sum(1,jpart))
   END DO
 END DO
 
!     (S ) = (GSUBE)(ALPHAS)
!       T
 
 CALL gmmats (gsube,3,3,0, alps,3,1,0, st)
 
!     MISC. DATA FOR PHASE-II
 
 id       = nest(1)
 isils(1) = nest(2)
 isils(2) = nest(3)
 isils(3) = nest(4)
 isils(4) = nest(5)
 elthik   = est(8)
 reftmp   = tsub0
 RETURN
 
!     ERROR CONDITION
 
 400 WRITE  (ioutpt,410) uwm,nest(1)
 410 FORMAT (a25,' 3101, SINGULARITY OR BAD GEOMETRY FOR QDMEM2 ELEM.',  &
     ' ID =',i9, /5X,'STRESS OR FORCES WILL BE INCORRECT.')
 RETURN
END SUBROUTINE sqdm21
