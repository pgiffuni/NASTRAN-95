SUBROUTINE qdmm2 (temps,pg)
     
!     THERMAL LOAD GENERATION FOR THE QDMEM2 ELEMENT.
 
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
 
 
 REAL, INTENT(IN)                         :: temps(1)
 REAL, INTENT(OUT)                        :: pg(1)
 LOGICAL :: planar
 INTEGER :: nest(7),map(4,3)
 REAL :: rmat(3,5),et(9),k5sum(9,5),isinth,kmat(27),  &
     itemp9(9),alpha(3),pmat(9),jtemp9(9),icosth, gsube(9), psum(3,5),it
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /matin / matid, inflag, eltemp, stress, sinth, costh
 COMMON /matout/ g11,g12,g13,g22,g23,g33, rho, alps(3), tsub0
 COMMON /system/ ksystm(65)
 COMMON /trimex/ est(26)
 COMMON /condas/ pi,twopi,radeg,degra,s4pisq
 EQUIVALENCE     (ksystm(2),ioutpt),(nest(1),est(1))
 DATA    map   / 1, 2, 3, 4, 2, 3, 4, 1,  &
     5, 5, 5, 5 /
 
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
 
!     FORM  ALPHA = ALPS *(T-T )  3X1 VECTOR USED IN SUB-TRIANGLE CALCS
!                       E     0
 
 tbar     = temps(1) - tsub0
 alpha(1) = alps(1)*tbar
 alpha(2) = alps(2)*tbar
 alpha(3) = alps(3)*tbar
 
!     NOTE THE ABOVE MAY BE MOVED TO BELOW AND COMPUTED USING THE
!     GRID TEMPS OF SUB-TRIANGLE.  (I.E. TOTAL AVERAGE FOR CENTER POINT
!     ONLY.)  AVERAGE OF WHOLE ELEMENT IS USED EXCLUSIVELY NOW.
 
!     BASIC WHOLE-ELEMENT CALCULATIONS
 
 CALL q2bcs (est,planar,rmat,et,ierror)
 IF (ierror > 0) THEN
   GO TO   140
 END IF
 
!     ZERO SUMMATION ARRAYS
 
 10 DO  i = 1,5
   DO  j = 1,9
     k5sum(j,i) = 0.0
   END DO
   psum(1,i) = 0.0
   psum(2,i) = 0.0
   psum(3,i) = 0.0
 END DO
 
!     SUB-TRIANGLE COMPUTATIONS AND SUMMATIONS.
 
 DO  i = 1,4
   ia = map(i,1)
   ib = map(i,2)
   ic = map(i,3)
   it = est(8)
   CALL q2trms (rmat(1,ia),rmat(1,ib),rmat(1,ic),alpha(1),isinth,  &
       icosth,gsube,it,ierror,2,kmat,pmat,dummy,dummy)
   IF (ierror > 0) THEN
     GO TO   140
   END IF
   
!     SUM IN KCA,KCB,KCC
   
   40 DO  k = 1,9
     k5sum(k,ia) = k5sum(k,ia) + kmat(k   )
     k5sum(k,ib) = k5sum(k,ib) + kmat(k+ 9)
     k5sum(k,ic) = k5sum(k,ic) + kmat(k+18)
   END DO
   
!     SUM IN PA,PB,PC
   
   DO  k = 1,3
     psum(k,ia) = psum(k,ia) + pmat(k  )
     psum(k,ib) = psum(k,ib) + pmat(k+3)
     psum(k,ic) = psum(k,ic) + pmat(k+6)
   END DO
   
 END DO
 
!     IF -PLANAR- MODIFY THE K5SUM MATRICES.
 
 IF (.NOT.planar) GO TO 90
 DO  i = 1,5
   k5sum(7,i) = 0.0
   k5sum(8,i) = 0.0
   k5sum(9,i) =-0.25
 END DO
 k5sum(9,5) = 1.0
 
!     INVERT K   AND NEGATE THE RESULT.
!             55
 
 90 CONTINUE
 
!     NO NEED TO COMPUTE DETERMINANT SINCE IT IS NOT USED SUBSEQUENTLY.
 
 ising = -1
 CALL invers (3,k5sum(1,5),3,dummy,0,determ,ising,itemp9)
 IF (ising == 2) GO TO 140
 
 DO  i = 1,9
   k5sum(i,5) = -k5sum(i,5)
 END DO
 
!     4 (3X1) LOAD VECTORS ARE COMPUTED AND ADDED INTO THE P-VECTOR IN
!     CORE
 
!       G        T   T                  -1      T
!     (P ) = (T ) (E) ((PSUM ) + ((-K  ) (K  )) (PSUM ))
!       I      I            I        55    5I        5
 
 DO  i = 1,4
   CALL gmmats (k5sum(1,5),3,3,0,k5sum(1,i),3,3,0,itemp9)
   CALL gmmats (itemp9,3,3,1,psum(1,5),3,1,0,jtemp9)
   DO  j = 1,3
     psum(j,i) = psum(j,i) + jtemp9(j)
   END DO
   CALL gmmats (et,3,3,1,psum(1,i),3,1,0,jtemp9)
   jtemp9(4) = 0.0
   jtemp9(5) = 0.0
   jtemp9(6) = 0.0
   k = 4*i + 6
   IF (nest(k) /= 0) CALL basglb (jtemp9,jtemp9,nest(k+1),nest(k))
   
!     ADD LOAD TO CORE FOR THIS GRID
!                                   I
   l = nest(i+1)
   DO  j = 1,3
     pg(l) = pg(l) + jtemp9(j)
     l = l + 1
   END DO
   
 END DO
 RETURN
 
!     ERROR CONDITIONS
 
 140 WRITE  (ioutpt,150) uwm,nest(1)
 150 FORMAT (a25,' 3100, ELEMENT THERMAL LOAD COMPUTATION FOR QDMEM2 ',  &
     'ELEMENT ID =',i9, /5X,'FINDS ILLEGAL GEOMETRY THUS NO ',  &
     'LOADS OUTPUT FOR ELEMENT-ID NOTED.')
 RETURN
END SUBROUTINE qdmm2
