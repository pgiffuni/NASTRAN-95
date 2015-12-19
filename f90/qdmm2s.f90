SUBROUTINE qdmm2s
     
!     THIS ROUTINE CALCULATES THE STIFFNESS, MASS AND DAMPING MATRICES
!     FOR THE QDMM2 ELEMENT.
 
!     SINGLE PRECISION VERSION
 
!     THIS ROUTINE USES SUBROUTINE E MAS TQ TO CALCULATE THE LUMPED
!     MASS USING THE SAME METHOD AS WITH THE QDMEM ELEMENT.
 
!     THIS ROUTINE MAY NOT BE CALLED IN A HEAT PROBLEM.
 
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
 
 LOGICAL :: planar,nogo,iheat
 INTEGER :: dict(11),elid,estid,ipart(4),nest(7),map(4,3)
 REAL :: rmat(3,5),et(9),k1sum(9,16),kij(1),isinth,  &
     kmat(63),k5sum(9,5),icosth,gsube(9),it,g(36),  &
     itemp9(9),k5mod(9,5),tmat(36),jtemp9(9),idetrm, ktemp9(9),kout(144)
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm
 COMMON /emgest/ est(26)
 COMMON /emgprm/ dumm(15),ismd(3),iprec,nogo,heat,icmbar
 COMMON /emgdic/ dum(2),ngrids,elid,estid
 COMMON /matin / matid,inflag,eltemp,stress,sinth,costh
 COMMON /matout/ g11,g12,g13,g22,g23,g33,rho,alps(3),tsub0,GE
 COMMON /system/ ksystm(65)
 COMMON /condas/ pi,twopi,radeg,degra,s4pisq
 EQUIVALENCE     (ksystm(2),ioutpt),(nest(1),est(1)),  &
     (dict(5),dict5),(k1sum(1,1),kij(1)), (ksystm(56),iheat)
 DATA    map   / 1, 2, 3, 4, 2, 3, 4, 1,  &
     5, 5, 5, 5  /
 
 
!     THIS ELEMENT NOT USED IN A HEAT PROBLEM
 
 IF (iheat) GO TO 320
 
!     CREATE AN ARRAY POINTING TO THE GRID POINTS ACCORDING TO
!     INCREASING SIL VALUE
 
 DO  i = 1,4
   ipart(i) = nest(i+1)
 END DO
 i = -4
 4 j = 0
 DO  k = 1,4
   IF (ipart(k) < j) CYCLE
   j = ipart(k)
   l = k
 END DO
 ipart(l) = i
 i = i + 1
 IF (i < 0) GO TO 4
 DO  i = 1,4
   ipart(i) = -ipart(i)
 END DO
 
!     IF STIFFNESS MATRIX NEEDED
!     SET UP DICT ARRAY AND FOR STIFFNESS MATRIX
!     CALCULATIONS, OTHERWISE SKIP
 
 IF (ismd(1) == 0) GO TO 400
 
!      COMPUTE BASIC SIN AND COSINE OF ELEMENT MATERIAL ANGLE.
 
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
 
 CALL q2bcs (est,planar,rmat,et,ierror)
 IF (ierror > 0) THEN
   GO TO   270
 END IF
 
!     ZERO SUMMATION ARRAYS
 
 10 DO  i = 1,9
   DO  j = 1,16
     k1sum(i,j) = 0.0
   END DO
   DO  j = 1,5
     k5sum(i,j) = 0.0
   END DO
 END DO
 
!     SUB-TRIANGLES ARE COMPUTED AND RESULTS SUMMED.
 
 DO  i = 1,4
   
!     CALL TRIANGLE CALCULATION ROUTINE TO GET (3X3) SUB-PARTITIONS
   
   ia = map(i,1)
   ib = map(i,2)
   ic = map(i,3)
   it = est(8)
   
   CALL q2trms (rmat(1,ia),rmat(1,ib),rmat(1,ic),dummy,isinth,icosth,  &
       gsube,it,ierror,1,kmat,dummy,dummy,dummy)
   IF (ierror > 0) THEN
     GO TO   270
   END IF
   
!     SUM IN KCA,KCB,KCC 3-(3X3)-S STORED FIRST IN KMAT
   
!     ALSO SUM IN KAA,KAB,KBA,KBB = LAST 4-(3X3)-S STORED IN KMAT.
!     THESE GO INTO 4 OF THE 16 POSSIBLE (3X3) SUM MATRICES = ,
   
!     K11,K12,K13,K14,K21,K22,K23,K24,K31,K32,K33,K34,K41,K42,K43,K44
   
!     J1,J2,J3,J4 WILL EACH POINT TO 1 OF THE 16 (3X3)-S.
   
   50 j1 = 5*ia - 4
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
   END DO
   
 END DO
 
!     FORMATION OF THE FOUR (3X3) G MATRICES.
!                     -1
!     (G ) = -(K5SUM  ) (K  )   NOTE.  IF -PLANAR- THEN MODIFIED
!       I           55    5I           K5SUM MATRICES ARE USED.
 
 IF (planar) GO TO 90
 DO  i = 1,5
   DO  j = 1,9
     k5mod(j,i) = k5sum(j,i)
   END DO
 END DO
 GO TO 110
 
 90 DO  i = 1,5
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
 
 110 CONTINUE
 
!     NO NEED TO COMPUTE DETERMINANT SINCE IT IS NOT USED SUBSEQUENTLY.
 
 ising = -1
 CALL invers (3,k5mod(1,5),3,dummy,0,idetrm,ising,itemp9)
 IF (ising == 2) GO TO 290
 
 DO  i = 1,9
   k5mod(i,5) = -k5mod(i,5)
 END DO
 
!     FORM G MATRICES
 
 DO  i = 1,4
   CALL gmmats (k5mod(1,5),3,3,0, k5mod(1,i),3,3,0, g(9*i-8))
 END DO
 
!     FORMATION OF THE 4 TRANSFORMATION MATRICES EACH (3X3)
 
 DO  i = 1,4
   iest = 4*i + 6
   IF (nest(iest) == 0) THEN
     GO TO   150
   END IF
   
!     GET TRANSFORMATION MATRIX
   
   140 CALL transs (nest(iest),itemp9)
   CALL gmmats (et,3,3,0, itemp9,3,3,0, tmat(9*i-8))
   CYCLE
   
   150 k = 9*i - 9
   DO  j = 1,9
     k = k + 1
     tmat(k) = et(j)
   END DO
   
 END DO
 
!     FORM STIFFNESS MATRIX BY ROW-PARTIONS.
 
 DO  i = 1,4
!                          T
!     IF -PLANAR- FORM (G ) (K  ) FOR USE IN COLUMN-PARTITIONS LOOP.
!                        I    55
   
   IF (.NOT.planar) GO TO 190
   CALL gmmats (g(9*i-8),3,3,1, k5sum(1,5),3,3,0, itemp9)
   
!     COLUMN-PARTITIONS-LOOP
   
   190 DO  j = 1,4
!                                   T
!     FORM (K  ) = (K5SUM  ) + (K  ) (G )
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
     
     IF (.NOT.planar) GO TO 220
     CALL gmmats (itemp9,3,3,0, g(9*j-8),3,3,0, jtemp9)
     CALL gmmats (g(9*i-8),3,3,1, k5sum(1,j),3,3,0, ktemp9)
     DO  k = 1,9
       k1sum(k,lpart) = k1sum(k,lpart) + ktemp9(k) + jtemp9(k)
     END DO
     
!     TRANSFORM THIS RESULTANT K   (3X3) STORED AT K1SUM(1,LPART)
!                               IJ
!     TO GLOBAL.
     
     220 CALL gmmats (tmat(9*i-8),3,3,1, k1sum(1,lpart),3,3,0, jtemp9)
     CALL gmmats (jtemp9,3,3,0, tmat(9*j-8),3,3,0, k1sum(1,lpart))
   END DO
 END DO
 
!     FOR THE MATRIX ASSEMBLER -EMG- THE 16 (3X3) PARTITIONS IN K1SUM
!     ARE REARRANGED TO STORE THEM BY ROWS TO A TOTAL OF
!     12X12 RATHER THAN 3X3.  BUT FIRST DICT MUST BE
!     SET UP.  THE SILS MUST BE SORTED SO THAT THE 12X12 WILL
!     BE BY INCREASING SIL VALUE
 
 dict(1) = estid
 dict(2) = 1
 dict(3) = 12
 dict(4) = 7
 dict5   = GE
 ip      = iprec
 
!     REORDER K1SUM INTO KOUT AS DESCRIBED ABOVE
 
!         ****          ****
!         * K   K   K   K  *
!         *  AA  AB  AC  AD*
!     K = * K   K   K   K  *
!         *  BA  BB  BC  BD*
!         * K   K   K   K  *
!         *  CA  CB  CC  CD*
!         * K   K   K   K  *
!         *  DA  DB  DC  DD*
!         ****          ****
 
!     WHERE SUBSCRIPTS ARE ARRANGED BY INCREASING SIL VALUE
 
 DO  i = 1,4
   ii = ipart(i)
   DO  j = 1,4
     jtt = ipart(j)
     jt = (i-1)*4  + j
     DO  k = 1,9
       modk = MOD(k,3)
       IF (modk == 0) modk = 3
       l = (ii-1)*36 + ((k-1)/3)*12 + (jtt-1)*3 + modk
       kout(l) = k1sum(k,jt)
     END DO
   END DO
 END DO
 
 CALL emgout (kout,kout,144,1,dict,1,ip)
 
!     CALCULATE THE MASS MATRIX HERE.  SUBROUTINE
!     E MAS TQ IS USED TO GENERATE A LUMPED
!     MASS MATRIX EXACTLY LIKE A QDMEM ELEMENT
 
 400 IF (ismd(2) == 0) RETURN
 
 CALL emastq (1,k1sum)
 
 dict(1) = estid
 dict(2) = 2
 dict(3) = 12
 dict(4) = 7
 dict(5) = 0
 
!     REARRANGE KIJ BY INCREASING SIL VALUE
 
 DO  i = 1,4
   ii = 1 + (ipart(i)-1)*3
   ij = (i-1)*3 + 1
   kout(ij  ) = kij(ii  )
   kout(ij+1) = kij(ii+1)
   kout(ij+2) = kij(ii+2)
 END DO
 
 CALL emgout (kout,kout,12,1,dict,2,ip)
 RETURN
 
!     ELEMENT ERRORS DETECTED.
 
 270 WRITE  (ioutpt,280) ufm,nest(1)
 280 FORMAT (a23,' 3098,  QDMEM2 ELEMENT STIFFNESS ROUTINE DETECTS ',  &
     'ILLEGAL GEOMETRY FOR ELEMENT ID =',i10)
 GO TO 310
 290 WRITE  (ioutpt,300) ufm,nest(1)
 300 FORMAT (a23,' 3099.  ELEMENT STIFFNESS COMPUTATION FOR QDMEM2 ',  &
     'ELEMENT ID =',i10, /5X,'IS IMPOSSIBLE DUE TO SINGULARITY',  &
     ' IN CONSTRAINT EQUATION.')
 310 nogo = .true.
 RETURN
 
 320 WRITE  (ioutpt,330) uwm,nest(1)
 330 FORMAT (a25,' 3115, QDMM2 FINDS ELEMENT NUMBER',i10,  &
     ' PRESENT IN A HEAT FORMULATION AND IS IGNORING SAME.')
 
 RETURN
END SUBROUTINE qdmm2s
