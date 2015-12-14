SUBROUTINE t3bgbs (ng,nb,gmat,bmat,kmat)
     
!     WITH ENTRY T3BGBD (NG,NB,GMAD,BMAD,KMAD)
 
!     ROUTINE FOR EFFICIENT TRIPLE-MULTPLICATION OF [B] AND [G] MATRICES
!     TO EVALUATE THE CONTRIBUTION TO THE ELEMENT STIFFNESS MATRIX FROM
!     THE CURRENT INTEGRATION POINT
 
 
!     INPUT :
!           NG          - NUMBER OF ROWS AND COLUMNS OF GMAT
!           NB          - NUMBER OF COLUMNS OF BMAT
!           GMAT/GMAD   - [G], FORCE-STRAIN RELATIONSHIP
!           BMAT/BMAD   - [B], STRAIN-DISPLACEMENT RELATIONSHIP
!     OUTPUT:
!           KMAT/KMAD   - CONTRIBUTION TO THE ELEMENT STIFFNESS MATRIX
!                         FROM THE CURRENT INTEGRATION POINT
 
!     ALGORITHM:
!           MATRICES ARE MULTIPLIED IN FULL WHEN MEMBRANE-BENDING
!           COUPLING IN PRESENT, OTHERWISE PARTIAL MULTIPLICATION
!           IS PERFORMED.
!           IN EACH TRIPLE MULTIPLY, THE RESULT IS ADDED TO KMAT.
 
 
 
 INTEGER, INTENT(IN OUT)                  :: ng
 INTEGER, INTENT(IN)                      :: nb
 REAL, INTENT(IN)                         :: gmat(9,1)
 REAL, INTENT(IN OUT)                     :: bmat(1)
 REAL, INTENT(IN OUT)                     :: kmat(1)
 LOGICAL :: membrn,bendng,shrflx,mbcoup,norpth
 REAL :: g1(3,3),gbmat(162)
 DOUBLE PRECISION :: gmad(9,1),bmad(1),kmad(1),g2(3,3),gbmad(162)
 COMMON /terms /  membrn,bendng,shrflx,mbcoup,norpth
 EQUIVALENCE      (g1(1,1),g2(1,1)),(gbmat(1),gbmad(1))
 
 
!     SINGLE PRECISION
 
 nd3 = nb*3
 nd6 = nb*6
 
!     IF [G] IS FULLY POPULATED, PERFORM STRAIGHT MULTIPLICATION AND
!     RETURN.
 
 IF (.NOT.mbcoup) GO TO 10
 CALL gmmats (gmat,ng,ng,0,  bmat,ng,nb,0,  gbmat)
 CALL gmmats (bmat,ng,nb,-1, gbmat,ng,nb,0, kmat )
 GO TO 60
 
!     MULTIPLY MEMBRANE TERMS WHEN PRESENT
 
 10 IF (.NOT.membrn) GO TO 30
 DO  i = 1,3
   DO  j = 1,3
     g1(i,j) = gmat(i,j)
   END DO
 END DO
 CALL gmmats (g1,3,3,0, bmat(1),3,nb,0, gbmat)
 CALL gmmats (bmat(1),3,nb,-1, gbmat,3,nb,0, kmat)
 
!     MULTIPLY BENDING TERMS WHEN PRESENT
 
 30 IF (.NOT.bendng) GO TO 60
 DO  i = 1,3
   ii = i + 3
   DO  j = 1,3
     jj = j + 3
     g1(i,j) = gmat(ii,jj)
   END DO
 END DO
 CALL gmmats (g1,3,3,0, bmat(nd3+1),3,nb,0, gbmat)
 CALL gmmats (bmat(nd3+1),3,nb,-1, gbmat,3,nb,0, kmat)
 
 DO  i = 1,3
   ii = i + 6
   DO  j = 1,3
     jj = j + 6
     g1(i,j) = gmat(ii,jj)
   END DO
 END DO
 CALL gmmats (g1,3,3,0, bmat(nd6+1),3,nb,0, gbmat)
 CALL gmmats (bmat(nd6+1),3,nb,-1, gbmat,3,nb,0, kmat)
 60 RETURN
 
 
 ENTRY t3bgbd (ng,nb,gmad,bmad,kmad)
!     ===================================
 
!     DOUBLE PRECISION
 
 nd3 = nb*3
 nd6 = nb*6
 
!     IF [G] IS FULLY POPULATED, PERFORM STRAIGHT MULTIPLICATION AND
!     RETURN.
 
 IF (.NOT.mbcoup) GO TO 100
 CALL gmmatd (gmad,ng,ng,0,  bmad,ng,nb,0,  gbmad)
 CALL gmmatd (bmad,ng,nb,-1, gbmad,ng,nb,0, kmad )
 GO TO 150
 
!     MULTIPLY MEMBRANE TERMS WHEN PRESENT
 
 100 IF (.NOT.membrn) GO TO 120
 DO  i = 1,3
   DO  j = 1,3
     g2(i,j) = gmad(i,j)
   END DO
 END DO
 CALL gmmatd (g2,3,3,0, bmad(1),3,nb,0, gbmad)
 CALL gmmatd (bmad(1),3,nb,-1, gbmad,3,nb,0, kmad)
 
!     MULTIPLY BENDING TERMS WHEN PRESENT
 
 120 IF (.NOT.bendng) GO TO 150
 DO  i = 1,3
   ii = i + 3
   DO  j = 1,3
     jj = j + 3
     g2(i,j) = gmad(ii,jj)
   END DO
 END DO
 CALL gmmatd (g2,3,3,0, bmad(nd3+1),3,nb,0, gbmad)
 CALL gmmatd (bmad(nd3+1),3,nb,-1, gbmad,3,nb,0, kmad)
 
 DO  i = 1,3
   ii = i + 6
   DO  j = 1,3
     jj = j + 6
     g2(i,j) = gmad(ii,jj)
   END DO
 END DO
 CALL gmmatd (g2,3,3,0, bmad(nd6+1),3,nb,0, gbmad)
 CALL gmmatd (bmad(nd6+1),3,nb,-1, gbmad,3,nb,0, kmad)
 
 150 RETURN
 
END SUBROUTINE t3bgbs
