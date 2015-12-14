SUBROUTINE pktrq2 (ntype)
! THIS ROUTINE CALCULATES PHASE II OUTPUT FOR PLA4
 
!     NTYPE = 1 TRI-MEMBRANE
!     NTYPE = 2 QUAD-MEMBRANE
 
!     PH1OUT CONTAINS THE FOLLOWING
!     *** NTYPE = 1 ***
!     ELEMENT ID
!     3 SILS
!     5 DUMMY-S
!     3 S ARRAYS EACH 3X3
 
!     *** NTYPE = 2 ***
!     ELEMENT ID
!     4 SILS
!     4 DUMMY-S
!     4 S ARRAYS EACH 3X3
 
 
 INTEGER, INTENT(IN)                      :: ntype
 DIMENSION nsil(4), si(36)
 
 COMMON /pla4uv/ ivec, z(24)
 COMMON /pla4es/ ph1out(300)
 COMMON /pla42s/ stress(3),vec(3),temp,delta,nsize,npoint, dum(315)
 
 EQUIVALENCE (nsil(1),ph1out(2))  &
     ,           (si(1),ph1out(10))
 
 
!                        I=NSIZE
!     STRESS VECTOR = (SUMMATION  (S ) (U ))
!                        I=1        I    I
 
 nsize = ntype + 2
 DO  i = 1,nsize
!     POINTER TO DISPLACEMENT VECTOR
   npoint = ivec + nsil(i) -1
   
   CALL gmmats( si(9*i-8),3,3,0,  z(npoint),3,1,0, vec(1))
   
   DO  j=1,3
     stress(j) = stress(j) + vec(j)
   END DO
 END DO
 
 RETURN
END SUBROUTINE pktrq2
