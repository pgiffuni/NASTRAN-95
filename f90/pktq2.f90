SUBROUTINE pktq2 (npts)
     
!     THIS ROUTINE CALCULATES PHASE II OUTPUT FOR PLA4 FOR COMBINATION
!     ELEMENTS
 
!     **** PHASE II OF STRESS DATA RECOVERY *********
 
!     NPTS = 3 IMPLIES STRIA1 OR STRIA2  (PHASE II)
!     NPTS = 4 IMPLIES SQUAD1 OR SQUAD2  (PHASE II)
 
 
 INTEGER, INTENT(IN)                      :: npts
 DIMENSION       nsil(4),nph1ou(1),si(36)
 COMMON /pla4uv/ ivec,z(24)
 COMMON /pla4es/ ph1out(300)
 COMMON /pla42s/ stress(3),temp,delta,npoint,i,j,npt1,vec(5),tem,  &
     z1ovri, z2ovri,dum1(308)
 EQUIVALENCE     (nsil(1),ph1out(2)),(nph1ou(1),ph1out(1)), (si(1),ph1out(9))
 
!     PHASE I OUTPUT FROM THE MEMBRANE IS THE FOLLOWING
!     NOTE..BEGIN = 30*NPTS+8
 
!     PH1OUT(BEGIN+ 1)               ELEMENT ID
!     PH1OUT(BEGIN+ 2 THRU BEGIN +5) 3 SILS AND DUMMY OR 4 SILS
!     PH1OUT(BEGIN+ 6 THRU BEGIN +9) DUMMY
!     PH1OUT(BEGIN+10 THRU BEGIN +9*NPTS+9) 3 OR 4 S SUB I 3X3 ARRAYS
 
 
!     FIND SIG X, SIG Y, SIG XY, FOR MEMBRANE CONSIDERATION
 
 IF (nph1ou(30*npts+10) == 0) RETURN
 
!                       I=NPTS
!     STRESS VECTOR = (SUMMATION(S )(U ))
!                       I=1       I   I
 
 DO  i = 1,npts
   
!     POINTER TO I-TH SIL IN PH1OUT
   
   npoint = 30*npts + 9 + i
   
!     POINTER TO DISPLACEMENT VECTOR IN VARIABLE CORE
   
   npoint = ivec + nph1ou(npoint) - 1
   
!     POINTER TO S SUB I 3X3
   
   npt1 = 30*npts + 9 + 9*i
   
   CALL gmmats (ph1out(npt1),3,3,0, z(npoint),3,1,0, vec(1))
   
   DO  j = 1,3
     stress(j) = stress(j) + vec(j)
   END DO
   
 END DO
 
 RETURN
END SUBROUTINE pktq2
