SUBROUTINE pstrq2 (ntype)
! THIS ROUTINE CALCULATES PHASE II OUTPUT FOR PLA3
 
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
 
 COMMON /pla3uv/ ivec, z(24)
 COMMON /pla3es/ ph1out(300)
 COMMON /pla32s/ stress(3),vec(3),temp,delta,nsize,npoint, dum(315)
 COMMON /sout/  stres(9)
 
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
 
 stres(1) = ph1out(1)
 stres(2) = stress(1)
 stres(3) = stress(2)
 stres(4) = stress(3)
 
!     ******************************************************************
 
!     PRINCIPAL STRESSES AND ANGLE OF ACTION PHI
 temp = stres(2) - stres(3)
 stres(8) = SQRT( (temp/2.0E0)**2 + stres(4)**2 )
 delta = (stres(2) + stres(3))/2.0E0
 stres(6) = delta + stres(8)
 stres(7) = delta - stres(8)
 delta = 2.0E0 * stres(4)
 IF( ABS(delta) < 1.0E-15 .AND. ABS(temp) < 1.0E-15)GO TO 101
 stres(5) = ATAN2( delta,temp ) * 28.6478898e00
 RETURN
 101 stres(5) = 0.0E0
 RETURN
END SUBROUTINE pstrq2
