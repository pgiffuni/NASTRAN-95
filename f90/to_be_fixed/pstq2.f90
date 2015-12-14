SUBROUTINE pstq2 (npts)
!  THIS ROUTINE CALCULATES PHASE II OUTPUT FOR PLA3
!  FOR COMBINATION ELEMENTS
 
!     ****PHASE II OF STRESS DATA RECOVERY*********
 
!     NPTS = 3 IMPLIES STRIA1 OR STRIA2  (PHASE II)
!     NPTS = 4 IMPLIES SQUAD1 OR SQUAD2  (PHASE II)
 
 
 INTEGER, INTENT(IN)                      :: npts
 DIMENSION  nsil(4), nph1ou(2), si(36)
 
 COMMON /pla3uv/ ivec,z(24)
 COMMON /pla3es/ ph1out(200),forvec(6),dummy(94)
 COMMON /pla32s/ stress(3),temp,delta,npoint,i,j,npt1,vec(5),tem,  &
     z1ovri, z2ovri,dum1(308)
 COMMON /sout/ str(18)
 EQUIVALENCE (nsil(1),ph1out(2))  &
     ,(nph1ou(1),ph1out(1)) ,(si(1),ph1out(9))
 
! **********************************************************************
! **********************************************************************
 
!     PHASE I OUTPUT FROM THE PLATE IS THE FOLLWOING
 
!     PH1OUT(1)                        ELEMENT ID
!     PH1OUT(2 THRU 5)                 3 SILS AND DUMMY OR 4 SILS
!     PH1OUT(6)                        I
!     PH1OUT(7 THRU 8)                 Z1 AND Z2
!     PH1OUT(9 THRU 30*NPTS+8)         3 OR 4 S SUB I  5X6 ARRAYS
 
! **********************************************************************
 
!     PHASE I OUTPUT FROM THE MEMBRANE IS THE FOLLOWING
!     NOTE..BEGIN = 30*NPTS+8
 
!     PH1OUT(BEGIN + 1)                ELEMENT ID
!     PH1OUT(BEGIN + 2 THRU BEGIN + 5) 3 SILS AND DUMMY OR 4 SILS
!     PH1OUT(BEGIN + 6)                T SUB 0
!     PH1OUT(BEGIN + 7 THRU BEGIN + 9) S SUB T  3X1 ARRAY
!     PH1OUT(BEGIN + 10 THRU BEGIN + 9*NPTS+9) 3 OR 4 S SUB I 3X3 ARRAYS
 
! **********************************************************************
! **********************************************************************
 
!     THE ABOVE ELEMENTS ARE COMPOSED OF PLATES AND MEMBRANES...
!     SOME MAY ONLY CONTAIN PLATES WHILE OTHERS MAY ONLY CONTAIN
!     MEMBRANES.
!     A CHECK FOR A ZERO FIRST SIL IN THE PHASE I OUTPUT, WHICH
!     INDICATES WHETHER ONE OR THE OTHER HAS BEEN OMITTED, IS MADE BELOW
 
 
 
!     FIRST GET FORCE VECTOR FOR THE PLATE CONSIDERATION
 
!     M ,  M ,  M  ,  V ,  V
!      X    Y    XY    X    Y
 
!                                NPTS
!     THE  5X1 FORCE VECTOR = SUMMATION  (S )(U )
!                                I=1       I   I
 
 
!     ZERO OUT LOCAL STRESSES
 
 sig x  1 = 0.0E0
 sig y  1 = 0.0E0
 sig xy 1 = 0.0E0
 sig x  2 = 0.0E0
 sig y  2 = 0.0E0
 sig xy 2 = 0.0E0
 
 IF( nsil(1) == 0 ) GO TO 30
 
!     FORM SUMMATION
 
 DO  i=1,npts
   
!     POINTER TO DISPLACEMENT VECTOR IN VARIABLE CORE
   
   npoint = ivec + nsil(i) - 1
   
   CALL gmmats( si(30*i-29),5,6,0,  z(npoint),6,1,0,  vec(1)  )
   
   DO  j=2,6
     forvec(j) = forvec(j) + vec(j-1)
   END DO
   
 END DO
 
!     FORCE VECTOR IS NOW COMPLETE
 
 z1 = ph1out(7)
 z2 = ph1out(8)
 
 z1 ovr i = - ph1out(7) / ph1out(6)
 z2 ovr i = - ph1out(8) / ph1out(6)
 
 sig x  1 = forvec(2) * z1 ovr i
 sig y  1 = forvec(3) * z1 ovr i
 sig xy 1 = forvec(4) * z1 ovr i
 sig x  2 = forvec(2) * z2 ovr i
 sig y  2 = forvec(3) * z2 ovr i
 sig xy 2 = forvec(4) * z2 ovr i
!     *******************************
 
 GO TO 40
 30 z1 = 0.0E0
 z2 = 0.0E0
 
!     FIND SIG X, SIG Y, SIG XY, FOR MEMBRANE CONSIDERATION
 40 IF( nph1ou(30*npts+10) == 0 ) GO TO 90
 
 
!                        I=NPTS
!     STRESS VECTOR = ( SUMMATION(S )(U ) )
!                        I=1       I   I
 
 DO  i=1,npts
   
!     POINTER TO I-TH SIL IN PH1OUT
   npoint = 30*npts + 9 + i
!     POINTER TO DISPLACEMENT VECTOR IN VARIABLE CORE
   npoint = ivec + nph1ou (npoint) - 1
   
!     POINTER TO S SUB I 3X3
   npt1 = 30 * npts + 9 + 9 * i
   
   CALL gmmats ( ph1out(npt1),3,3,0,  z(npoint),3,1,0,  vec(1)  )
   
   DO  j=1,3
     stress(j) = stress(j) + vec(j)
   END DO
   
 END DO
 
 
!     ADD MEMBRANE STRESSES TO PLATE STRESSES
 
 sig x  1 = sig x  1 + stress(1)
 sig y  1 = sig y  1 + stress(2)
 sig xy 1 = sig xy 1 + stress(3)
 sig x  2 = sig x  2 + stress(1)
 sig y  2 = sig y  2 + stress(2)
 sig xy 2 = sig xy 2 + stress(3)
 
!     STRESS OUTPUT VECTOR IS THE FOLLOWING
 
!      1) ELEMENT ID
!      2) Z1 = FIBER DISTANCE 1
!      3) SIG X  1
!      4) SIG Y  1
!      5) SIG XY 1
!      6) ANGLE OF ZERO SHEAR AT Z1
!      7) SIG P1 AT Z1
!      8) SIG P2 AT Z1
!      9) TAU MAX = MAXIMUM SHEAR STRESS AT Z1
 
!     10) ELEMENT ID
!     11) Z2 = FIBER DISTANCE 2
!     12) SIG X  2
!     13) SIG Y  2
!     14) SIG XY 2
!     15) ANGLE OF ZERO SHEAR AT Z2
!     16) SIG P1 AT Z2
!     17) SIG P2 AT Z2
!     S7) SIG P2 AT Z2
!     18) TAU MAX = MAXIMUM SHEAR STRESS AT Z2
 
 
 90 IF( nph1ou(2) == 0 .AND. nph1ou(30*npts+10) == 0 ) GO TO 120
 
!     COMPUTE PRINCIPAL STRESSES
 
 str( 1) = ph1out(1)
 str( 2) = z1
 str( 3) = sig x  1
 str( 4) = sig y  1
 str( 5) = sig xy 1
 str(10) = ph1out(1)
 str(11) = z2
 str(12) = sig x  2
 str(13) = sig y  2
 str(14) = sig xy 2
 
 DO  i=3,12,9
   temp = str(i) - str(i+1)
   str(i+6) = SQRT( (temp/2.0E0)**2 + str(i+2)**2 )
   delta = (  str(i)  +  str(i+1)  )  /  2.0E0
   str(i+4) = delta + str(i+6)
   str(i+5) = delta - str(i+6)
   delta = 2.0E0 * str(i+2)
   IF( ABS(delta) < 1.0E-15 .AND. ABS(temp) < 1.0E-15)GO TO 100
   str(i+3) = ATAN2( delta,temp ) * 28.6478898E0
   CYCLE
   100 str(i+3) = 0.0E0
 END DO
 
 GO TO 140
 120 DO  i=2,18
   str(i) = 0.0E0
 END DO
 140 str(1) = ph1out(1)
 str(10) = ph1out(1)
 
 
!     ADDITION TO ELIMINATE 2ND ELEMENT ID IN OUTPUT
 
 DO  i=10,17
   str(i) = str(i+1)
 END DO
 
 RETURN
END SUBROUTINE pstq2
