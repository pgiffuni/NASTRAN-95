SUBROUTINE strm62 (ti)
     
 
!     PHASE II OF STRESS DATA RECOVERY FOR TRIANGULAR MEMBRANE ELEMENT
!     TRIM6
 
!     PHASE I OUTPUT IS THE FOLLOWING
 
!     PH1OUT(1)               ELEMENT ID
!     PH1OUT(2, THRU 7)       6 S1L5
!     PH1OUT(8 THRU 10)       THICKNESSES AT CORNER GRID POINT
!     PH1OUT(11)              REFERENCE TEMPERATURE
!     PH1OUT(12)-(227)        S SUB I MATRICES FOR 4 POINTS
!     PH1OUT(228)-(230)       THERMAL VECTOR - G TIMES ALPHA
 
 
 
 REAL, INTENT(IN)                         :: ti(6)
 INTEGER :: tloads
 DIMENSION  ns1l(6),nph1ou(990),str(18),si(36), stout(99),stress(3)
 COMMON /zzzzzz/ z(1)
 COMMON /sdr2x4/ dummy(35),ivec,ivecn,ldtemp,deform,dum8(8),tloads
 COMMON /sdr2x7/ ph1out(250)
 COMMON /sdr2x8/ temp,delta,npoint,ij1,ij2,npt1,vec(5),tem
 EQUIVALENCE     (ns1l(1),ph1out(2)),(nph1ou(1),ph1out(1)),  &
     (si(1),ph1out(11)),(ldtemp,ftemp)
 
 DO  ii=1,4
   
!     ZERO OUT LOCAL STRESSES
   
   sig x  1 =0.0
   sig y  1 =0.0
   sig xy 1 =0.0
   sig x  2 =0.0
   sig y  2 =0.0
   sig xy 2 =0.0
   IF (ns1l(1) == 0) GO TO 90
   
!     ZERO STRESS VECTOR STORAGE
   
   DO  i=1,3
     stress(i)=0.0
   END DO
   
!                        I=6
!     STRESS VECTOR =(SUMMATION (5 )(U ) ) - (S )(TEMP      - TEMP   )
!                        I=1      I   I        T      POINT       REF
   
   DO  i=1,6
     
!     POINTER TO I-TH SIL IN PH1OUT
     
     npoint = ivec + nph1ou(i+1) - 1
     
!     POINTER TO  3X3 S SUB I MATRIX
     
     npt1=12+(i-1)*9+(ii-1)*54
     
     CALL gmmats (ph1out(npt1),3,3,0,z(npoint),3,1,0,vec(1))
     DO  j=1,3
       stress(j)=stress(j)+vec(j)
       str(j)=stress(j)
     END DO
   END DO
   IF (ldtemp == (-1)) GO TO 80
   ii12=ii*2-1
   IF (ii /= 4) tem=ti(ii12)-ph1out(11)
   IF( ii == 4) tem=(ti(1)+ti(2)+ti(3)+ti(4)+ti(5)+ti(6))/6.0- ph1out(11)
   DO  i=1,3
     stress(i)=stress(i)-ph1out(227+i)*tem
     str(i)=stress(i)
   END DO
   80 CONTINUE
   90 IF (nph1ou(2) == 0) GO TO 120
   
!     COMPUTE PRINCIPAL STRESSES
   
   
!     8 LOCATIONS FOR STRESS AT A POINT AS FOLLOWS
   
!      1. ELEMENT ID
!      2. SIGMA X1
!      3. SIGMA Y1
!      4. SIGMA XY1
!      5. ANGLE OF ZERO SHEAR
!      6. SIGMA PRINCIPAL STRESS 1
!      7. SIGMA PRINCIPAL STRESS 2
!      8. TAU MAX
   
!     FOR EACH POINT, THESE VALUES ARE STORED IN STOUT(1-8,9-16,
!     17-24,25-32) ALSO IN LOCATIONS STR(1-7) EXCEPT THE ELEMENT ID
!     FINALLY, THESE VALUES ARE STORED IN PH1OUT(101-108,109-115,
!     116-122,123-129)
   
   temp = stress(1)-stress(2)
   temp1= SQRT ((temp/2.0E0)**2 + stress(3)**2)
   str(7)= temp1
   delta= (stress(1)+stress(2))/2.0
   str(5)=delta+temp1
   str(6)=delta-temp1
   delta= 2.0E0 * stress(3)
   IF (ABS(delta) < 1.0E-15.AND.ABS(temp) < 1.0E-15) GO TO 100
   str(4)=ATAN2(delta,temp)*28.6478898E0
   GO TO 110
   100 str(4)=0.0
   110 CONTINUE
   GO TO 140
   120 DO  i=1,9
     str(i)=0.0E0
   END DO
   140 CONTINUE
   ijk=(ii-1)*8
   stout(ijk+1)=ph1out(1)
   DO  i=2,8
     stout(ijk+i)=str(i-1)
   END DO
 END DO
 DO  i=1,8
   ph1out(100+i)=stout(i)
 END DO
 DO  j=1,3
   DO  i=1,7
     j1=108+(j-1)*7+i
     j2=j*8+i+1
     ph1out(j1)=stout(j2)
   END DO
 END DO
 RETURN
END SUBROUTINE strm62
