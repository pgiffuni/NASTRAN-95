SUBROUTINE sqdm12
     
!      PHASE TWO STRESS DATA RECOVERY QUADRILATERAL MEMBRANE
 
!      ELEMENT ID
!      4 SILS
!      T SUB 0
!      S SUB T 3X1
!      4 S ARRAYS EACH 3X3
 
 
!     STRES(1) - PH1OUT(1)
!     STRES(2) - SIGMA X
!     STRES(3) - SIGMA Y
!     STRES(4) - SIGMA XY
!     STRES(5) - PHI 1 ANGLE OF PRINCIPAL DIRECTION OF STRESS
!     STRES(6) - SIGMA 1
!     STRES(7) - SIGMA 2
!     STRES(8) - TAU MAXIMUM SHEAR STRESS
 
 DIMENSION nsil(4),s(36),st(3),frlast(2),ph1out(45)
 INTEGER :: eject,ished(7),istyp(2)
 
 COMMON   /system/  ibfsz    ,nout     ,idm(9)   ,line
 COMMON   /zzzzzz/ z(1)
 COMMON   /sdr2x4/ dummy(35),ivec,ivecn,ldtemp,deform
 COMMON   /sdr2x7/ est(100),stres(100),forvec(25)
 COMMON   /sdr2x8/ stress(3),vec(3),tem,temp,npoint,delta,nsize,  &
     cstrs(4),cvc(3)
 COMMON /sdr2x9/ nchk,isub,ild,frtmei(2),twotop,fnchk
 
 EQUIVALENCE (ph1out(1),est(1)) , (nsil(1),ph1out(2)) ,  &
     (tsub0,ph1out(6)) , (st(1),ph1out(7)) , (s(1),ph1out(10))  &
     , (ftemp,ldtemp) , (ished(1),lsub) , (ished(2),lld) , (ished(6),frlast(1))
 DATA istyp / 4HQDME, 2HM1 /
 DATA lsub,lld,frlast / 2*-1, -1.0E30, -1.0E30 /
 
!      ZERO OUT THE STRESS VECTOR
 
 stress(1)=0.
 stress(2)=0.
 stress(3)=0.
 cstrs(2) = 0.0E0
 cstrs(3) = 0.0E0
 cstrs(4) = 0.0E0
 
!                           I=4                      -
!         STRESS VECTOR =(SUMMATION (S )(U )) - (S )(T - T)
!                           I=1       I   I       T       0
 DO  i=1,4
   npoint=ivec+nsil(i)-1
   CALL smmats (s(9*i-8),3,3,0, z(npoint),3,1,0, vec(1),cvc(1))
   DO  j=1,3
     IF (nchk <= 0)GO TO 1
     cstrs(j+1) = cstrs(j+1) + cvc(j)
     1 stress(j) = stress(j) + vec(j)
   END DO
 END DO
 stres(1) = ph1out(1)
 stres(2) = stress(1)
 stres(3) = stress(2)
 stres(4) = stress(3)
 cstrs(1) = stres(1)
 
!      ADD IN TEMPERATURE EFFECTS
 
 IF(ldtemp == (-1)) GO TO 200
 tem = ftemp-tsub0
 DO  i=2,4
   stres(i)=stres(i)-st(i-1)*tem
 END DO
 
!      STRESS VECTOR COMPLETE AND CONTAINS SIGMA X ,  SIGMA Y ,  SIGMA X
 
!      PRINCIPAL STRESSES AND ANGLE OF ACTION PHI
 
 200 temp=stres(2)-stres(3)
 
!     COMPUTE TAU
 
 stres(8)=SQRT((temp/2.0E0)**2+stres(4)**2)
 delta=(stres(2)+stres(3))/2.0E0
 
!     COMPUTE SIGMA 1 AND SIGMA 2
 
 stres(6)=delta+stres(8)
 stres(7)=delta-stres(8)
 delta=2.0E0*stres(4)
 IF (ABS(delta) < 1.0E-15.AND.ABS(temp) < 1.0E-15) GO TO 5
 IF(ABS(temp) < 1.0E-15) GO TO 6
 
!     COMPUTE PHI 1 DEPENDING ON WHETHER OR NOT SIGMA XY AND/OR
!               (SIGMA 1 - SIGMA 2) ARE ZERO
 
 stres(5)=ATAN2(delta,temp)*28.6478898E00
 GO TO 7
 5 stres(5)=0.0E0
 GO TO 7
 6 stres(5)=45.
 7 IF (nchk <= 0) GO TO 150
 
!  . STRESS PRECISION CHECK...
 
 k = 0
 CALL sdrchk (stres(2),cstrs(2),3,k)
 IF (k == 0) GO TO 150
 
!  . LIMITS EXCEEDED...
 j = 0
 IF (lsub == isub .AND. frlast(1) == frtmei(1) .AND.  &
     lld  == ild  .AND. frlast(2) == frtmei(2)) GO TO 120
 
 lsub = isub
 frlast(1) = frtmei(1)
 frlast(2) = frtmei(2)
 lld = ild
 j = 1
 CALL page1
 100 CALL sd2rhd (ished,j)
 WRITE(nout,110)
 line = line + 1
 110 FORMAT (7X,4HTYPE,5X,3HEID,5X,2HSX,5X,2HSY,4X,3HSXY)
 GO TO 130
 120 IF (eject (2) /= 0) GO TO 100
 
 130 WRITE(nout,140) istyp,cstrs
 140 FORMAT (1H0,5X,a4,a2,i7,4F7.1)
 
 150 CONTINUE
 RETURN
END SUBROUTINE sqdm12
