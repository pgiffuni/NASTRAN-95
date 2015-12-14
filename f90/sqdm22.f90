SUBROUTINE sqdm22
     
!     PHASE-II STRESS-DATA-RECOVERY ROUTINE FOR THE -QDMEM2- ELEMENT.
 
!     THIS ROUTINE USES DATA PREPARED BY -SQDM21-, THE PHASE-I ROUTINE,
!     TOGETHER WITH THE DISPLACEMENT VECTOR AND TEMPERATURE DATA
!     TO ARRIVE AT STRESS AND FORCE OUTPUTS.
 
 INTEGER :: iforce(1),istr(1),isils(4),eject,ished(7),istyp(2)
 REAL :: stress(8),force(17),kij(9,16),sg(36),pt(3,4), st(3),rg(4),frlast(2)
 COMMON /zzzzzz/ z(1)
 COMMON /sdr2x4/ dummy(35),ivec,ivecn,ldtemp,deform,dum8(8),tloads
 COMMON /sdr2x7/ id(217)
 COMMON /sdr2x8/ vec(4),sigxyz(3),f(4,3),shears(4),cvc(4),ff(4,3),  &
     cfrvec(20),cshars(4)
 COMMON /sdr2x9/ nchk,isub,ild,frtmei(2),twotop,fnchk
 COMMON /system/ ibfsz,nout,idm(9),line
 EQUIVALENCE     (temp,ldtemp)
 EQUIVALENCE     (id(2),isils(1)),(id(7),tsub0),(id(8),kij(1,1)),  &
     (id(101),stress(1),istr(1)), (id(152),sg(1)),(id(188),pt(1,1)),  &
     (id(200),st(1)),(id(201),force(1),iforce(1)), (id(203),rg(1))
 EQUIVALENCE     (f2,force( 2)),(f1,force( 3)),(f3,force( 4)),  &
     (f4,force( 5)),(f6,force( 6)),(f5,force( 7)),  &
     (f7,force( 8)),(f8,force( 9)),(fk1,force(10)),  &
     (q1,force(11)),(fk2,force(12)),(q2,force(13)),  &
     (fk3,force(14)),(q3,force(15)),(fk4,force(16)), (q4,force(17))
 EQUIVALENCE     (ished(1),lsub),(ished(2),lld), (ished(6),frlast(1))
 DATA    istyp / 4HQDME, 2HM2 /
 DATA    lsub  , lld,frlast   / 2*-1, -1.0E30,-1.0E30 /
 
!     SIG , SIG , TAU   = SUMMATION((S )(U )) - (S )(TEMP-T )
!        X     Y     XY               I   I       T        0
 
 sigxyz(1) = 0.0
 sigxyz(2) = 0.0
 sigxyz(3) = 0.0
 cfrvec(2) = 0.0
 cfrvec(3) = 0.0
 cfrvec(4) = 0.0
 
 DO  i = 1,4
   j = ivec + isils(i)
   CALL smmats (sg(9*i-8),3,3,0, z(j-1),3,1,0, vec,cvc)
   DO  j = 1,3
     sigxyz(j  ) = sigxyz(j  ) + vec(j)
     cfrvec(j+1) = cfrvec(j+1) + cvc(j)
   END DO
 END DO
 
 IF (ldtemp == -1) GO TO 40
 tbar = temp - tsub0
 DO  j = 1,3
   sigxyz(j) = sigxyz(j) - st(j)*tbar
 END DO
 
!     FORCES
!          I                             T
!        (F ) = SUMMATION((K  )(U )) - (P )(TEMP-T )
!                           IJ   I       I        0
 
 40 ipart   = 0
 DO  i = 1,4
   f(i,1)  = 0.0
   f(i,2)  = 0.0
   f(i,3)  = 0.0
   ff(i,1) = 0.0
   ff(i,2) = 0.0
   ff(i,3) = 0.0
   DO  j = 1,4
     k       = ivec  + isils(j)
     ipart   = ipart + 1
     CALL smmats (kij(1,ipart),3,3,0, z(k-1),3,1,0, vec,cvc)
     f(i,1)  = f(i,1)  + vec(1)
     f(i,2)  = f(i,2)  + vec(2)
     f(i,3)  = f(i,3)  + vec(3)
     ff(i,1) = ff(i,1) + cvc(1)
     ff(i,2) = ff(i,2) + cvc(2)
     ff(i,3) = ff(i,3) + cvc(3)
   END DO
   IF (ldtemp == -1) CYCLE
   tbar   = temp - tsub0
   f(i,1) = f(i,1) - pt(1,i)*tbar
   f(i,2) = f(i,2) - pt(2,i)*tbar
   f(i,3) = f(i,3) - pt(3,i)*tbar
 END DO
 
!     SHEARS = SUMMATION (R )(U )
!                          I   I
 DO  i = 1,4
   ip1 = i + 1
   IF (ip1 == 5) ip1 = 1
   shears(i) = (f(ip1,2)-f(i,1))/rg(i)
   cshars(i) = (ff(ip1,2)-ff(i,1))/ABS(rg(i))
 END DO
 
!     ALL COMPUTATIONS COMPLETE.
 
 q1 = -shears(1)
 q2 =  shears(2)
 q3 = -shears(3)
 q4 =  shears(4)
 cfrvec(14) = -cshars(1)
 cfrvec(16) = +cshars(2)
 cfrvec(18) = -cshars(3)
 cfrvec(20) = +cshars(4)
 
 istr(1)   = id(1)
 cfrvec(1) = stress(1)
 stress(2) = sigxyz(1)
 stress(3) = sigxyz(2)
 stress(4) = sigxyz(3)
 
 iforce(1) = id(1)
 f1 = f(1,1)
 f2 = f(1,2)
 f3 = f(2,2)
 f4 = f(2,1)
 f5 = f(3,1)
 f6 = f(3,2)
 f7 = f(4,2)
 f8 = f(4,1)
 cfrvec( 6) = ff(1,1)
 cfrvec( 5) = ff(1,2)
 cfrvec( 7) = ff(2,2)
 cfrvec( 8) = ff(2,1)
 cfrvec(10) = ff(3,1)
 cfrvec( 9) = ff(3,2)
 cfrvec(11) = ff(4,2)
 cfrvec(12) = ff(4,1)
 
 fk1 = f(1,3)
 fk2 = f(2,3)
 fk3 = f(3,3)
 fk4 = f(4,3)
 cfrvec(13) = ff(1,3)
 cfrvec(15) = ff(2,3)
 cfrvec(17) = ff(3,3)
 cfrvec(19) = ff(4,3)
 
 temp = stress(2) - stress(3)
 
!     COMPUTE TAU
 
 stress(8) = SQRT((temp/2.0)**2+stress(4)**2)
 delta = (stress(2)+stress(3))/2.0
 
!     COMPUTE SIGMA 1 AND SIGMA 2
 
 stress(6) = delta + stress(8)
 stress(7) = delta - stress(8)
 delta = 2.0*stress(4)
 
!     COMPUTE PHI 1 DEPENDING ON WHETHER OR NOT SIGMA XY AND/OR
!               (SIGMA 1 - SIGMA 2) ARE ZERO
 
 IF (ABS(temp) < 1.0E-15) GO TO 5
 stress(5) = ATAN2(delta,temp)*28.64788980
 GO TO 7
 5 IF (ABS(delta) < 1.0E-15) GO TO 6
 stress(5) = 0.0
 GO TO 7
 6 stress(5) = 45.0
 7 IF (nchk <= 0) GO TO 150
 
!     STRESS/FORCE PRECISION CHECK
 
 k = 0
 
!     STRESSES
 
 CALL sdrchk (stress(2),cfrvec(2),3,k)
 
!     FORCES
 
 CALL sdrchk (force(2),cfrvec(5),16,k)
 IF (k == 0) GO TO 150
 
!     LIMITS EXCEEDED
 
 j = 0
 IF (lsub == isub .AND. frlast(1) == frtmei(1) .AND.  &
     lld == ild  .AND. frlast(2) == frtmei(2)) GO TO 120
 
 lsub = isub
 lld  = ild
 frlast(1) = frtmei(1)
 frlast(2) = frtmei(2)
 j = 1
 CALL page1
 100 CALL sd2rhd (ished,j)
 line = line + 1
 WRITE  (nout,110)
 110 FORMAT (3X,4HTYPE,5X,3HEID,4X,2HSX,4X,2HSY,3X,3HSXY,11H  f1-4  f1-  &
     ,60H2  f2-1  f2-3  f3-2  f3-4  f4-3  f4-1   k-1  sh12   k-2  sh2  &
     ,25H3   k-3  sh34   k-4  sh41)
 GO TO 130
 120 IF (eject(2) /= 0) GO TO 100
 
 130 WRITE  (nout,140) istyp,cfrvec
 140 FORMAT (2H0 ,a4,a2,i7,19F6.1)
 
 150 RETURN
END SUBROUTINE sqdm22
