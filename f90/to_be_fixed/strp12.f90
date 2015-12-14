SUBROUTINE strp12 (ti)
     
!     PHASE II OF STRESS DATA RECOVERY
 
 
 REAL, INTENT(IN)                         :: ti(6)
 LOGICAL :: flag
 INTEGER :: tloads
 REAL :: sdelta(3)
 DIMENSION       nsil(6),str(18),nph1ou(990),si(36),stout(68), reali(4)
 COMMON /zzzzzz/ z(1)
 COMMON /sdr2x4/ dummy(35),ivec,ivecn,ldtemp,deform,dum8(8), tloads,maxsiz
 COMMON /sdr2x7/ ph1out(990),forvec(24)
 COMMON /sdr2x8/ temp,delta,npoint,ij1,ij2,npt1,vec(5),tem,  &
     z1 ovr i,z2 ovr i,stress(18)
 EQUIVALENCE     (nsil(1),ph1out(2)),(nph1ou(1),ph1out(1)),  &
     (si(1),ph1out(19)),(ldtemp,ftemp),(f1,n1)
 
!     FIRST GET FORCE VECTOR FOR THE PLATE CONSIDERATION
 
!     M ,  M ,  M ,  V ,  V    FOR ALL SIX GRID POINTS
!      X   Y   XY   X   Y
 
!                                NPTS
!     THE 5X1 FORCE VECTOR = SUMMATION  (S )(U )   FOR EACH POINT
!                                I=1       I   I
 
 npts = 6
 DO  i = 1,24
   forvec( i) = 0.0
 END DO
 forvec( 1) = ph1out(1)
 forvec( 7) = ph1out(1)
 forvec(13) = ph1out(1)
 forvec(19) = ph1out(1)
 
!     DO 155 II = 1,4
 
 ii = 0
 17 ii = ii + 1
 IF (ii > 4) GO TO 155
 
!     ZERO OUT LOCAL STRESSES
 
 sig x  1 = 0.0
 sig y  1 = 0.0
 sig xy 1 = 0.0
 sig x  2 = 0.0
 sig y  2 = 0.0
 sig xy 2 = 0.0
 
 IF (nsil(1) == 0) GO TO 30
 
!     FORM SUMMATION
 
 DO  i = 1,6
   
!     POINTER TO DISPLACEMENT VECTOR IN VARIABLE CORE
   
   npoint = ivec + nsil(i) - 1
   
   ii1 = (ii-1)*180 + 30*i - 29
   CALL gmmats (si(ii1),5,6,0, z(npoint),6,1,0, vec(1))
   
   DO  j = 2,6
     ij = (ii-1)*6 + j
     forvec(ij) = forvec(ij) + vec(j-1)
   END DO
 END DO
 
 IF (tloads == 0) GO TO 23
 jst  = (ii-1)*3 + 738
 i1   = (ii-1)*6
 flag = .false.
 f1   = ti(6)
 IF (n1 == 1) GO TO 22
 forvec(i1+2) = forvec(i1+2) - ti(2)
 forvec(i1+3) = forvec(i1+3) - ti(3)
 forvec(i1+4) = forvec(i1+4) - ti(4)
 IF (ti(5) == 0.0 .AND. ti(6) == 0.0) flag = .true.
 GO TO 23
 22 forvec(i1+2) = forvec(i1+2) + ti(2)*ph1out(jst+1)
 forvec(i1+3) = forvec(i1+3) + ti(2)*ph1out(jst+2)
 forvec(i1+4) = forvec(i1+4) + ti(2)*ph1out(jst+3)
 IF (ti(3) == 0.0 .AND. ti(4) == 0.0) flag = .true.
 23 CONTINUE
 
!     FORCE VECTOR IS NOW COMPLETE
 
 IF (ii == 4) GO TO 24
 i1 = ii*2 + 9
 i2 = i1 + 1
 z1 ovr i = -12.0*ph1out(i1)/ph1out(7+ii)**3
 z2 ovr i = -12.0*ph1out(i2)/ph1out(7+ii)**3
 GO TO 25
 24 zi ovr i = -1.5/ph1out(17)**2
 z2 ovr i = -z1 ovr i
 25 CONTINUE
 ii1 = (ii-1)*6
 
 k1 = 0
 ASSIGN 26 TO iretrn
 GO TO 170
 
 26 sig x  1 = forvec(ii1+2)*z1 ovr i - sdelta(1)
 sig y  1 = forvec(ii1+3)*z1 ovr i - sdelta(2)
 sig xy 1 = forvec(ii1+4)*z1 ovr i - sdelta(3)
 
 k1 = 1
 ASSIGN 27 TO iretrn
 GO TO 170
 
 27 sig x  2 = forvec(ii1+2)*z2 ovr i - sdelta(1)
 sig y  2 = forvec(ii1+3)*z2 ovr i - sdelta(2)
 sig xy 2 = forvec(ii1+4)*z2 ovr i - sdelta(3)
 
 GO TO 40
 30 z1 = 0.0
 z2 = 0.0
 
 40 CONTINUE
 
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
 
 IF (nph1ou(2) == 0) GO TO 120
 
!     COMPUTE PRINCIPAL STRESSES
 
 str( 1) = ph1out(1)
 str( 2) = ph1out(ii*2+9)
 str( 3) = sig x 1
 str( 4) = sig y 1
 str( 5) = sig xy 1
 str(10) = ph1out(1)
 str(11) = ph1out(ii*2+10)
 str(12) = sig x  2
 str(13) = sig y  2
 str(14) = sig xy 2
 
 DO  i = 3,12,9
   temp     = str(i)-str(i+1)
   str(i+6) = SQRT((temp/2.0)**2+str(i+2)**2)
   delta    = (str(i)+str(i+1))/2.0
   str(i+4) = delta+str(i+6)
   str(i+5) = delta-str(i+6)
   delta    = 2.0*str(i+2)
   IF (ABS(delta) < 1.0E-15 .AND. ABS(temp) < 1.0E-15) GO TO 100
   str(i+3) = ATAN2(delta,temp)*28.6478898E0
   CYCLE
   100 str(i+3) = 0.0
 END DO
 GO TO 140
 120 DO  i = 2,18
   str( i) = 0.0
 END DO
 140 str( 1) = ph1out(1)
 str(10) = ph1out(1)
 
!     ADDITION TO ELIMINATE 2ND ELEMENT ID IN OUTPUT
 
 ijk = (ii-1)*17
 stout(ijk+1) = ph1out(1)
 DO  i = 2,9
   stout(ijk+i) = str(i)
 END DO
 DO  i = 10,17
   stout (ijk+i) = str(i+1)
 END DO
 
 GO TO 17
 155 CONTINUE
 DO  i = 1,17
   ph1out(100+i) = stout(i)
 END DO
 DO  j = 1,3
   DO  i = 1,16
     j1 = 117 + (j-1)*16 + i
     j2 = (j-1)*17 + i + 18
     ph1out(j1) = stout(j2)
   END DO
 END DO
 DO  i = 1,6
   ph1out(200+i) = forvec(i)
 END DO
 DO  i = 1,5
   ph1out(206+i) = forvec(i+ 7)
   ph1out(211+i) = forvec(i+13)
 END DO
 RETURN
 
!     INTERNAL SUBROUTINE
 
 170 IF (tloads == 0 .OR. flag) GO TO 200
 jst = 738 + (ii-1)*3
 reali(1) = ph1out(8)**3/12.0
 reali(2) = ph1out(9)**3/12.0
 reali(3) = ph1out(10)**3/12.0
 centhk   = ph1out(17)*2.0
 reali(4) = centhk**3/12.0
 IF (n1 == 1) GO TO 190
 ff = ti(k1+5) - ti(1)
 IF (ABS(ph1out(k1+9+2*ii)) <= 1.0E-07) GO TO 200
 sdelta(1) = (ph1out(jst+1)*ff + ti(2)*ph1out(k1+9+2*ii))/reali(ii)
 sdelta(2) = (ph1out(jst+2)*ff + ti(3)*ph1out(k1+9+2*ii))/reali(ii)
 sdelta(3) = (ph1out(jst+3)*ff + ti(4)*ph1out(k1+9+2*ii))/reali(ii)
 GO TO 210
 190 CONTINUE
 IF (ABS(ph1out(k1+9+2*ii)) <= 1.0E-07) GO TO 200
 ff1 = (ti(k1+3) - ph1out(k1+9+2*ii)*ti(2)-ti(1))/reali(ii)
 ff2 = (ti(k1+3) - ph1out(k1+9+2*ii)*ti(2)-ti(1))/reali(ii)
 ff3 = (ti(k1+3) - ph1out(k1+9+2*ii)*ti(2)-ti(1))/reali(ii)
 sdelta(1) = ph1out(jst+1)*ff1
 sdelta(2) = ph1out(jst+2)*ff2
 sdelta(3) = ph1out(jst+3)*ff3
 GO TO 210
 200 sdelta(1) = 0.0
 sdelta(2) = 0.0
 sdelta(3) = 0.0
 210 GO TO iretrn, (26,27)
END SUBROUTINE strp12
