SUBROUTINE spanl2(iarg)
!*****
! THIS ROUTINE IS PHASE II OF STRESS DATA RECOVERY FOR THE SHEAR AND
! TWIST PANEL ELEMENTS.
!*****
 
 
 
 INTEGER, INTENT(IN)                      :: iarg
 REAL :: frlast(2)
 INTEGER :: eject    ,ishd(7)  ,istyp(2) ,typ(4)   ,ifor(1)
 
 
! SDR2 VARIABLE CORE
 
 COMMON   /zzzzzz/  zz(1)
 
! BLOCK FOR POINTERS, LOADING TEMPERATURE AND ELEMENT DEFORMATION.
 
 COMMON   /sdr2x4/ dummy(33)          ,icstm  &
     ,                  ncstm              ,ivec  &
     ,                  ivecn              ,templd ,                  eldefm
 
! SDR2 PHASE II INPUT AND OUTPUT BLOCK.
 
 COMMON   /sdr2x7/ ielid              ,isilno(4)  &
     ,                  s(3,4)             ,a(2)  &
     ,                  t                  ,ratio(3)  &
     ,                  sigs               ,rq(4)  &
     ,                  rk(4)              ,xxxxxx(68)
 COMMON   /sdr2x7/ jselid             ,stres(3)  &
     ,                  yyyyyy(96)
 COMMON   /sdr2x7/ jfelid             ,forces(16)  &
     ,                  zzzzzz(8)
 
! SDR2 SCRATCH BLOCK
 
 COMMON   /sdr2x8/ s1bar              ,term  &
     ,                  tau(4)             ,idisp  &
     ,                  iu                 ,ctu(4) ,                  cfrvec(19)
 
! OUTPUT PRECISION CHECK BLOCK
 
 COMMON /sdr2x9/ nchk,isub,ild,frtmei(2),twotop,fnchk
 
 COMMON   /system/  ibfsz    ,nout     ,idm(9)   ,line
 
 EQUIVALENCE(stres(1),taumax)
 EQUIVALENCE(stres(2),tauavg)
 EQUIVALENCE(stres(3),marsaf,safmar)
 EQUIVALENCE(forces(1),ifor(1),p13)
 EQUIVALENCE(forces(2),p24)
!////////  FOLLOWING 8 FORCES MAY NOT BE EQUIVALENCED CORRECTLY YET/////
 EQUIVALENCE(forces(1),f1 )
 EQUIVALENCE(forces(2),f2 )
 EQUIVALENCE(forces(3),f3 )
 EQUIVALENCE(forces(4),f4 )
 EQUIVALENCE(forces(5),f5 )
 EQUIVALENCE(forces(6),f6 )
 EQUIVALENCE(forces(7),f7 )
 EQUIVALENCE(forces(8),f8 )
 EQUIVALENCE(forces( 9),rk1)
 EQUIVALENCE(forces(10),q1)
 EQUIVALENCE(forces(11),rk2)
 EQUIVALENCE(forces(12),q2)
 EQUIVALENCE(forces(13),rk3)
 EQUIVALENCE(forces(14),q3)
 EQUIVALENCE(forces(15),rk4)
 EQUIVALENCE(forces(16),q4)
 EQUIVALENCE(ishd(1),lsub)
 EQUIVALENCE(ishd(2),lld)
 EQUIVALENCE(ishd(6),frlast(1))
 EQUIVALENCE(cfrvec(1),ifrvec)
 
 DATA lsub,lld,frlast / 2*-1, -1.0E30, -1.0E30 /
 DATA typ / 4HSHEA,1HR, 4HTWIS,1HT /
 DATA larg / 0 /
 
 idisp = ivec - 1
 
! COMPUTE AVERAGE STRESS ALONG SIDE 1 IF WE ARE DEALING WITH A SHEAR
! PANEL OR MEAN FIBRE SHEAR STRESS IF WE HAVE A TWIST PANEL.
 
 cs1br = 0.0
 s1bar = 0.0
 DO  i = 1,4
   iu = idisp + isilno(i)
   IF (iarg == 5) iu = iu + 3
   CALL smmats(s(1,i),3,1,1,zz(iu),3,1,0,term,ctrm)
   cs1br = cs1br + ctrm
   s1bar = s1bar + term
 END DO
 
! COMPUTE STRESSES AT THE CORNERS
 
 tau(1) = ratio(1) * s1bar
 tau(2) = s1bar / ratio(1)
 tau(3) = ratio(2) * s1bar
 tau(4) = ratio(3) * s1bar
 ctu(1) = ABS (ratio(1)) * cs1br
 ctu(2) = cs1br / ABS (ratio(1))
 ctu(3) = ABS (ratio(2)) * cs1br
 ctu(4) = ABS (ratio(3)) * cs1br
 
! COMPUTE AVERAGE STRESS
 
 tauavg = 0.25 * (tau(1) + tau(2) + tau(3) + tau(4))
 cfrvec(3) = 0.25E0 * (ctu(1) + ctu(2) + ctu(3) + ctu(4) )
 
! COMPUTE MAXIMUM STRESS
 
 taumax = ABS(tau(1))
 cfrvec(2) = taumax
 DO  i = 2,4
   IF (ABS(tau(i)) > taumax) taumax = ABS(tau(i))
   IF (ctu(i) > cfrvec(2))  cfrvec(2) = ctu(i)
 END DO
 
! COMPUTE MARGIN OF SAFETY
 
 IF(sigs <= 0.0)GO TO 100
 IF(taumax == 0.0)GO TO 100
 safmar=sigs/taumax-1.0
 GO TO 101
 100 marsaf=1
 101 CONTINUE
 
! FOR A SHEAR PANEL COMPUTE LOADS, FOR A TWIST PANEL COMPUTE STRESSES.
 
 IF( iarg /= 4 ) GO TO 70
 
!     SHEAR PANEL FORCES
 
 q1 = s1bar*t /  SQRT( 1.0 + ( rq(4)/rk(1) )**2)
 q2 = s1bar * rq(1) /  SQRT( 1.0 + ( rq(4)/rk(2) )**2)
 q3 = s1bar * rq(2) /  SQRT( 1.0 + ( rq(4)/rk(3) )**2)
 q4 = s1bar * rq(3) /  SQRT( 1.0 + ( rq(4)/rk(4) )**2)
 cfrvec(13) = cs1br * ABS(t) / SQRT (1.0E0 + (rq(4)/rk(1) ) **2 )
 DO  i = 1,3
   f     = SQRT (1.0E0 + ( rq(4)/rk(i+1) ) **2  )
   forces(2*i+10) = s1bar * rq(i) / f
   cfrvec(2*i+13) = cs1br * ABS(rq(i)) / f
 END DO
 
 f     = ABS (rq(4))
 rk1 = -( q1 + q4 ) * rq(4)
 rk2 = -( q1 + q2 ) * rq(4)
 rk3 = -( q2 + q3 ) * rq(4)
 rk4 = -( q3 + q4 ) * rq(4)
 cfrvec(12) = (cfrvec(13) + cfrvec(19)) * f
 cfrvec(14) = (cfrvec(13) + cfrvec(15)) * f
 cfrvec(16) = (cfrvec(15) + cfrvec(17)) * f
 cfrvec(18) = (cfrvec(17) + cfrvec(19)) * f
 f1 = q4 * rk(4)
 f2 = q1 * rk(1)
 f5 = q2 * rk(2)
 f6 = q3 * rk(3)
 cfrvec(4) = cfrvec(19) * ABS (rk(4) )
 cfrvec(5) = cfrvec(13) * ABS (rk(1) )
 cfrvec(8) = cfrvec(15) * ABS (rk(2) )
 cfrvec(9) = cfrvec(17) * ABS (rk(3) )
 f3 = -f2
 f4 = -f5
 f7 = -f6
 f8 = -f1
 cfrvec( 6) = cfrvec(5)
 cfrvec( 7) = cfrvec(8)
 cfrvec(10) = cfrvec(9)
 cfrvec(11) = cfrvec(4)
 GO TO 80
 
!     TWIST STRESSES
 
 70 p13 = a(1) * s1bar * t
 p24 = a(2) * s1bar * t
 term = t / 6.0
 cfrvec(4) = a(1) * cs1br * t
 cfrvec(5) = a(2) * cs1br * t
 p13  = p13 * term
 p24  = p24 * term
 cfrvec(4) = ABS (cfrvec(4) * term)
 cfrvec(5) = ABS (cfrvec(5) *term)
 
! STORE ELEMENT ID IN OUTPUT SLOTS.
 
 80 jselid = ielid
 jfelid = ielid
 IF (nchk <= 0) GO TO 260
 
!  . CHECK PRECISION...
 
 k = 0
 
!  . STRESSES...
 CALL sdrchk (stres(1),cfrvec(2),2,k)
 
!  . FORCES...
 i = 16
 IF (iarg /= 4) i = 2
 CALL sdrchk (forces(1),cfrvec(4),i,k)
 IF (k == 0) GO TO 260
 
!  . LIMITS EXCEEDED...
 ifrvec = ielid
 i = 1
 IF (iarg /= 4) i = 3
 istyp(1) = typ(i)
 istyp(2) = typ(i+1)
 j = 0
 
 IF  (lsub == isub .AND. frlast(1) == frtmei(1) .AND. larg == iarg  &
     .AND. lld == ild  .AND. frlast(2) == frtmei(2) ) GO TO 230
 lsub = isub
 larg = iarg
 lld = ild
 frlast(1) = frtmei(1)
 frlast(2) = frtmei(2)
 j = 2
 CALL page1
 200 CALL sd2rhd (ishd,j)
 line = line + 1
 IF (iarg == 4) WRITE(nout,210)
 IF (iarg /= 4) WRITE(nout,220)
 210 FORMAT (7X,4HTYPE,5X,42HEID  smax  SAVE  f1-4  f1-2  f2-1  f2-3  f  &
     ,60H3-2  f3-4  f4-3  f4-1   k-1  sh12   k-2  sh23   k-3  sh34 , 9HK-4  sh41)
 220 FORMAT (7X,4HTYPE,5X,27HEID  smax  SAVE  m1-3  m2-4)
 GO TO 240
 230 IF (eject(2) /= 0) GO TO 200
 240 i = 19
 IF (iarg /= 4) i = 5
 WRITE(nout,250) istyp,(cfrvec(j),j=1,i)
 250 FORMAT (1H0,6X,a4,a1,i7,18F6.1)
 
 260 CONTINUE
 RETURN
END SUBROUTINE spanl2
