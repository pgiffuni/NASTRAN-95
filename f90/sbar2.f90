SUBROUTINE sbar2( ti )
!******
! THIS ROUTINE IS THE PHASE II SUBROUTINE OF STRESS DATA RECOVERY FOR
! THE BEAM ELEMENT.
!******
  
 REAL, INTENT(IN OUT)                     :: ti(14)
 REAL :: i1       ,i2       ,l        ,m1a      ,m2a      ,m1b,  &
         m2b      ,i12      ,k1a      ,k2a      ,k1b      ,k2b , frlast(2)
 INTEGER :: tloads   ,eject    ,ished(7)
 
 
 COMMON   /system/  ibfsz    ,nout     ,idm(9)   ,line
 
! SDR2 VARIABLE CORE
 
 COMMON   /zzzzzz/  zz(1)
 
! BLOCK FOR POINTERS, LOADING TEMPERATURE AND ELEMENT DEFORMATION.
 
 COMMON   /sdr2x4/ xxxxxx(33)         ,icstm  &
     ,                  ncstm              ,ivec  &
     ,                  ivecn              ,ldtemp  &
     ,                  eldefm             ,dum8(8),tloads
 
! THE FIRST 100 LOCATIONS OF THE SDR2X7 BLOCK ARE RESERVED FOR INPUT
! PARAMETERS, THE SECOND 100 FOR STRESS OUTPUT PARAMETERS, AND FORCE
! OUTPUT PARAMETERS BEGIN AT LOCATION 201.
 
 COMMON   /sdr2x7/ jelid              ,jsilno(2)  &
     ,                  sa(36)             ,sb(36)  &
     ,                  st                 ,sdelta  &
     ,                  a                  ,fj  &
     ,                  i1                 ,i2  &
     ,                  i12                ,c1  &
     ,                  c2                 ,d1  &
     ,                  d2                 ,f1  &
     ,                  f2                 ,g1  &
     ,                  g2                 ,t_sub_0  &
     ,                  sigmat             ,sigmac  &
     ,                  l                  ,therm(6)
 
!     THERM ACTUALLY HAS 30 VALUES
 
 COMMON   /sdr2x7/ iselid             ,sig1a  &
     ,                  sig2a              ,sig3a  &
     ,                  sig4a              ,sigax  &
     ,                  sigamx             ,sigamn  &
     ,                  msten              ,sig1b  &
     ,                  sig2b              ,sig3b  &
     ,                  sig4b              ,sigbmx  &
     ,                  sigbmn             ,mscom ,                  yyyyyy(84)
 COMMON   /sdr2x7/ ifelid             ,m1a  &
     ,                  m2a                ,m1b  &
     ,                  m2b                ,v1  &
     ,                  v2                 ,fx ,                  t
 
! SDR2 SCRATCH BLOCK
 
 COMMON   /sdr2x8/ fa(6)              ,fb(6)  &
     ,                  idisp              ,iua  &
     ,                  iub                ,p1  &
     ,                  k1a                ,k2a  &
     ,                  k1b                ,k2b  &
     ,                  q                  ,w  &
     ,                  cfa(6)             ,cfb(6)  &
     ,                  cfrvec(10)         ,frvec(10)
 
!  STRESS/FORCE PRECISION CHECK
 
 COMMON   /sdr2x9/ nchk               ,isub  &
     ,                  ild                ,frtmei(2)  &
     ,                  twotop             ,fnchk
 
 EQUIVALENCE (ldtemp,templd)    ,(msten,smten)  &
     ,                  (mscom,smcom)      ,(ished(6),frlast(1))  &
     ,                  (ieid,cfrvec(1))  &
     ,                  (ished(1),lsub)    ,(ished(2),lld)
 
 DATA lld, lsub, frlast / 2*-100, -1.0E30, -1.0E30 /
 
 idisp = ivec - 1
 iua = idisp + jsilno(1)
 CALL smmats (sa(1),6,6,0, zz(iua),6,1,0,  fa,cfa )
 iub = idisp + jsilno(2)
 CALL smmats (sb(1),6,6,0, zz(iub),6,1,0,  fb,cfb )
 p1  =  fa(1) + fb(1)
 v1  = -fa(2) - fb(2)
 v2  = -fa(3) - fb(3)
 t   = -fa(4) - fb(4)
 m2a =  fa(5) + fb(5)
 m1a = -fa(6) - fb(6)
 fx  = -p1 - sdelta * eldefm
 cfrvec(2) = cfa(6) + cfb(6)
 cfrvec(3) = cfa(5) + cfb(5)
 cfrvec(9) = cfa(4) + cfb(4)
 cfrvec(7) = cfa(3) + cfb(3)
 cfrvec(6) = cfa(2) + cfb(2)
 cfrvec(8) = cfa(1) + cfb(1)
 
! IF LDTEMP = -1, THE LOADING TEMPERATURE IS UNDEFINED
 
 IF( tloads == 0 ) GO TO 10
 tsave = ti(2)
 ti(2) = (ti(1) + ti(2))/2.0  -  t_sub_0
 CALL gmmats( therm,6,5,0,  ti(2),5,1,0,  fa(1) )
 ti(2) = tsave
 fx = fx - fa(1)
 v1 = v1 - fa(2)
 v2 = v2 - fa(3)
 t  = t - fa(4)
 m2a = m2a + fa(5)
 m1a = m1a - fa(6)
 10 m1b = m1a - v1*l
 m2b = m2a - v2*l
 cfrvec(4) = cfrvec(2) + cfrvec(6) * l
 cfrvec(5) = cfrvec(3) + cfrvec(7) * l
 frvec(2) = m1a
 frvec(3) = m2a
 frvec(4) = m1b
 frvec(5) = m2b
 frvec(6) = v1
 frvec(7) = v2
 frvec(8) = fx
 frvec(9) = t
!*****
! COMPUTE ELEMENT STRESSES AT 4 POINTS
!*****
 
! COMPUTE K1A AND K2A
 
 IF (i12 /= 0.0) GO TO 30
 IF (i1 /= 0.0) GO TO 20
 k1a = 0.0
 GO TO 40
 20 k1a = -m1a / i1
 GO TO 40
 30 k1a = (m2a * i12  -  m1a * i2) / (i1 * i2  -  i12**2)
 k2a = (m1a * i12  -  m2a * i1) / (i1 * i2  -  i12**2)
 GO TO 60
 40 IF (i2 /= 0.0) GO TO 50
 k2a = 0.0
 GO TO 60
 50 k2a = -m2a / i2
 
! COMPUTE SIG1A, SIG2A, SIG3A AND SIG4A
 
 60 sig1a = k1a * c1  +  k2a * c2
 sig2a = k1a * d1  +  k2a * d2
 sig3a = k1a * f1  +  k2a * f2
 sig4a = k1a * g1  +  k2a * g2
 
! COMPUTE K1B AND K2B
 
 IF (i12 /= 0.0) GO TO 80
 IF (i1 /= 0.0) GO TO 70
 k1b = 0.0
 GO TO 90
 70 k1b = -m1b / i1
 GO TO 90
 80 k1b = (m2b * i12  -  m1b * i2) / (i1 * i2  -  i12**2)
 k2b = (m1b * i12  -  m2b * i1) / (i1 * i2  -  i12**2)
 GO TO 110
 90 IF (i2 /= 0.0) GO TO 100
 k2b = 0.0
 GO TO 110
 100 k2b = -m2b / i2
 
! COMPUTE SIG1B, SIG2B, SIG3B AND SIG4B
 
 110 sig1b = k1b * c1  +  k2b * c2
 sig2b = k1b * d1  +  k2b * d2
 sig3b = k1b * f1  +  k2b * f2
 sig4b = k1b * g1  +  k2b * g2
 IF( tloads == 0 ) GO TO 115
 
!     TEST IF AT LEAST ONE POINT TEMPERATURE IS GIVEN
 
 DO  i = 7,14
   IF( ti(i) /= 0.0 ) GO TO 112
 END DO
 GO TO 115
 112 IF( a == 0.0 ) GO TO 115
 ealf =-st / a
 sig1a = sig1a + ealf*(ti(7) - ti(3)*c1 - ti(5)*c2 - ti(1))
 sig2a = sig2a + ealf*(ti(8) - ti(3)*d1 - ti(5)*d2 - ti(1))
 sig3a = sig3a + ealf*(ti(9) - ti(3)*f1 - ti(5)*f2 - ti(1))
 sig4a = sig4a + ealf*(ti(10) - ti(3)*g1 - ti(5)*g2 - ti(1))
 sig1b = sig1b + ealf*(ti(11) - ti(4)*c1 - ti(6)*c2 - ti(2))
 sig2b = sig2b + ealf*(ti(12) - ti(4)*d1 - ti(6)*d2 - ti(2))
 sig3b = sig3b + ealf*(ti(13) - ti(4)*f1 - ti(6)*f2 - ti(2))
 sig4b = sig4b + ealf*(ti(14) - ti(4)*g1 - ti(6)*g2 - ti(2))
 115 CONTINUE
 
! COMPUTE AXIAL STRESS
 
 cfrvec(10) = 0.0
 sigax = 0.0
 IF (a /= 0.0) sigax = fx / a
 IF (a /= 0.0) cfrvec(10) = cfrvec(8) / a
 frvec(10) = sigax
 
! COMPUTE MAXIMA AND MINIMA
 
 sigamx = sigax + AMAX1(sig1a,sig2a,sig3a,sig4a)
 sigbmx = sigax + AMAX1(sig1b,sig2b,sig3b,sig4b)
 sigamn = sigax + AMIN1(sig1a,sig2a,sig3a,sig4a)
 sigbmn = sigax + AMIN1(sig1b,sig2b,sig3b,sig4b)
 
! COMPUTE MARGIN OF SAFETY IN TENSION
 
 IF(sigmat <= 0.0)GO TO 620
 IF(AMAX1(sigamx,sigbmx) <= 0.0) GO TO 620
 q=sigmat/AMAX1(sigamx,sigbmx)
 smten=q-1.0
 GO TO 630
 620 msten=1
 
!      COMPUTE MARGIN OF SAFETY IN COMPRESSION
 
 630 IF(sigmac <= 0.0) GO TO 640
 IF(AMIN1(sigamn,sigbmn) >= 0.0) GO TO 640
 w = -sigmac/AMIN1(sigamn,sigbmn)
 smcom=w-1.0
 GO TO 150
 640 mscom=1
 150 iselid = jelid
 ifelid = jelid
 
!  . STRESS CHECK...
 
 IF (nchk <= 0) GO TO 230
 ieid = jelid
 k = 0
 CALL sdrchk (frvec(2),cfrvec(2),9,k)
 
 IF (k == 0) GO TO 230
 
!  . LIMITS EXCEEDED...
 j = 0
 IF (lsub == isub .AND. frlast(1) == frtmei(1) .AND.  &
     lld == ild  .AND. frlast(2) == frtmei(2) ) GO TO 200
 lsub = isub
 lld = ild
 frlast(1) = frtmei(1)
 frlast(2) = frtmei(2)
 j = 1
 CALL page1
 180 CALL sd2rhd (ished,j)
 line = line + 1
 WRITE(nout,190)
 190 FORMAT(7X,47HTYPE     eid    m1a    m2a    m1b    m2b     v1,5X,  &
     23HV2     fa      t     sa)
 GO TO 210
 
 200 IF (eject(2) /= 0) GO TO 180
 210 WRITE(nout,220) ieid,(cfrvec(i),i=2,10)
 220 FORMAT (1H0,7X,3HBAR,i8,9F7.1)
 
 230 CONTINUE
 RETURN
END SUBROUTINE sbar2
