SUBROUTINE srod2
!*****
! THIS ROUTINE IS PHASE II OF STRESS DATA RECOVERY FOR THE ROD.
!*****
 REAL :: cfrvec(4),frlast(2)
 INTEGER :: eject    ,ishd(7)  ,typ(4)
 
 COMMON   /system/  ibfsz    ,nout     ,idm(9)   ,line
 COMMON /sdr2de/ skp2de(8),ieltyp
 
! SDR2 VARIABLE CORE
 
 COMMON   /zzzzzz/  zz(1)
 
! BLOCK FOR POINTERS, LOADING TEMPERATURE AND ELEMENT DEFORMATION.
 
 COMMON   /sdr2x4/ dummy(33)          ,icstm  &
     ,                  ncstm              ,ivec  &
     ,                  ivecn              ,templd ,                  eldefm
 
! SDR2 INPUT AND OUTPUT BLOCK
 
 COMMON   /sdr2x7/ ielid              ,isilno(2)  &
     ,                  sat(3)             ,sbt(3)  &
     ,                  sar(3)             ,sbr(3)  &
     ,                  st                 ,sdelta  &
     ,                  area               ,fjovrc  &
     ,                  t_subc_0           ,sigmat  &
     ,                  sigmac             ,sigmas  &
     ,                  dummy2(77) ,                  jselid             ,sigma  &
     ,                  smsig              ,tau  &
     ,                  smtau              ,dummy3(95)  &
     ,                  jfelid             ,p  &
     ,                  torque             ,dummy4(22)
 
! SCRATCH BLOCK
 
 COMMON   /sdr2x8/ trana              ,tranb  &
     ,                  rota               ,rotb  &
     ,                  iuta               ,iutb  &
     ,                  iura               ,iurb  &
     ,                  ifrvec(7)          ,chkvec(4)
 
 COMMON /sdr2x9/ nchk,isub,ild,frtmei(2),twotop,fnchk
 
 EQUIVALENCE (templd,ldtemp)    ,(smsig,mssig)  &
     ,                  (smtau,mstau)  &
     ,      (cfrvec(1),csiga) , (cfrvec(2),ctau) , (cfrvec(3),cp)  &
     ,      (cfrvec(4),ctrque), (ifrvec(4),cfrvec(1))  &
     ,      (ishd(1),lsub), (ishd(2),lld), (ishd(6),frlast(1))
 
 DATA lld,lsub,frlast / 2*-1, -1.0E30, -1.0E30 /
 DATA typ / 4H con , 4HROD  , 4HTUBE, 1H   /
 
 idisp = ivec - 1
 iuta  = idisp + isilno(1)
 CALL smmats (sat(1),3,1,1, zz(iuta),3,1,0, trana,ctrna)
 iutb  = idisp + isilno(2)
 CALL smmats (sbt(1),3,1,1, zz(iutb),3,1,0, tranb,ctrnb)
 sigma = trana + tranb + sdelta * eldefm
 csiga = ctrna + ctrnb
 IF (ldtemp == (-1) ) GO TO 10
 sigma = sigma + st * (templd - t_subc_0)
 10 iura  = iuta + 3
 chkvec(1) = sigma
 CALL smmats (sar(1),3,1,1, zz(iura),3,1,0, rota,crta)
 iurb  = iutb + 3
 CALL smmats (sbr(1),3,1,1, zz(iurb),3,1,0, rotb,crtb)
 torque = rota + rotb
 cp = area * csiga
 chkvec(3) = p
 ctau = ABS (fjovrc) * ctrque
 chkvec(2) = tau
 ctrque = crta + crtb
 chkvec(4) = torque
 
! COMPUTE AXIAL FORCE, P, AND TORQUE
 
 p = area * sigma
 tau = fjovrc * torque
 
! COMPUTE MARGIN OF SAFETY IN EXTENSION
 
 IF(sigma <= 0.0)GO TO 101
 IF(sigmat <= 0.0)GO TO 102
 smsig=sigmat/sigma-1.0
 GO TO 180
 101 IF(sigma /= 0.0) GO TO 103
 GO TO 102
 103 IF(sigmac <= 0.0) GO TO 102
 smsig = -sigmac/sigma - 1.0
 GO TO 180
 102 mssig=1
 
!     COMPUTE MARGIN OF SAFETY IN TORSION
 
 180 IF(sigmas <= 0.0) GO TO 190
 IF(tau == 0.0)GO TO 190
 smtau= sigmas/ABS(tau) - 1.0
 GO TO 110
 190 mstau=1
 110 jselid = ielid
 jfelid = ielid
 IF (nchk <= 0 ) GO TO 260
 
!  . CHECK PRECISION...
 
 ifrvec(3) = ielid
 k = 0
 CALL sdrchk (chkvec,cfrvec,4,k)
 
 IF (k == 0) GO TO 260
 
!  . LIMITS EXCEEDED...
 
 j = 0
 ifrvec(1) = typ(4)
 IF (ieltyp == 10) ifrvec(1) = typ(1)
 ifrvec(2) = typ(2)
 IF (ieltyp == 3) ifrvec(2) = typ(3)
 
 IF (lsub == isub .AND. frlast(1) == frtmei(1) .AND.  &
     lld == ild  .AND. frlast(2) == frtmei(2) ) GO TO 240
 
 lsub = isub
 lld = ild
 frlast(1) = frtmei(1)
 frlast(2) = frtmei(2)
 j = 1
 CALL page1
 
 220 CALL sd2rhd (ishd,j)
 line = line + 1
 WRITE(nout,230)
 230 FORMAT(7X,4HTYPE,5X,3HEID,5X,2HSA,5X,2HST,5X,9HAF torque )
 GO TO 245
 
 240 IF (eject(2) /= 0) GO TO 220
 245 WRITE(nout,250) ifrvec
 250 FORMAT (1H0,3X,2A4,i7,4F7.1)
 260 CONTINUE
 RETURN
END SUBROUTINE srod2
