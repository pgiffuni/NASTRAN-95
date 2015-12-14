SUBROUTINE selas2
!*****
! THIS ROUTINE IS PHASE II OF STRESS DATA RECOVERY FOR THE SCALAR SPRING
! ELEMENTS ELAS1, ELAS2, ELAS3 AND ELAS4.
!*****
 
 
 
 
! SDR2 VARIABLE CORE
 
 COMMON   /zzzzzz/  zz(1)
 
! BLOCK FOR POINTERS, LOADING TEMPERATURE AND ELEMENT DEFORMATION.
 
 COMMON   /sdr2x4/ dummy(33)          ,icstm  &
     ,                  ncstm              ,ivec  &
     ,                  ivecn              ,templd ,                  eldefm
 
! SDR2 INPUT AND OUTPUT BLOCK
 
 COMMON   /sdr2x7/ jelid              ,isilno(2)  &
     ,                  stiff              ,scoeff  &
     ,                  xxxxxx(95)  &
     ,                  jselid             ,stress  &
     ,                  yyyyyy(98) ,                  jfelid             ,force  &
     ,                  zzzzzz(23)
 EQUIVALENCE (scoeff,icoeff)
 
 
 
 idisp = ivec - 1
 disp1 = 0.0
 disp2 = 0.0
 IF (isilno(1) <= 0) GO TO 10
 iu = idisp + isilno(1)
 disp1 = zz(iu)
 10 IF (isilno(2) <= 0) GO TO 20
 iu = idisp + isilno(2)
 disp2 = zz(iu)
 20 jfelid = jelid
 force = stiff * (disp1 - disp2)
 IF (icoeff == (-1)) RETURN
 stress = scoeff * force
 jselid = jelid
 RETURN
END SUBROUTINE selas2
