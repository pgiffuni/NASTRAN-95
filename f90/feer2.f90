SUBROUTINE feer2 (iret)
     
!     FEER2 INITIALIZES THEN CALLS  SDCOMP
 
 
 INTEGER, INTENT(OUT)                     :: iret
 INTEGER :: filea     ,filel    ,fileu    ,sr1fle   ,  &
     sr2fle    ,sr3fle   ,sr4fle   ,sr5fle   ,  &
     sr6fle    ,sr7fle   ,sr8fle   ,rdp      , uprtri    ,prec
 DOUBLE PRECISION :: det       ,detc     ,mindd
 
 COMMON   /opinv /  mcblt(7)  ,mcbsma(7)
 COMMON   /sfact /  filea(7)  ,filel(7) ,fileu(7) ,isr1fl   ,  &
     isr2fl    ,nz       ,det      ,detc     , power     ,isr3fl   ,mindd    ,ichl
 COMMON   /feerxx/  dumm(12)  ,ifset
 COMMON   /feercx/  ifkaa(7)  ,ifmaa(7) ,iflelm(7),iflvec(7),  &
     sr1fle    ,sr2fle   ,sr3fle   ,sr4fle   ,  &
     sr5fle    ,sr6fle   ,sr7fle   ,sr8fle   ,  &
     dmpfle    ,nord     ,xlmbda   ,neig     ,  &
     mord      ,ibk      ,critf    ,northo   , iflrva    ,iflrvc
 COMMON   /zzzzzz/  z(1)
 COMMON   /names /  ij(8)     ,rdp      ,ik(5)    ,lowtri   , uprtri
 COMMON   /system/  ksystm(54),prec
 
 iret     = 0
 
 filea(1) = iflelm(1)
 filel(1) = iflvec(1)
 fileu(1) = sr3fle
 isr1fl   = sr4fle
 isr2fl   = sr5fle
 isr3fl   = sr6fle
 ichl     = 0
 IF (ibk == 1 .OR. ifset == 1) ichl = 1
 filea(2) = ifkaa(2)
 filea(3) = ifkaa(3)
 filea(4) = ifkaa(4)
 filea(5) = prec
 filea(6) = 0
 filea(7) = 0
 filel(5) = prec
 
!     SYMMETRIC DECOMPOSITION
 
 nz = korsz(z)
 CALL sdcomp (*30,z,z,z)
 10 filel(3) = filel(2)
 filel(4) = lowtri
 CALL wrttrl (filel)
 DO  i = 1,7
   mcblt(i) = filel(i)
 END DO
 RETURN
 
 30 iret = 1
 GO TO 10
END SUBROUTINE feer2
