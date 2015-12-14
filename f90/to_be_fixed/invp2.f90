SUBROUTINE invp2 (*)
     
!     INVP2 INITIALIZES THEN CALLS EITHER SDCOMP OR DECOMP DEPENDING ON
!     THE OPTION SELECTED ON THE EIGR CARD
 
 
 , INTENT(IN OUT)                         :: *
 INTEGER :: filea     ,filel    ,fileu    ,scr1     ,  &
     scr2      ,scr3     ,scr4     ,scr5     ,  &
     sr1fil    ,sr2fil   ,dum      ,scr6     , rdp       ,uprtri   ,  &
     switch    ,scr7     ,scr8     ,option   , opt2      ,prec     ,q(1)
 DOUBLE PRECISION :: det       ,detdet   ,detc     ,mindd
 COMMON   /sfact /  filea(7)  ,filel(7) ,fileu(7) ,sr1fil   ,  &
     sr2fil    ,nz       ,det      ,detc     , power     ,isr3fl   ,mindd    ,ichl
 COMMON   /invpxx/  dumm(12)  ,switch
 COMMON   /invpwx/  dum(14)   ,scr1(7)  ,scr2(7)  ,scrx     ,  &
     scrxx     ,scr3     ,scr4     ,scr5     , scr6      ,scr7     ,scr8
 COMMON   /names /  ij(8)     ,rdp      ,ik(5)    ,lowtri   , uprtri
 COMMON   /dcompx/  ia(7)     ,il(7)    ,iu(7)    ,iscr1    ,  &
     iscr2     ,iscr3    ,detdet   ,ipowr    , mz        ,mind
 COMMON   /system/  ksystm(63)
 COMMON   /reigkr/  option
 COMMON   /zzzzzz/  z(1)
 EQUIVALENCE        (q(1),z(1))
 EQUIVALENCE        (ksystm(55),prec)
 DATA      opt2  /  4HUINV/
 
 filea(1) = scr1(1)
 IF (switch == 1) GO TO 10
 filel(1) = scr2(1)
 fileu(1) = scr3
 GO TO 20
 10 filel(1) = scr7
 fileu(1) = scr8
 20 CONTINUE
 sr1fil   = scr4
 sr2fil   = scr5
 isr3fl   = scr6
 ichl     = 0
 filea(2) = dum(2)
 filea(3) = dum(3)
 filea(4) = dum(4)
 filea(5) = prec
 filea(6) = 0
 filea(7) = 0
 filel(5) = prec
 IF (option == opt2) GO TO 40
 
!     SYMMETRIC DECOMPOSITION SELECTED.
 
 nz       = korsz(z)
 CALL sdcomp (*30,z,z,z)
 filel(3) = filel(2)
 filel(4) = lowtri
 CALL wrttrl (filel)
 RETURN
 30 RETURN 1
 
!     UNSYMMETRIC DECOMPOSITION SELECTED.
 
 40 DO  i = 1,21
   ia(i) = filea(i)
 END DO
 iscr1 = scr4
 iscr2 = scr5
 iscr3 = scr6
 mz    = korsz(q)
 CALL decomp (*30,q,q,q)
 il(3) = il(2)
 il(4) = lowtri
 CALL wrttrl (il)
 iu(3) = iu(2)
 iu(4) = uprtri
 iu(5) = il(5)
 CALL wrttrl (iu)
 RETURN
END SUBROUTINE invp2
