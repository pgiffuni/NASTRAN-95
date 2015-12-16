SUBROUTINE cinvp2 (*)
     
!     CINVP2 INITIALIZES AND CALLS CDCOMP FOR CINVPR
 
 
 INTEGER :: filea     ,filel    ,fileu    ,scr1     ,  &
     scr2      ,scr3     ,scr4     ,scr5     ,  &
     scr6      ,scr7     ,scr8     ,switch   ,  &
     sr1fil    ,sr2fil   ,sr3fil   ,sysbuf   , dumm      ,cdp      ,scr9
 DOUBLE PRECISION :: det       ,mindia
 COMMON   /cdcmpx/  filea(7)  ,filel(7) ,fileu(7), sr1fil   ,  &
     sr2fil    ,sr3fil   ,det(2)   ,power    , nz        ,mindia   ,ib
 COMMON   /cinvxx/  dum(4)    ,switch
 COMMON   /cinvpx/  dumm(36)  ,scr1     ,scr2     ,scr3     ,  &
     scr4      ,scr5     ,scr6     ,scr7     , scr8      ,scr9
 COMMON   /zzzzzz/  z(1)
 COMMON   /names /  ij(10)    ,cdp
 COMMON / system /  sysbuf
 
 ioff     = fileu(7)
 filea(1) = scr1
 IF (switch == 0) GO TO 10
 filel(1) = scr8
 fileu(1) = scr9
 IF (switch < 0) filea(1) = -filea(1)
 IF (switch == -204) GO TO 20
 switch   = 1
 GO TO 20
 10 filel(1) = scr3
 fileu(1) = scr4
 20 sr1fil   = scr5
 sr2fil   = scr6
 sr3fil   = scr7
 filea(2) = dumm(3)
 filea(3) = dumm(3)
 filea(4) = dumm(4)
 filea(5) = cdp
 filea(6) = 0
 filea(7) = 0
 filel(5) = cdp
 nz       = korsz(z)
 IF (switch == -204) nz = nz - 2*sysbuf
 ib       = 0
 CALL cdcomp (*30,z,z,z)
 IF (switch  /= 0) fileu(7) = ioff
 RETURN
 
 30 RETURN 1
END SUBROUTINE cinvp2
