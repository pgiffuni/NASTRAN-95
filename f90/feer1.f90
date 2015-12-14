SUBROUTINE feer1
     
!     FEER1 INITIALIZES AND CALLS SUBROUTINE ADD FOR FEER
 
 INTEGER :: filea     ,fileb    ,filec    ,filek    ,  &
     filem     ,scr1     ,typalp   ,typbta   , sqr       ,rdp
 DOUBLE PRECISION :: lambda    ,dalpha   ,dbeta
 COMMON   /feercx/ filek(7)  ,filem(7) ,scr1
 COMMON   /feerxx/ lambda
 COMMON   /saddx / nomat     ,nz       ,filea(7) ,typalp   ,  &
     dalpha(2) ,fileb(7) ,typbta   ,dbeta(2) , dum(36)   ,filec(7)
 COMMON   /zzzzzz/ z(1)
 COMMON   /names / ij(8)     ,rdp      ,ik(2)    ,sqr
 COMMON   /system/ ksystm(56)
 EQUIVALENCE       (ksystm(55),iprec)
 
!     SET UP CALL TO ADD
 
 DO   i = 1,7
   filea(i) = filem(i)
   fileb(i) = filek(i)
 END DO
 dalpha(1)= lambda
 dbeta(1) = 1.0D+0
 typalp   = iprec
 typbta   = iprec
 nz       = korsz(z)
 filec(1) = scr1
 filec(2) = filek(2)
 filec(3) = filek(3)
 filec(4) = sqr
 filec(5) = iprec
 nomat    = 2
 IF (fileb(1) == 0) nomat = 1
 CALL sadd (z,z)
 CALL wrttrl (filec)
 RETURN
END SUBROUTINE feer1
