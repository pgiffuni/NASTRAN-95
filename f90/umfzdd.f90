SUBROUTINE umfzdd
     
!     SUBROUTINE TO INITIALIZE COMMON /UMFZZZ/ USED BY UMFEDT.
 
 INTEGER :: kt,tid1,tid2,pid1,no,ll
 LOGICAL :: again,t1,t2,end1
 COMMON   /umfzzz/ again,t1,t2,kt,tid1,tid2,pid1,end1,no,ll
 
 again = .false.
 t1    = .false.
 t2    = .false.
 kt = 0
 tid1 = -1
 tid2 = -1
 pid1 = -1
 end1  = .false.
 no = 0
 ll = 0
 
 RETURN
END SUBROUTINE umfzdd
