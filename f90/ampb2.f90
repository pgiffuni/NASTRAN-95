SUBROUTINE ampb2(a,a11,a12,a21,a22,rp,cp,n1,n2)
     
    !     THIS SUBROUTINE IS A GENERAL DRIVER FOR PARTN
 
 
    INTEGER, INTENT(IN)                      :: a
    INTEGER, INTENT(IN)                      :: a11
    INTEGER, INTENT(IN)                      :: a12
    INTEGER, INTENT(IN)                      :: a21
    INTEGER, INTENT(IN)                      :: a22
    INTEGER, INTENT(IN)                      :: rp
    INTEGER, INTENT(IN)                      :: cp
    INTEGER, INTENT(IN OUT)                  :: n1
    INTEGER, INTENT(IN OUT)                  :: n2
    INTEGER :: rule,mcb(20),mcb1(20)
 
    COMMON /parmeg/mcba(7),mcba11(7),mcba21(7),mcba12(7),mcba22(7), nx,rule
    COMMON /zzzzzz/ iz(1)
 
    !-----------------------------------------------------------------------
 
    mcb(1)=rp
    CALL rdtrl(mcb)
    mcb1(1)=cp
    CALL rdtrl(mcb1)
    nx=korsz(iz)
    rule=0
    mcba11(1)=a11
    IF(a11 == 0)GO TO 10
    CALL rdtrl(mcba11)
    IF(mcba11(1) <= 0)mcba11(1)=0
10 CONTINUE
   mcba21(1)=a21
   IF(a21 <= 0)GO TO 20
   CALL rdtrl(mcba21)
   IF(mcba21(1) <= 0)mcba21(1)=0
20 CONTINUE
   mcba12(1)=a12
   IF(a12 == 0)GO TO 30
   CALL rdtrl(mcba12)
   IF(mcba12(1) <= 0)mcba12(1)=0
30 CONTINUE
   mcba22(1)=a22
   IF(a22 == 0)GO TO 40
   CALL rdtrl(mcba22)
   IF(mcba22(1) <= 0)mcba22(1)=0
40 CONTINUE
   mcba(1)=a
   CALL rdtrl(mcba)
   mcba11(2) = mcba(2) - mcb(6)
   mcba11(3) = mcba(3) -mcb1(6)
   mcba12(2) = mcba(2) -mcba11(2)
   mcba12(3) = mcba11(3)
   mcba21(2) = mcba11(2)
   mcba21(3) = mcba(3) -mcba11(3)
   mcba22(2) = mcb(6)
   mcba22(3) = mcb1(6)
   mcba11(4)=2
   mcba21(4)=2
   mcba12(4)=2
   mcba22(4)=2
   mcba11(5)=mcba(5)
   mcba21(5)=mcba(5)
   mcba12(5)=mcba(5)
   mcba22(5)=mcba(5)
   CALL partn(mcb,mcb1,iz)
   IF(mcba11(1) > 0)CALL wrttrl(mcba11)
   IF(mcba21(1) > 0)CALL wrttrl(mcba21)
   IF(mcba12(1) > 0)CALL wrttrl(mcba12)
   IF(mcba22(1) > 0)CALL wrttrl(mcba22)
   RETURN
END SUBROUTINE ampb2
