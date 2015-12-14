SUBROUTINE merged (a11,a12,a21,a22,a,rp,cp,n1,n2)
     
 
 INTEGER, INTENT(IN)                      :: a11
 INTEGER, INTENT(IN)                      :: a12
 INTEGER, INTENT(IN)                      :: a21
 INTEGER, INTENT(IN)                      :: a22
 INTEGER, INTENT(IN)                      :: a
 INTEGER, INTENT(IN)                      :: rp
 INTEGER, INTENT(IN)                      :: cp
 INTEGER, INTENT(IN)                      :: n1
 INTEGER, INTENT(IN)                      :: n2
 INTEGER :: rule,mcb(20),mcb1(20)
 COMMON /parmeg/ mcba(7),mcba11(7),mcba21(7),mcba12(7),mcba22(7), nx,rule
 COMMON /zzzzzz/ iz(1)
 
 IF (rp /= 0) GO TO 10
 mcb(1) = 0
 mcb(2) = 1
 mcb(3) = n1
 mcb(4) = 2
 mcb(5) = 1
 GO TO 20
 
 10 mcb(1) = rp
 CALL rdtrl (mcb)
 20 IF (cp /= 0) GO TO 30
 mcb1(1) = 0
 mcb1(2) = 1
 mcb1(3) = n2
 mcb1(4) = 2
 mcb1(5) = 1
 GO TO 40
 
 30 mcb1(1) = cp
 CALL rdtrl (mcb1)
 40 nx    = korsz (iz)
 rule  = 0
 iotyp = 0
 mcba11(1) = a11
 IF (a11 == 0) GO TO 50
 CALL rdtrl (mcba11)
 IF (mcba11(1) <= 0) mcba11(1) = 0
 
 50 mcba21(1) = a21
 IF (a21 == 0) GO TO 60
 CALL rdtrl (mcba21)
 IF (mcba21(1) <= 0) mcba21(1) = 0
 
 60 mcba12(1) = a12
 IF (a12 == 0) GO TO 70
 CALL rdtrl (mcba12)
 IF (mcba12(1) <= 0) mcba12(1) = 0
 
 70 mcba22(1) = a22
 IF (a22 == 0) GO TO 80
 CALL rdtrl (mcba22)
 IF (mcba22(1) <= 0) mcba22(1) = 0
 
 80 mcba(1) = a
 mcba(2) = mcb(3)
 mcba(3) = mcb1(3)
 DO  i = 1,28,7
   IF (mcba11(i) == 0) CYCLE
   iotyp   = MAX0(iotyp,mcba11(i+4))
 END DO
 mcba(4) = 2
 mcba(5) = iotyp
 IF (mcba(2) == mcba(3)) mcba(4) = 1
 CALL merge (mcb,mcb1,iz)
 CALL wrttrl (mcba)
 RETURN
END SUBROUTINE merged
