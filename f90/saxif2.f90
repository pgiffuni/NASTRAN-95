SUBROUTINE saxif2 (iopt,ipart,branch,eigen)
     
!     THIS ROUTINE CALCULATES FLUID VELOCITIES DUE TO HARMONIC
!     PRESSURES IN AN AXISYMMETRIC FLUID
 
!     THE OPTIONS FOR IOPT ARE
!         IOPT    ELEMENT
!           0     CAXIF2
!           1     CAXIF3
!           2     CAXIF4
!     IPART-  FIRST = 1,   SECOND = 2
!     BRANCH-  SDR2 PROCESS CODE WORD
 
 
 INTEGER, INTENT(IN)                      :: iopt
 INTEGER, INTENT(IN OUT)                  :: ipart
 INTEGER, INTENT(IN OUT)                  :: branch
 REAL, INTENT(IN)                         :: eigen(3)
 INTEGER :: sil
 
 COMMON /condas/ consts(5)
 COMMON /zzzzzz/ zz(1)
 COMMON /sdr2x4/ dumy(35)  ,ivec
 COMMON /sdr2x7/ ide       ,sil(4)   ,sv(95)    ,id1      ,  &
     velr(11)  ,id2      ,veli(11)
 EQUIVALENCE     (consts(2),twopi)
 
 
 IF (ipart == 2) GO TO 20
 DO  i = 1,11
   velr(i) = 0.0
   veli(i) = 0.0
 END DO
 20 x = 1.0
 y = 0.0
 IF (branch == 2) x = SQRT(ABS(eigen(2)))
 IF (branch == 5) x = twopi*eigen(1)
 IF (x /= 0.0) x = 1.0/x
 IF (branch /= 9) GO TO 30
 em = eigen(2)**2 + eigen(3)**2
 IF (em == 0.0) GO TO 30
 x = eigen(2)/em
 y =-eigen(3)/em
 30 IF (ipart /= 2) GO TO 40
 em = x
 x  =-y
 y  = em
 40 id1= ide
 id2= ide
 kc = iopt + 2
 kr = 3 + 2*kc
 IF (iopt == 0) kr = 6
 DO  i = 1,kc
   k = ivec + sil(i) - 1
   IF (x == 0.0) GO TO 65
   
   DO  j = 1,kr
     ij = kc*(j-1) + i
     velr(j) = sv(ij)*zz(k)*x + velr(j)
   END DO
   65 CONTINUE
   IF (y == 0.0) CYCLE
   
   DO  j = 1,kr
     ij = kc*(j-1) + i
     veli(j) = sv(ij)*zz(k)*y + veli(j)
   END DO
 END DO
 RETURN
END SUBROUTINE saxif2
