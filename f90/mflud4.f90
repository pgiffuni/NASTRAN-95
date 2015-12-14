SUBROUTINE mflud4
!*****
!     THIS ROUTINE IS USED FOR THE 4-SIDED FLUID ELEMENT. IT REARRANGES
!      THE DATA AND  CALLS THE MFLUD3 ROUTINE FOR EACH SUBELEMENT.
!****
!     THE ECPT DATA FOR THE ELEMENT AND ITS SUBELEMENTS ARE
 
!        FIELD      SYMBOL(FLUID4)      SYMBOL(FLUID3)
!           1            ID                  ID
!           2            SIL1                SIL1
!           3            SIL2                SIL2
!           4            SIL3                SIL3
!           5            SIL4                RHO
!           6            RHO                 BULK
!           7            BULK                N
!           8            N                   CSF
!           9            CSF                 R1
!          10            R1                  Z1
!          11            Z1                  -
!          12            -                   CSF
!          13            CSF                 R2
!          14            R2                  Z2
!          15            Z2                  -
!          16            -                   CSF
!          17            CSF                 R3
!          18            R3                  Z3
!          19            Z3                  -
!          20            -
!          21            CSF
!          22            R4
!          23            Z4
!          24            -
!          25            -
!****
 INTEGER :: necpt(100)
 COMMON/sma2io/ dum1(10),ifmgg
 COMMON /sma2cl/    iopt1,k1ggsw,npvt
 COMMON /sma2et/    ecpt(100)
 EQUIVALENCE   (ecpt(1),necpt(1))
 IF(ecpt(7) == 0.0) GO TO 120
 ecpt(7)=ecpt(7)*2.0
 DO  i=1,24
   ecpt(i+50) =ecpt(i)
 END DO
 DO  i= 5,19
   ecpt(i)= ecpt(i+1)
 END DO
 iret =1
 GO TO 100
 70 ecpt(4) = ecpt(55)
 ecpt(17)= ecpt(72)
 ecpt(18)= ecpt(73)
 iret =2
 GO TO 100
 80 ecpt(13)= ecpt(68)
 ecpt(14)= ecpt(69)
 ecpt(3)= ecpt(54)
 iret=3
 GO TO 100
 90 ecpt(9) = ecpt(64)
 ecpt(10)= ecpt(65)
 ecpt(2)= ecpt(53)
 iret=4
!*****
 
 100 IF((necpt(2) /= npvt).AND.(necpt(3) /= npvt).AND.  &
     (necpt(4) /= npvt))  GO TO 110
!*****
 CALL mflud3
 110 SELECT CASE ( iret )
   CASE (    1)
     GO TO 70
   CASE (    2)
     GO TO 80
   CASE (    3)
     GO TO 90
   CASE (    4)
     GO TO 120
 END SELECT
 120 RETURN
END SUBROUTINE mflud4
