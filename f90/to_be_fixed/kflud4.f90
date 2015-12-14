SUBROUTINE kflud4
     
!     THIS ROUTINE IS USED FOR THE 4-SIDED FLUID ELEMENT. IT REARRANGES
!     THE DATA AND  CALLS THE KFLUD3 ROUTINE FOR EACH SUBELEMENT.
 
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
 
 LOGICAL :: nogo
 INTEGER :: out,necpt(100)
 REAL :: ki
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /system/ sysbuf,out,nogo,skip(34),iaxif
 COMMON /sma1dp/ r(3),z(3),nneg,nj,nptj,ki
 COMMON /sma1cl/ iopt1,k1ggsw,npvt
 COMMON /sma1et/ ecpt(100)
 EQUIVALENCE     (ecpt(1),necpt(1))
 
 IF (necpt(6) <= 0.0) RETURN
 
!     TEST FOR INTERIOR ANGLES GREATER THAN 180 DEGREES
 
 nneg = 0
 ip   = 0
 DO  i = 1,4
   DO  j = 1,3
     nj   = i + j - 1
     IF (nj > 4) nj = nj - 4
     nptj = 4*(nj-1) + 10
     r(j) = ecpt(nptj  )
     z(j) = ecpt(nptj+1)
   END DO
   IF (npvt == necpt(i+1)) ip = ip + 1
   ki   = (r(2)-r(1))*(z(3)-z(1)) - (r(3)-r(1))*(z(2)-z(1))
   IF (ki < 0) THEN
     GO TO    25
   ELSE IF (ki == 0) THEN
     GO TO  2000
   ELSE
     GO TO    30
   END IF
   25 nneg = nneg + 1
 END DO
 IF (nneg == 1 .OR. nneg == 3) GO TO 2000
 IF (ip /= 1) GO TO 2000
 ecpt(6) = ecpt(6)*2.0
 DO  i = 1,24
   ecpt(i+50) = ecpt(i)
 END DO
 DO  i = 5,24
   ecpt(i) = ecpt(i+1)
 END DO
 iret = 1
 GO TO 100
 70 ecpt( 4) = ecpt(55)
 ecpt(17) = ecpt(72)
 ecpt(18) = ecpt(73)
 iret = 2
 GO TO 100
 80 ecpt(13) = ecpt(68)
 ecpt(14) = ecpt(69)
 ecpt( 3) = ecpt(54)
 iret = 3
 GO TO 100
 90 ecpt( 9) = ecpt(64)
 ecpt(10) = ecpt(65)
 ecpt( 2) = ecpt(53)
 iret = 4
 
 100 IF (necpt(2) /= npvt .AND. necpt(3) /= npvt .AND.  &
     necpt(4) /= npvt)  GO TO 110
 CALL kflud3
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
 
 2000 CONTINUE
 nj = necpt(1)
 IF (iaxif == 0) GO TO 2001
 nj = nj/1000
 2001 CONTINUE
 WRITE  (out,3000) ufm,nj
 3000 FORMAT (a23,' 5002, INTERIOR ANGLE GREATER THAN OR EQUAL TO 180 ',  &
     'DEGREES FOR ELEMENT',i12)
 nogo = .true.
 RETURN
END SUBROUTINE kflud4
