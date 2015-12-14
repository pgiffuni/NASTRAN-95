SUBROUTINE xylog( v1, v2, cycles )
     
 
 REAL, INTENT(IN OUT)                     :: v1
 REAL, INTENT(IN OUT)                     :: v2
 INTEGER, INTENT(OUT)                     :: cycles
 INTEGER :: power1, power2
!*****
!  THIS SUBROUTINE TAKES V1 AND V2 REGARDLESS OF THEIR VALUES
!  AND COMPUTES A LOG SCALE OF AT LEAST 1 CYCLE...
!*****
 IF( v1 > 0.0E0 ) GO TO 20
 IF( v2 > 0.0E0 ) GO TO 10
 
!     V1 AND V2 ARE BOTH NEGATIVE OR ZERO.  SET ARBITRARY LIMITS
 
 5 v1 = 1.0E-5
 v2 = 1.0E+5
 cycles = 10
 RETURN
 
!     V2 IS POSITIVE BUT V1 IS NEGATIVE OR 0
 
 10 v1 = v2 * 1.0E-5
 GO TO 40
 
 20 IF( v2 > 0.0E0 ) GO TO 30
 
!     V1 IS POSITIVE BUT V2 IS NEGATIVE OR 0
 
 v2 = v1 * 1.0E+5
 GO TO 40
 
 30 IF( v2 > v1 ) GO TO 40
 temp = v1
 v1 = v2
 v2 = temp
 
!     RAISE V2 TO POWER OF 10,  LOWER V1 TO POWER OF 10
 
 40 power1 = 0
!WKBR 9/93  50 IF( V1 .LT. 1.0E0      ) GO TO 70
 50 IF( v1 < 0.99999) GO TO 70
!WKBR 9/93  60 IF( V1 .LT. 10.0E0) GO TO 80
 60 IF( v1 <= 10.0001) GO TO 80
 v1 = v1 / 10.0E0
 power1 = power1 + 1
 GO TO 60
 70 v1 = v1 * 10.0E0
 IF( v1 <= 0.0E0 ) GO TO 5
 power1 = power1 - 1
 GO TO 50
 
 80 v1 = 10.0E0 ** power1
 
 power2 = 1
 90 IF(v2 <= 1.0E0) GO TO 110
!WKBR 9/93 100 IF( V2 .LT. 10.00001E0) GO TO 120
 100 IF( v2 <= 10.0001   ) GO TO 120
 v2 = v2 / 10.0E0
 power2 = power2 + 1
 GO TO 100
 110 v2 = v2 * 10.0E0
 IF( v2 <= 0.0E0 ) GO TO 5
 power2 = power2 - 1
 GO TO 90
 
 120 v2 = 10.0 ** power2
 
 cycles = power2 - power1
 RETURN
END SUBROUTINE xylog
