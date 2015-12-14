SUBROUTINE INT 2 a8 (*,x,a8)
     
 
 , INTENT(IN OUT)                         :: *
 REAL, INTENT(IN)                         :: x(1)
 INTEGER, INTENT(IN OUT)                  :: a8(2)
 INTEGER :: jx,       power
 REAL :: rx
 CHARACTER (LEN=1) :: a(10),    ip,       im,       ib,       pt, alp(10)
 CHARACTER (LEN=8) :: k8(1),    temp,     zero,     zerox
 CHARACTER (LEN=10) :: alp10,    temp10
 COMMON /machin/ mach
 COMMON /system/ dummy(38),nbpc,     nbpw,     ncpw
 EQUIVALENCE     (temp,temp10,a(1)), (jx,rx),  (alp10,alp(1))
 DATA ip,  im,  ib,  pt,  temp, zero, zerox, nn, ll, alp10        /  &
     '+', '-', ' ', '.', 'T',  '0',  '0.0', 0,  0,  '1234567890' /
 
!     THESE ROUTINES ENCODE AN INTEGER OR F.P. NUMBER IN X, TO AN 8-BYTE
!     BCD WORD IN A8, OR AN 8-CHARACTER WORD IN K8, LEFT ADJUSTED.
!     WITH MAXIMUM NUMBERS OF DIGITS SQUEEZED INTO THE 8-BYTE FORMAT.
 
!     ENTRY POINT    INT 2 A8  (INTEGER-BCD VERSION)
!                    INT 2 K8  (INTEGER-CHARACTER VERSION)
!                    FP  2 A8  (REAL-BCD VERSION)
!                    FP  2 K8  (REAL-CHARACTER VERSION)
 
!     WRITTEN BY G.CHAN/UNISYS IN AUG. 1985
!     PARTICULARLY FOR XREAD ROUTINE, IN SUPPORT OF ITS NEW FREE-FIELD
!     INPUT FORMAT.
!     THIS ROUTINE IS MACHINE INDEPENDENT
 
 nt = +1
 GO TO 100
 
 ENTRY int2k8 (*,x,k8)
!     =====================
 
 nt = -1
 GO TO 100
 
 ENTRY fp2a8 (*,x,a8)
!     ====================
 
 nt = +2
 GO TO 100
 
 ENTRY fp 2 k8 (*,x,k8)
!     ======================
 
 nt = -2
 
 100  INT = IABS(nt)
 DO  j = 1,8
   a(j) = ip
 END DO
 a( 9) = ib
 a(10) = ib
 IF (INT /= 1) GO TO 200
 
!     INTEGER
 
 lu = 8
 n  = 0
 rx = x(1)
 ix = IABS(jx)
 xll = FLOAT(ix) + .01
 absx = ABS(xll)
 nn = 0
 IF (jx >= 0 .AND. ix < 10**8) GO TO 140
 IF (jx < 0 .AND. ix < 10**7) GO TO 140
 RETURN 1
 140  IF (jx < 0) THEN
   GO TO   210
 ELSE IF (jx == 0) THEN
   GO TO   150
 ELSE
   GO TO   220
 END IF
 
 150  temp = zero
 GO TO 310
 160  temp = zerox
 GO TO 310
 
!     F.P. NUMBER
 
 200  absx = ABS(x(1))
 IF (absx < 1.e-20) GO TO 160
 absx = absx*(1.0+1.e-20)
 lu = 7
 ll =-3
 n  = 0
 IF (x(1) > 0.) GO TO 220
 lu = lu - 1
 ll = ll + 1
 210  n  = 1
 a(1) = im
 220  n1 = n
 IF (INT == 1) GO TO 240
 xll = ALOG10(absx)
 IF (xll < 0.) xll = xll - .99998
 IF (xll > 0.) xll = xll + .00002
 power = IFIX(xll)
 np1 = power + 1
 ip1 = IABS(np1)
 xlu = 10.**lu
 xll = 10.**ll
 IF (absx < xll .OR. absx > xlu) GO TO 400
 
!     F.P. NUMBER IS SQUEEZED INTO AN EIGHT DIGIT F FORMAT, IF
!     X IS BETWWEN 10**-3 AND 10**7 AND X IS POSITUVE, OR
!          BETWWEN 10**-2 AND 10**6 AND X IS NEGATIVE,
 
 230  IF (ip1 >= 10) lu = lu - 1
 IF (np1 == -1) lu = lu + 1
 nn = lu - np1
 IF (INT == 2 .AND. nn > 7) nn = 7
 ix = IFIX(absx*10.**nn)
 240  lu = lu - 1
 IF (lu < 0 .AND. INT == 3) GO TO 420
 IF (lu < 0 .AND.   n == 7) GO TO 260
 power = 10**lu
 IF (power == 0) power = 1
 j  = ix/power
 IF (j >= 10) GO TO 240
 ix = MOD(ix,power)
 IF (lu-nn+1 < 0) THEN
   GO TO   280
 ELSE IF (lu-nn+1 == 0) THEN
   GO TO   250
 ELSE
   GO TO   270
 END IF
 250  IF (INT == 3) GO TO 420
 260  n  = n + 1
 a(n) = pt
 IF (n >= 8) GO TO 290
 270  IF (j == 0 .AND. n <= n1) SELECT CASE ( INT )
   CASE (    1)
     GO TO 240
   CASE (    2)
     GO TO 280
   CASE (    3)
     GO TO 280
 END SELECT
 280  IF (j == 0) j = 10
 n  = n + 1
 a(n) = alp(j)
 IF (lu == 0 .AND. INT == 1) GO TO 350
 IF (n < 8) GO TO 240
 290  DO  j = 1,8
   IF (a(n) == pt) GO TO 310
   IF (a(n) /= alp(10)) GO TO 310
   a(n) = ib
   n  = n - 1
 END DO
 
 310  IF (nt < 0) THEN
   GO TO   320
 ELSE IF (nt == 0) THEN
   GO TO   440
 ELSE
   GO TO   330
 END IF
 320  k8(1) = temp
 GO TO 440
 330  IF (mach /= 4) CALL khrbc2 (temp,a8(1))
!WKBD IF (MACH .EQ. 4) A8(1) = ISWAP(TEMP10)
!     IF (NCPW .GE. 8) A8(2) = LSHIFT(A8(1),4*NBPC)
 GO TO 440
 
 350  n = n + 1
 IF (n > 8) GO TO 310
 DO  j = n,8
   a(j) = ib
 END DO
 GO TO 310
 
!     F.P. NUMBER IN .XXXXX+X, .XXXX-XX, -.XXXX-X, OR -.XXX+XX FORMS
!     FOR MAXIMUM NOS. OF DIGITS POSSIBLE IN AN A8 WROD.
 
 400  INT = 3
 n   = n + 1
 a(n)= pt
 lu  = lu - 2
 GO TO 230
 
 420  n = n + 1
 IF (np1 >= 0) a(n) = ip
 IF (np1 < 0) a(n) = im
 IF (ip1 >= 10) GO TO 430
 a(n+1) = alp(ip1)
 GO TO 310
 430  j = ip1/10
 a(n+1) = alp(j)
 j = MOD(ip1,10)
 IF (j == 0) j = 10
 a(n+2) = alp(j)
 GO TO 310
 
 440  RETURN
END SUBROUTINE INT 2 a8
