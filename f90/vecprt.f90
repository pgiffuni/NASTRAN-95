SUBROUTINE vecprt (*,*,px,nx,a,ox)
     
 INTEGER, INTENT(IN)                      :: px
 INTEGER, INTENT(IN)                      :: nx
 REAL, INTENT(IN OUT)                     :: a(nx)
 INTEGER, INTENT(IN)                      :: ox
 INTEGER :: p, o, count,eject,pm,tra,rsp,rdp,csp,cdp
 
 COMMON /system/ skp1,mo,skp2(6),maxlin,skp3(2),count
 DATA            rsp,rdp,csp,cdp / 1,2,3,4 /
 
!     PX = VECTOR TYPE + PRECISION.
!     NX = VECTOR LENGTH.
!     A  = VECTOR LOCATION.
 
!     THE VECTOR COMPONENTS WILL BE PRINTED 6 PER LINE IF REAL OR
!                IMAGINARY, AND 3 PER LINE IF COMPLEX.
!          O = 0 IF ALL THE VECTOR COMPONENTS ARE TO BE PRINTED, AND IF
!                THEY ARE TO BE PRINTED STARTING ON A NEW PAGE IF THEY
!                WILL NOT FIT ON THE CURRENT PAGE.
!          O = 1 IF ONLY THOSE LINES WHICH HAVE AT LEAST ONE NON-ZERO
!                COMPONENT ARE TO BE PRINTED, AND IF THE VECTOR IS TO BE
!                PRINTED STARTING ON A NEW PAGE IF IT WILL NOT FIT ON
!                THE CURRENT PAGE.
!          O =-1 IF ONLY THOSE LINES WHICH HAVE AT LEAST ONE NON-ZERO
!                COMPONENT ARE TO BE PRINTED, AND IF THE VECTOR IS TO BE
!                PRINTED ON THE CURRENT PAGE UNLESS TWO LINES WILL NOT
!                FIT.
 
!     RETURN 1 - PRINT SUBTITLE + VECTOR IDENTIFICATION.
!     RETURN 2 - PRINT VECTOR IDENTIFICATION ONLY.
!                PRTVEC = RETURN ENTRY POINT.
 
 p = px
 n = nx
 o = ox
 
 pm = p
 IF (p == rdp) pm = rsp
 IF (p == cdp) pm = csp
 kk = 1
 IF (pm == csp) kk = 2
 IF (p == rdp .OR. p == cdp) kk = 2*kk
 kn = kk*n
 IF (pm == csp) kk = kk/2
 k6 = kk*6
 IF (o == 0) GO TO 40
 
 m = 1
 DO  k = 1,kn,k6
   l = k + k6 - kk
   IF (l > kn)  l = kn
   DO  i = k,l,kk
     IF (a(i) /= 0.) GO TO 20
   END DO
   CYCLE
   20 m = m + 1
 END DO
 IF (m == 1) GO TO 160
 
 IF (o < 0) m = 2
 GO TO 50
 40 m = (n+5)/6 + 1
 IF (pm == csp) m = (n+2)/3 + 2
 50 ASSIGN 60 TO tra
 IF (eject(m) == 0.0) THEN
   GO TO   180
 ELSE
   GO TO   170
 END IF
 60 count = count - m
 knkk  = kn/kk
 IF (knkk > 6) GO TO 70
 CALL FORMAT (a,1,kn,kk,-1,n)
 count = count + 1
 GO TO 140
 
 70 ASSIGN 110 TO tra
 k = 1
 80 l = k + k6 - kk
 IF (l > kn) l = kn
 IF (o ==  0) GO TO 100
 DO  i = k,l,kk
   IF (a(i) /= 0.) GO TO 100
 END DO
 GO TO 130
 100 IF (eject(1) /= 0) GO TO 170
 110 k1 = (k + kk - 1)/kk
 k2 = (l + kk - 1)/kk
 IF (pm /= csp) GO TO 120
 k1 = (k1+1)/2
 k2 = k2/2
 120 CALL FORMAT (a,k,l,kk,k1,k2)
 130 k = k + k6
 IF (k <= kn) GO TO 80
 
 140 WRITE  (mo,150)
 150 FORMAT (1X)
 count = count + 1
 160 RETURN
 
 170 RETURN 1
 180 RETURN 2
 
 
 ENTRY prtvec (*,*)
!     ==================
 
 count = count + 1
 IF (pm /= csp) GO TO 260
 count = count + 1
 IF (knkk-4 < 0) THEN
   GO TO   200
 ELSE IF (knkk-4 == 0) THEN
   GO TO   220
 ELSE
   GO TO   240
 END IF
 200 WRITE  (mo,210)
 210 FORMAT (51X,4HREAL,11X,9HIMAGINARY)
 GO TO 260
 220 WRITE  (mo,230)
 230 FORMAT (21X,2(12X,4HREAL,11X,9HIMAGINARY))
 GO TO 260
 240 WRITE  (mo,250)
 250 FORMAT (3X,3(12X,4HREAL,11X,9HIMAGINARY))
 260 GO TO tra, (60,110)
END SUBROUTINE vecprt
