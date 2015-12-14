SUBROUTINE mbbslj (arg,n,bsl)
     
!     SUBROUTINE TO COMPUTE EVEN ORDERED BESSEL FUNCTIONS OF FIRST KIND
 
!     UNDERFLOW MAY OCCUR IN THIS ROUTINE. THE RESULTS ARE NOT AFFECTED
 
 
 REAL, INTENT(IN)                         :: arg
 INTEGER, INTENT(OUT)                     :: n
 REAL, INTENT(OUT)                        :: bsl(4)
 
 CHARACTER (LEN=25) :: sfm
 CHARACTER (LEN=29) :: uim
 CHARACTER (LEN=25) :: uwm
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm,uwm,uim,sfm
 COMMON /system/ sysbuf,n6
 
 DO  i = 1,20
   bsl(i) = 0.0
 END DO
 asq = arg**2
 IF (asq < 0.01) GO TO 60
 n  = AMIN1(17.0,(arg+10.0))
 f  = 2*n + 4
 bsl(n+3) = 0.0
 pf = (4.0*f*(f-1.0)/asq-(f-1.0)/f)*0.3
 IF (pf <= 1.e-08) GO TO 70
 bsl(n+2) = pf*1.e-30
 pf = 0.0
 j  = n + 1
 DO  i = 1,j
   m  = n - i + 2
   f  = 2*m + 1
   bsl(m) = ((4.*(f-1.)/asq-1./f-1./(f-2.))*bsl(m+1)-bsl(m+2)/f)* (f-2.0)
   pf = pf + 2.0*bsl(m+1)
 END DO
 pf = pf + bsl(1)
 f  = 0.0
 IF (ABS(pf) <= 1.0) GO TO 20
 f  = ABS(pf)*1.e-10
 20 n  = n + 2
 DO  i = 1,n
   IF (f >= ABS(bsl(i))) bsl(i) = 0.0
   bsl(i) = bsl(i)/pf
 END DO
 m  = n
 DO  i = 1,m
   IF (ABS(bsl(n)) > 1.0E-07) RETURN
   n  = n - 1
 END DO
 RETURN
 
 60 bsl(2) = 0.125*asq
 bsl(1) = 1.0 - 2.0*bsl(2)
 n  = 2
 GO TO 90
 
 70 CALL page2 (3)
 WRITE  (n6,80) sfm,arg
 80 FORMAT (a25,' 2435, MBBSLJ SUBROUTINE FAILED BECAUSE THE ARGUMEN',  &
     'T IS TOO LARGE FOR THE BSL ARRAY', /5X,'ARG =',1P,e13.5)
 CALL mesage (-61,0,0)
 90 RETURN
END SUBROUTINE mbbslj
