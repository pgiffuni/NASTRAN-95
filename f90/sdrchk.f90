SUBROUTINE sdrchk (forvec,cfrvec,lvec,kont)
!   THIS ROUTINE IS USED BY ELEMENT SUBROUTINES THAT DETERMINE IF THE
! REQUESTED STRESS/FORCE PRECISION IS AVAILABLE...
!-----
 
 REAL, INTENT(IN OUT)                     :: forvec(lvec)
 REAL, INTENT(IN OUT)                     :: cfrvec(lvec)
 INTEGER, INTENT(IN)                      :: lvec
 INTEGER, INTENT(OUT)                     :: kont
 
 COMMON /sdr2x9/ nchk(5),twotop,fnchk
 
 DO  i = 1,lvec
   IF (cfrvec(i) == 0.0) r = 1.0E0
   IF (cfrvec(i) /= 0.0) r = ABS (forvec(i)/cfrvec(i) )
   IF (r > 1.001) r = 1.0E0
   IF (r == 0.0) rj = twotop
   IF (r /= 0.0) rj = twotop + ALOG10 (r)
   IF (rj < 0.0) rj = 0.0
   cfrvec(i) = rj
   IF (rj < fnchk) kont = kont + 1
 END DO
!-----
 RETURN
END SUBROUTINE sdrchk
