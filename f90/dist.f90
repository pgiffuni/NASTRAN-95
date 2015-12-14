SUBROUTINE dist (ideg,hist,median,modd)
     
!     COMPUTE THE DISTRIBUTION OF NODAL DEGREES WITH MEDIAN AND MODE
!     THIS ROUTINE IS USED ONLY IN BANDIT MODULE
 
 
 INTEGER, INTENT(IN)                      :: ideg(1)
 INTEGER, INTENT(OUT)                     :: hist(1)
 INTEGER, INTENT(OUT)                     :: median
 INTEGER, INTENT(OUT)                     :: modd
 
 COMMON /system/ isys,     nout
 COMMON /bands / nn,       mm
 
!     IDEG(I) = DEGREE OF NODE I
!     HIST(I) = NUMBER OF NODES OF DEGREE I
 
!     COMPUTE HISTOGRAM.
 
 mm1 = mm + 1
 DO  i = 1,mm1
   hist(i) = 0
 END DO
 DO  i = 1,nn
   k = ideg(i) + 1
   hist(k) = hist(k) + 1
 END DO
 
!     COMPUTE MODE (MODD).
 
 modd = 0
 MAX  = 0
 DO  i = 1,mm1
   k = hist(i)
   IF (k <= MAX) CYCLE
   MAX  = k
   modd = i - 1
 END DO
 
!     COMPUTE CUMULATIVE DISTRIBUTION, AND MEDIAN.
 
 DO  i = 2,mm1
   hist(i) = hist(i) + hist(i-1)
 END DO
 nn2 = nn/2
 DO  i = 1,mm1
   IF (hist(i) > nn2) EXIT
 END DO
 60 median = i - 1
 RETURN
END SUBROUTINE dist
