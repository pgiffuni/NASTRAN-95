SUBROUTINE strap2 (ti)
     
 
!*****
! THIS ROUTINE IS PHASE II OF STRESS DATA RECOVERY FOR THE TRAPEZOIDAL
! CROSS SECTION RING
!*****
 
 
 
 
 REAL, INTENT(IN)                         :: ti(4)
 
 DIMENSION          dum3(225)
 DIMENSION          stres(100),    force(25)
 DIMENSION          istres(100),   iforce(25)
 
 
! SDR2 VARIABLE CORE
 
 COMMON   /zzzzzz/  zz(1)
 
 
! SDR2 BLOCK FOR POINTERS AND LOADING TEMPERATURES
 
 COMMON   /sdr2x4/ dum1(33)  &
     ,                  icstm,    ncstm,    ivec,     ivecn  &
     ,                  templd,   eldefm
 
 
! SDR2 INPUT AND OUTPUT BLOCK
 
 COMMON   /sdr2x7/ idel,     igp(4),   tz  &
     ,                  sel(240), ts(4),    ak(144)
 
 
! SCRATCH BLOCK
 
 COMMON   /sdr2x8/ disp(12), eforc(12),estres(20)
 
 
 EQUIVALENCE (dum3(1) , idel)
 EQUIVALENCE  (dum3(101) , stres(1) , istres(1))
 EQUIVALENCE  (dum3(201) , force(1) , iforce(1))
 EQUIVALENCE (ldtemp, templd)
 
 
! INITIALIZE COUNTERS
 
 ndof  = 3
 numpt = 4
 n = ndof * numpt
 nsp   = 5
 ncomp =  4
 ns = nsp * ncomp
 
 
! LOCATE THE DISPLACEMENTS
 
 k = 0
 DO  i = 1,numpt
   iloc = ivec + igp(i) - 2
   DO  j = 1,ndof
     iloc = iloc + 1
     k = k + 1
     disp(k) = zz(iloc)
   END DO
 END DO
 
 
! COMPUTE THE GRID POINT FORCES
 
 CALL gmmats ( ak(1) , n, n, 0, disp(1) , n, 1, 0, eforc(1) )
 
 
! COMPUTE THE STRESSES
 
 CALL gmmats ( sel(1), ns, n, 0, disp(1) , n, 1, 0, estres(1) )
 
 
! COMPUTE THERMAL STRESS IF THERMAL LOAD EXISTS
! AND SUBTRACT FROM APPARENT STRESS
 
 IF (ldtemp == (-1)) GO TO 300
 
 k = 0
 DO  i = 1,nsp
   dt = ti(i) - tz
   IF (i == 5) dt = (ti(1)+ti(2)+ti(3)+ti(4)) / 4.0E0  -  tz
   DO  j = 1,ncomp
     k = k + 1
     estres(k) = estres(k) - dt * ts(j)
   END DO
 END DO
 
 300 CONTINUE
 
 
! STORE RESULTS FOR OUTPUT PRINT
 
 k = 0
 j = 1
 istres(1)   = idel
 DO  kk = 1,nsp
   DO  i = 1,ncomp
     j = j + 1
     k = k + 1
     stres(j) = estres(k)
   END DO
 END DO
 
 
 k = 0
 j = 1
 iforce(1)   = idel
 DO  i = 1,numpt
   DO  kk= 1,ndof
     j = j + 1
     k = k + 1
     force(j) = eforc(k)
   END DO
 END DO
 
 RETURN
END SUBROUTINE strap2
