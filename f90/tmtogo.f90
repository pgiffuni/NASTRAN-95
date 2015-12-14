SUBROUTINE tmtogo (togo)
     
!     TO COMPUTE THE TIME (IN SECONDS) REMAINING
 
 
 INTEGER, INTENT(OUT)                     :: togo
 INTEGER :: tbegin, tprob,tnow
 COMMON /system/ xsys(17),tmbegn
 COMMON /stime / tprob
 
!     GET PRESENT TIME
 
 CALL klock (tnow)
 
!     COMPUTE TIME TO GO
 
 tbegin = tmbegn
 togo = tprob - (tnow - tbegin)
 RETURN
END SUBROUTINE tmtogo
