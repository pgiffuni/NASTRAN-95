SUBROUTINE waltim (walsec)
     
!     THIS ROUTINE OBTAINS THE CURRENT WALL CLOCK TIME IN SECONDS,
!     PASS MID-NIGHT
 
 
 
 INTEGER, INTENT(OUT)                     :: walsec
 INTEGER :: time(3)
 
 CALL itime ( time )
 walsec = time(1) * 3600 + time(2) * 60 + time(3)
 RETURN
END SUBROUTINE waltim
