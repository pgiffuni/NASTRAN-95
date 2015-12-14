SUBROUTINE nastim (ihr, imn, isc, cpusec)
     
 INTEGER, INTENT(OUT)                     :: ihr
 INTEGER, INTENT(OUT)                     :: imn
 INTEGER, INTENT(OUT)                     :: isc
 REAL, INTENT(OUT)                        :: cpusec
 REAL :: array(2)
 REAL :: result
 
 CALL etime(array,result)
 secs   = array(2)
 ihr    = secs / 3600.
 imn    = ( secs - 3600.*ihr ) / 60.
 isc    = secs - ( 3600.*ihr ) - ( 60.*imn )
 cpusec = secs
 RETURN
END SUBROUTINE nastim
