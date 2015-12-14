SUBROUTINE dscpos ( FILE, icblk, iclr, icbp )
     
! RETURNS THE CURRENT BLOCK NUMBER "ICBLK", CURRENT LOGICAL RECORD
! POINTER "ICLR" AND CURRENT BUFFER POINT "ICBP" FOR "FILE"
 
 INCLUDE 'DSIOF.COM'
 
 
 INTEGER, INTENT(IN)                      :: FILE
 INTEGER, INTENT(OUT)                     :: icblk
 INTEGER, INTENT(OUT)                     :: iclr
 INTEGER, INTENT(OUT)                     :: icbp
 
 NAME   = FILE
 CALL dsgefl
 icblk  = fcb( 4, ifilex )
 iclr   = indclr - indbas + 1
 icbp   = indcbp - indbas + 1
 RETURN
END SUBROUTINE dscpos
