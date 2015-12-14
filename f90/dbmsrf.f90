SUBROUTINE dbmsrf (  NAME, ifilex )
!********************************************************************
!  DBMNAM RETURNS THE UNIT NUMBER ASSOCIATED WITH A DMAP FILE NAME
!     ARGUMENTS
!       NAME    (INPUT) (2A4) DMAP FILE NAME
!       IFILEX  (OUTPUT) (I)  UNIT ASSOCIATED WITH FILE NAME IN FIAT
!********************************************************************
 
 INTEGER, INTENT(IN)                      :: NAME(2)
 INTEGER, INTENT(OUT)                     :: ifilex
 INTEGER :: fist(100), fiat(100)
 COMMON / xfist / fist
 COMMON / xfiat / fiat
 
 lim  = fiat(2)*11 + 3
 DO  i = 4, lim, 11
   IF ( NAME( 1 ) /= fiat( i+1 ) .OR. NAME( 2 ) /= fiat( i+2 ) ) CYCLE
   ifilex = IAND( fiat( i ), 32767 )
   GO TO 700
 END DO
 ifilex = 0
 700   RETURN
END SUBROUTINE dbmsrf
