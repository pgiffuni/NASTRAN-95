SUBROUTINE dbmnam ( igname, NAME, ifilex )
!********************************************************************
!  DBMNAM RETURNS THE DMAP NAME AND TRAILER FOR A GIVEN GINO FILE
!     ARGUMENTS
!       IGNAME  (INPUT )   GINO FILE NAME (E.G., 101,201,301)
!       FIST    (INPUT )   COMMON BLOCK /XFIST/
!       FIAT    (INPUT )   COMMON BLOCK /XFIAT/
!       NAME    (OUTPUT)   (2A4) DMAP FILE NAME
!********************************************************************
 
 INTEGER, INTENT(IN)                      :: igname
 INTEGER, INTENT(OUT)                     :: NAME(2)
 INTEGER, INTENT(IN)                      :: ifilex
 INTEGER :: fist(100), fiat(100)
 INTEGER :: BLANK, pool
 COMMON / xfist / fist
 COMMON / xfiat / fiat
 DATA     pool / 4HPOOL /, BLANK / 4H    /
 
 IF ( igname <= 100 .OR. igname >= 400 ) GO TO 100
 lim  = fist(2)*2 - 1
 DO  i = 1, lim, 2
   IF ( igname /= fist(2+i) ) CYCLE
   INDEX     = fist(3+i)
   NAME(1)   = fiat( INDEX+2 )
   NAME(2)   = fiat( INDEX+3 )
   GO TO 700
 END DO
 NAME(1)   = 0
 NAME(2)   = 0
 GO TO 700
 100   NAME(1) = igname
 NAME(2) = BLANK
 IF ( ifilex == 22 ) NAME(1) = pool
 700   RETURN
END SUBROUTINE dbmnam
