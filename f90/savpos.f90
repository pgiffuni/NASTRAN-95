SUBROUTINE savpos ( FILE, ipos )
     INCLUDE 'DSIOF.COM'
 
 INTEGER, INTENT(IN)                      :: FILE
 INTEGER, INTENT(OUT)                     :: ipos
 COMMON / ddiosv / iflpos( 2,80 )
 
 
 NAME = FILE
 CALL dsgefl
 ipos = iflpos( 1,ifilex )*mulq2 + iflpos( 2, ifilex )
 IF (iprvop == 0) ipos = fcb(3,ifilex)*mulq2 + fcb(4,ifilex)
 RETURN
END SUBROUTINE savpos
