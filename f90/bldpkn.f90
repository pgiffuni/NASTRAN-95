SUBROUTINE bldpkn ( FILE, BLOCK,  mcb )
     INCLUDE 'PAKBLK.COM'
 INCLUDE 'DSIOF.COM'
 
 
 INTEGER, INTENT(IN)                      :: FILE
 INTEGER, INTENT(IN)                      :: BLOCK( 15 )
 INTEGER, INTENT(IN OUT)                  :: mcb( 7 )
 
 NAME = FILE
 IF ( BLOCK( 1) == 0 ) GO TO 10
 CALL dsbpnk ( BLOCK, mcb )
 GO TO 20
 10    CALL dsbpnk ( iblka, mcb )
 20    CONTINUE
 RETURN
END SUBROUTINE bldpkn
