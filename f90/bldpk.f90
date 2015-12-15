SUBROUTINE bldpk ( itypin, itypot, FILE, BLOCK, iflag )
     INCLUDE 'PAKBLK.COM'
 INCLUDE 'DSIOF.COM'
 
 
 INTEGER, INTENT(IN)                      :: itypin
 INTEGER, INTENT(IN)                      :: itypot
 INTEGER, INTENT(IN)                      :: FILE
 INTEGER, INTENT(IN OUT)                  :: BLOCK(15)
 INTEGER, INTENT(IN)                      :: iflag
 
 itrail = iflag
 itypi  = itypin
 itypo  = itypot
 NAME   = FILE
 IF ( itypi < 1 .OR. itypi > 4 ) GO TO 40
 IF ( itypo < 1 .OR. itypo > 4 ) GO TO 40
 IF ( iflag == 0 ) GO TO 20
 CALL dsblpk ( BLOCK )
 GO TO 30
 20      itrail = 0
 CALL dsblpk ( iblka )
 30      GO TO 700
 40      IF ( iflag == 0 ) CALL dsmsg1 ( iblka )
 IF ( iflag /= 0 ) CALL dsmsg1 ( BLOCK )
 CALL dsmsg( 118 )
 700     RETURN
END SUBROUTINE bldpk
