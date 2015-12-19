SUBROUTINE fwdrec ( *, FILE )

  INCLUDE 'DSIOF.COM'
 
  INTEGER, INTENT(IN)                      :: FILE
 
 NAME = FILE
 CALL dsgefl
 CALL dsfwr1
 CALL dssdcb
 IF ( iretrn == 1 ) RETURN 1
 
 RETURN
END SUBROUTINE fwdrec
