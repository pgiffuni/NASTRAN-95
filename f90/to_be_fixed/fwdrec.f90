SUBROUTINE fwdrec ( *, FILE )
     INCLUDE 'DSIOF.COM'
 
 
 , INTENT(IN OUT)                         :: *
 INTEGER, INTENT(IN)                      :: FILE
 
 NAME = FILE
 CALL dsgefl
 CALL dsfwr1
 CALL dssdcb
 IF ( iretrn == 1 ) RETURN 1
 RETURN
END SUBROUTINE fwdrec
