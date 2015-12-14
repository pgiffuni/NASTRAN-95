SUBROUTINE bckrec ( FILE )
     INCLUDE 'DSIOF.COM'
 
 
 INTEGER, INTENT(IN)                      :: FILE
 
 NAME = FILE
 CALL dsgefl
 CALL dsbrc1
 CALL dssdcb
 RETURN
END SUBROUTINE bckrec
