SUBROUTINE rectyp ( FILE, itype )
     INCLUDE 'DSIOF.COM'
 INCLUDE 'XNSTRN.COM'
 
 
 INTEGER, INTENT(IN)                      :: FILE
 INTEGER, INTENT(OUT)                     :: itype
 
 NAME = FILE
 CALL dsgefl
 5      id   = IAND( ibase( indclr ), maskq1 )
 IF ( id == idssb ) GO TO 10
 IF ( id == idseb ) GO TO 20
 itype = 0
 GO TO 7000
 10      itype = 1
 GO TO 7000
 20      CALL dsrdnb
 CALL dssdcb
 GO TO 5
 7000    RETURN
END SUBROUTINE rectyp
