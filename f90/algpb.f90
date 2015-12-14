SUBROUTINE algpb (idat,ntype)
     
 
 INTEGER, INTENT(IN OUT)                  :: idat
 INTEGER, INTENT(OUT)                     :: ntype
 INTEGER :: na(4)
 DATA    na / 2 ,     2  ,   3 ,    1 /
!                ZERO, INTEGER, REAL, ALPHA
 
!     RETURN FROM NUMTYP IS            SET NTYPE TO
!       0 -  ZERO                        1 - ALPHA
!       1 -  INTEGER                     2 - INTEGER
!       2 -  REAL                        3 - REAL
!       3 -  BCD
 
!     BLANK IS ALPHA,  ZERO IS INTEGER UNLESS NUMTYP SET IT TO REAL
 
 itype = numtyp(idat) + 1
 ntype = na(itype)
 RETURN
END SUBROUTINE algpb
