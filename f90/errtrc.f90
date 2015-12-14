SUBROUTINE errtrc ( NAME, ival )
     
 CHARACTER (LEN=*), INTENT(OUT)           :: NAME
 INTEGER, INTENT(IN)                      :: ival
 
 COMMON / system / isysbf, nout
 WRITE ( nout, * ) ' ERRTRC CALLED'
 
 WRITE ( nout, * ) ' NAME=',NAME
 WRITE ( nout, * ) ' IVAL=',ival
 RETURN
END SUBROUTINE errtrc
