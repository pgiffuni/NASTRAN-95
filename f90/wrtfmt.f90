SUBROUTINE wrtfmt ( iout, nwds, ifmt )
     
 
 INTEGER, INTENT(IN OUT)                  :: iout
 INTEGER, INTENT(IN OUT)                  :: nwds
 CHARACTER (LEN=1), INTENT(IN OUT)        :: ifmt(*)
 
 CALL forwrt ( ifmt, iout, nwds )
 RETURN
END SUBROUTINE wrtfmt
