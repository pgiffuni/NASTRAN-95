SUBROUTINE dsopen ( dsname, iunit, iop )
     
 CHARACTER (LEN=80), INTENT(IN OUT)       :: dsname
 INTEGER, INTENT(IN OUT)                  :: iunit
 INTEGER, INTENT(IN)                      :: iop
 
 INCLUDE          'DSIOF.COM'
!      print *,' dsopen,iunit,iop,dsname=',iunit,iop,dsname
 IF ( iop /= 1 ) CALL dsopff ( dsname, iunit, iccer )
 IF ( iop == 1 ) CALL dsocff ( dsname, iunit, iccer )
 numopn = numopn + 1
 700  RETURN
END SUBROUTINE dsopen
