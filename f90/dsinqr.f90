SUBROUTINE dsinqr ( dsn, istat, isize )
!        DSINQR DETERMINES THE EXISTANCE OF A FILE:
!            DSN   ( INPUT  )   FILE NAME
!            ISTAT ( OUTPUT )   =0, IF NOT EXIST; =1, IF EXIST
!            ISIZE ( OUTPUT )   = FILE SIZE IN GINO BLOCKS
 
 
 CHARACTER (LEN=*), INTENT(OUT)           :: dsn
 INTEGER, INTENT(OUT)                     :: istat
 INTEGER, INTENT(OUT)                     :: isize
 LOGICAL :: avail
 
 
 INQUIRE( FILE=dsn, EXIST=avail, NEXTREC = nrec )
 istat = 0
 IF ( avail ) istat = 1
 isize = nrec - 1
 RETURN
END SUBROUTINE dsinqr
