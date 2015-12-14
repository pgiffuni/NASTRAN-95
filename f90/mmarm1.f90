SUBROUTINE mmarm1 ( zi, zr, mempcol )
     
!  MMARM1 - This routine will store matrix columns in memory in compact
!           form and in real single precision.  The input matrix is
!           assumed to be stored as real single precision.
!           The column is stored in memory according to the following scheme:
 
!  MEMPCOL  = Input, extra memory needed for each column that is stored
!             in memory in compact form.  This is needed for methods 40
!             and 41 where for each column of "B" stored in compact form
!             in memory, there needs to be space available for a column
!             of the "D" matrix.
 
!  1st word = column number (negative)
!  2nd word = index to next column within this array
!  3st word = row position of first element in following string
!  4nd word = number of terms in string (ntms)
!  5rd word           }
!     |               }
!     |               } = actual
!     |               }   matrix
!     |               }   string
!     |               }   data
!     |               }
!     |               }
!  5+(ntms*prec)      } (where prec=1 for s.p.;  =2 for d.p. )
!     n               } Last value of last string for this column
 
!  Words 3 through 5+(ntms*prec) above data repeat for all strings
!  within a column.  Words 1 through n repeat for all columns that are
!  read into memory.
 
!  Argument list :
!     ZI  - Memory for storage of data (integer)
!     ZR  - Same location as ZI but real single reference
 
 
 INTEGER, INTENT(OUT)                     :: zi(1)
 REAL, INTENT(OUT)                        :: zr(1)
 INTEGER, INTENT(IN)                      :: mempcol
 
 
 INTEGER :: iblk(15), module( 2 )
 INCLUDE          'MMACOM.COM'
 COMMON  /system/ ibfsiz, iwr
 COMMON  /zzzzzz/ rxl(1)
 DATA             module / 4HMMAR, 4HM1   /
 
 mem       = 1
 DO  i   = 1,15
   iblk(i)   = 0
 END DO
 iblk(1)   = irfile
 
! IRCOL1, FIRST COLUMN EXPECTED FOR THIS PASS
! IRCOLN, ON INPUT, THIS IS THE LAST COLUMN THAT IS NEEDED
!         ON OUTPUT, THIS IS THE LAST COLUMN READ
! LASMEM, LAST AVAILABLE MEMORY INDEX TO THE "ZI" ARRAY
 
 icol      = ircol1
 100   CONTINUE
 iblk(8)   = -1
 lasindm   = mem - 1
 CALL dscpos ( irfile, icblk, iclr, icbp )
 CALL getstr ( *900, iblk )
!      IF ( ICOL .NE. IBLK( 12 ) ) GO TO 7001
 zi(mem  ) = -icol
 mem1      = mem + 1
 mem       = mem + 2
 105   CONTINUE
 ntms      = iblk( 6 )
 IF ( ( mem + 2 + ntms ) > lasmem ) GO TO 2000
 jrow      = iblk( 4 )
 INDEX     = iblk( 5 )
 zi(mem)   = jrow
 zi(mem+1) = ntms
 mem       = mem + 1
 DO  ii = 1,ntms
   zr(mem+ii)= SIGN*rxl(INDEX+ii-1)
 END DO
 mem       = mem + 1 + ntms
 CALL endget ( iblk )
 CALL getstr ( *1000, iblk )
 GO TO 105
 900   CONTINUE
 zi( mem )   = -icol
 mem1        = mem + 1
 mem         = mem + 2
 1000  CONTINUE
 
! CHECK If SPACE AVAILABLE FOR A FULL COLUMN OF "D" MATRIX, IF NECESSARY
 
 IF ( mem > ( lasmem-mempcol ) ) GO TO 2000
 lasmem      = lasmem - mempcol
 zi( mem1 )  = mem
 icol        = icol + 1
 IF ( icol > ircoln ) GO TO 7000
 GO TO 100
 2000  lasindm   = mem1 - 2
 
! SAVE I/O LOCATION OF LAST COLUMN FOR NEXT PASS
 
 irpos( 1 ) = icblk
 irpos( 2 ) = iclr
 irpos( 3 ) = icbp
 ircoln     = icol - 1
 IF ( ircoln < ircol1 ) CALL mesage ( -8, mem+mempcol, module )
 GO TO 7777
 7000  CONTINUE
 lasindm    = mem - 1
!      GO TO 7777
!7001  WRITE( IWR, 9001 ) ICOL, IBLK(12), IRFILE
!9001  FORMAT(' ERROR OCCURRED IN MMARM1, EXPECTED COLUMN =',I10
!     &,/,    ' BUT READ COLUMN =',I10,' FROM FILE =',I5 )
!      CALL DSMSG( 777 )
!      CALL MESAGE ( -61, 0, 0 )
 7777  RETURN
END SUBROUTINE mmarm1
