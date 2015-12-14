SUBROUTINE mmarc1 ( zi, zr )
!  MMARC1 - This routine will store a matrix column in memory in compact
!           form and in real single precision.  The input matrix is
!           assumed to be stored as real single precision.
!           The column is stored in memory according to the following scheme:
 
 
!  1st word = row position of first element in following string
!  2nd word = number of terms in string (ntms)
!  3rd word           }
!     |               }
!     |               } = actual
!     |               }   matrix
!     |               }   string
!     |               }   data
!     |               }
!     |               }
!  3+(ntms*prec)      } (where prec=1 for s.p.;  =2 for d.p. )
 
!  The above data repeats for all strings within a column
 
!  Argument list :
!     ZI  - Memory for storage of data (integer)
!     ZR  - Same location as ZI but real single reference
 
 
 INTEGER, INTENT(OUT)                     :: zi(1)
 REAL, INTENT(OUT)                        :: zr(1)
 
 
 INTEGER :: iblk(15)
 INCLUDE          'MMACOM.COM'
 COMMON  /zzzzzz/ rxl(1)
 
 mem       = 1
 DO  i   = 1,15
   iblk(i)   = 0
 END DO
 iblk(1)   = irfile
 iblk(8)   = -1
 lasind    = mem - 1
 zi( mem ) = 0
 100   CALL getstr ( *7000, iblk )
 jrow      = iblk( 4 )
 INDEX     = iblk( 5 )
 ntms      = iblk( 6 )
 zi(mem)   = jrow
 zi(mem+1) = ntms
 mem       = mem + 1
 DO  ii = 1,ntms
   zr(mem+ii)= SIGN*rxl(INDEX+ii-1)
 END DO
 mem       = mem + 1 + ntms
 CALL endget ( iblk )
 GO TO 100
 7000 CONTINUE
 lasind    = mem - 1
 RETURN
END SUBROUTINE mmarc1
