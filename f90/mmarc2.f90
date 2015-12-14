SUBROUTINE mmarc2 ( zi, zd )
!  MMARC2 - This routine will store a matrix column in memory in compact
!           form and in real double precision.  The input matrix is
!           assumed to be stored as either real single or double precision.
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
!     ZD  - Same location as ZI but real double reference
 
 
 INTEGER, INTENT(OUT)                     :: zi(1)
 DOUBLE PRECISION, INTENT(OUT)            :: zd(1)
 
 INTEGER :: iblk(15)
 REAL :: rxl(1)
 DOUBLE PRECISION :: dxl
 INCLUDE          'MMACOM.COM'
 COMMON  /zzzzzz/ dxl(1)
 EQUIVALENCE      ( dxl,rxl )
 
 mem      = 1
 DO  i  = 1,15
   iblk(i)  = 0
 END DO
 iblk(1)  = irfile
 iblk(8)  = -1
 lasind   = mem - 1
 zi( mem ) = 0
 100 CALL getstr ( *7000, iblk )
 itype    = iblk( 2 )
 jrow     = iblk( 4 )
 INDEX    = iblk( 5 )
 ntms     = iblk( 6 )
 zi(mem)  = jrow
 zi(mem+1)= ntms
 SELECT CASE ( itype )
   CASE (    1)
     GO TO  110
   CASE (    2)
     GO TO  120
 END SELECT
 110 CONTINUE
 mindex   = mem/2 + 1
 DO  ii = 1,ntms
   zd( mindex+ii ) = SIGN*rxl( INDEX+ii-1 )
 END DO
 GO TO 180
 120 CONTINUE
 mindex     = mem/2+1
 DO  ii  = 1,ntms
   zd( mindex+ii ) = SIGN*dxl( INDEX+ii-1 )
 END DO
 180 CONTINUE
 mem        = mem + 2 + ntms*2
 CALL endget ( iblk )
 GO TO 100
 7000 CONTINUE
 lasind     = mem - 1
 RETURN
END SUBROUTINE mmarc2
