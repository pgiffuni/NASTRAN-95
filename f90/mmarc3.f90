SUBROUTINE mmarc3 ( zi, zr )
     
!  MARRC3 - This routine will store a matrix column in memory in compact
!           form and in complex single precision.  The input matrix is
!           assumed to be stored as either real or complex single precision.
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
 COMMON  /system/ ibfsiz, iwr
 
 mem      = 1
 DO  i  = 1,15
   iblk(i)  = 0
 END DO
 iblk(1)  = irfile
 iblk(8)  = -1
 lasind   = mem - 1
 zi( mem )= 0
 100 CALL getstr ( *1000, iblk )
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
   CASE (    3)
     GO TO  130
 END SELECT
 110 CONTINUE
 mindex   = mem + 2
 DO  ii = 1,ntms
   zr( mindex   ) = SIGN*rxl( INDEX+ii-1 )
   zr( mindex+1 ) = 0.
   mindex = mindex+2
 END DO
 GO TO 180
 
! THE FOLLOWING LINE SHOULD NEVER BE REFERENCED
 
 120 CONTINUE
 WRITE( iwr, * ) ' ERROR IN MMARC3'
 STOP
 130 CONTINUE
 mindex   = mem + 1
 ntms2    = ntms*2
 DO  ii = 1,ntms2
   zr( mindex+ii ) = SIGN*rxl( INDEX+ii-1 )
 END DO
 180 CONTINUE
 mem        = mem + 2 + ntms*2
 CALL endget ( iblk )
 GO TO 100
 1000 CONTINUE
 lasind     = mem - 1
 RETURN
END SUBROUTINE mmarc3
