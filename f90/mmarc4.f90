SUBROUTINE mmarc4 ( zi, zd )
!  MMARC4 - This routine will store a matrix column in memory in compact
!           form and in complex double precision.  The input matrix can
!           be stored in any precision or type.
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
 DOUBLE PRECISION :: dxl(1)
 INCLUDE          'MMACOM.COM'
 COMMON  /zzzzzz/ rxl(1)
 EQUIVALENCE      ( rxl, dxl )
 
 mem      = 1
 DO  i  = 1,15
   iblk(i)  = 0
 END DO
 iblk(1)  = irfile
 iblk(8)  = -1
 lasind   = mem - 1
 zi( mem) = 0
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
   CASE (    4)
     GO TO  140
 END SELECT
 110 CONTINUE
 mindex   = mem/2 + 1
 DO  ii = 1,ntms
   zd( mindex+1 ) = SIGN*rxl( INDEX+ii-1 )
   zd( mindex+2 ) = 0.d0
   mindex = mindex+2
 END DO
 GO TO 180
 120 CONTINUE
 mindex = mem/2 + 1
 DO  ii = 1,ntms
   zd( mindex+1 ) = SIGN*dxl( INDEX+ii-1 )
   zd( mindex+2 ) = 0.d0
   mindex  = mindex + 2
 END DO
 GO TO 180
 130 CONTINUE
 mindex  = mem/2 + 1
 ntms2   = ntms*2
 DO  ii = 1,ntms2
   zd( mindex+ii ) = SIGN*rxl( INDEX+ii-1 )
 END DO
 GO TO 180
 140 CONTINUE
 mindex  = mem/2 + 1
 ntms2   = ntms*2
 DO  ii = 1,ntms2
   zd( mindex+ii ) = SIGN*dxl( INDEX+ii-1 )
 END DO
 180 CONTINUE
 mem     = mem + 2 + ntms*4
 CALL endget ( iblk )
 GO TO 100
 1000 CONTINUE
 lasind  = mem - 1
 RETURN
END SUBROUTINE mmarc4
