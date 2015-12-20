SUBROUTINE fbsrdm ( mcb   , icore , rcore , dcore,  &
        memtot, buff, lasind, ipos )
     
!  FBSRDM - This routine will store an entire matrix in memory
!           if sufficient memory exists.  The matrix
!           is stored in memory according to the following scheme:
!           (Subroutine FERRDM is very similiar to this subroutine)
 
!  1st word = current column number
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
!  3+(ntms*prec)+1 = row position of first element in above string
!  3+(ntms*prec)+2 = number of terms in ABOVE string (ntms)
 
!  The above data repeats for all strings within a column and then
!  for all columns in the matrix.
 
!  Argument list :
!     MCB    - Matrix control block for input matrix
!     ICORE  - Memory for storage of data (integer)
!     RCORE  - Same location as ICORE but real single reference
!     DCORE  - Same location as ICORE but real double reference
!     MEMTOT - Total amount of memory available for this data
!     BUFF   - Buffer allocation for input matrix
!     LASIND - Memory index of last string stored in memory
!     IPOS   - 6 word array with the following information
!              (1) = last column read into memory
!              (2) = block number of following column not read into memory
!              (3) = current logical record pointer for following column
!                    not read into memory
!              (4) = current buffer pointer for following record not read
!                    into memory
!              (5) = last block number in file
!              (6) = current logical record pointer for last record in file
!              (7) = current buffer pointer for last record in file
 
 
 INTEGER, INTENT(IN)                      :: mcb(7)
 INTEGER, INTENT(OUT)                     :: icore(1)
 REAL, INTENT(OUT)                        :: rcore(1)
 DOUBLE PRECISION, INTENT(OUT)            :: dcore(1)
 INTEGER, INTENT(IN)                      :: memtot
 INTEGER, INTENT(IN OUT)                  :: buff(2)
 INTEGER, INTENT(OUT)                     :: lasind
 INTEGER, INTENT(OUT)                     :: ipos(7)
 DOUBLE PRECISION :: dxl
 REAL :: rxl(1)
 INTEGER :: rd, rdrew, wrt, wrtrew, rew
 
 INTEGER :: iblk(20)
 COMMON  /zzzzzz/ dxl(1)
 COMMON  /system/ ksystm(65)
 COMMON  /names / rd, rdrew, wrt, wrtrew, rew
 EQUIVALENCE      ( ksystm( 2), nout  )
 EQUIVALENCE      ( dxl,rxl )
 
 mem          = 1
 ncol         = mcb( 2 )
 ntype        = mcb( 5 )
 incr         = 1
 IF ( ntype == 2 .OR. ntype == 3 ) incr = 2
 IF ( ntype == 4 ) incr = 4
 ntwds        = 0
 ipos( 1 )    = ncol
 DO   i  = 2,7
   ipos( i ) = 0
 END DO
 DO  i  = 1,20
   iblk(i)      = 0
 END DO
 iblk(1)      = mcb( 1 )
 iblk(9)      = 1
 iblk(10)     = 1
 CALL gopen  ( mcb, buff, rdrew )
 CALL REWIND ( mcb)
 CALL skprec ( mcb, 1 )
 DO  jcol = 1,ncol
   iblk(8)      = -1
   lasind       = mem - 1
   CALL dscpos  ( mcb, iblock, iclr, icbp )
   100 CALL getstr(*1000,iblk(1))
   INDEX        = iblk( 5 )
   ntms         = iblk( 6 )
   jrow         = iblk( 4 )
   ntwds        = ntwds + 4 + ntms*incr
   IF ( ntwds > memtot ) GO TO 2000
   icore(mem)   = jcol
   icore(mem+1) = ntms
   SELECT CASE ( ntype )
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
   mindex     = mem + 1
   DO  ii  = 1,ntms
     rcore(mindex+ii) = rxl(INDEX+ii-1)
   END DO
   mem        = mem + 2 + ntms
   GO TO 180
   120 CONTINUE
   mindex     = mem/2+1
   DO  ii  = 1,ntms
     dcore(mindex+ii) = dxl(INDEX+ii-1)
   END DO
   mem        = mem + 2 + ntms*2
   GO TO 180
   130 CONTINUE
   mindex     = mem + 1
   ntms2      = ntms*2
   DO  ii  = 1,ntms2
     rcore(mindex+ii) = rxl(INDEX+ii-1)
   END DO
   mem        = mem + 2 + ntms2
   GO TO 180
   140 CONTINUE
   mindex     = mem/2+1
   ntms2      = ntms*2
   DO  ii = 1,ntms2
     dcore(mindex+ii) = dxl(INDEX+ii-1)
   END DO
   mem        = mem + 2 + ntms*4
   GO TO 180
   180 CONTINUE
   icore(mem  ) = jrow
   icore(mem+1) = ntms
   mem          = mem + 2
   185 CALL endget (iblk( 1 ) )
   GO TO 100
   1000 CONTINUE  
 END DO
 lasind    = mem - 1
 GO TO 7000
 2000 ipos( 1 ) = jcol - 1
 ipos( 2 ) = iblock
 ipos( 3 ) = iclr
 ipos( 4 ) = icbp
 CALL skprec ( mcb, ncol-jcol+1 )
 CALL dscpos ( mcb, iblock, iclr, icbp )
 ipos( 5 ) = iblock
 ipos( 6 ) = iclr
 ipos( 7 ) = icbp
 7000 CONTINUE
 CALL CLOSE ( mcb , rew )
 
 RETURN
END SUBROUTINE fbsrdm
