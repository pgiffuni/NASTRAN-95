SUBROUTINE ferrdm ( mcb,nidx,memtot,ibuffi,lasind,ipos )
     
!  FERRDM - This routine will store an entire matrix in memory
!           if sufficient memory exists.  The matrix
!           is stored in memory according to the following scheme:
 
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
!  3+(ntms*prec)      } (where prec=1 for s.p.; =2 for d.p. )
!  3+(ntms*prec)+1 = row position of first element in above string
!  3+(ntms*prec)+2 = number of terms in ABOVE string (ntms)
 
!  The above data repeats for all strings within a column and then
!  for all columns in the matrix.
 
!  Argument list :
!     MCB    - Matrix control block for input matrix
!     NIDX   - Memory index for storing matrix data
!     MEMTOT - Total amount of memory available for this data
!     IBUFFI - Buffer allocation for input matrix
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
 INTEGER, INTENT(IN)                      :: nidx
 INTEGER, INTENT(IN)                      :: memtot
 INTEGER, INTENT(IN OUT)                  :: ibuffi
 INTEGER, INTENT(OUT)                     :: lasind
 INTEGER, INTENT(OUT)                     :: ipos(7)
 DOUBLE PRECISION :: dcore(1), dxl(1)
 REAL :: rcore(1), rxl(1)
 INTEGER :: rd, rdrew, wrt, wrtrew, rew, ixl(1)
 
 DIMENSION        iblk(20)
 COMMON  /system/ ksystm(65)
 COMMON  /zzzzzz/ icore(1)
 COMMON  /names / rd, rdrew, wrt, wrtrew, rew
 EQUIVALENCE      ( ksystm( 2), nout              )
 EQUIVALENCE      ( ksystm(55), iprec             )
 EQUIVALENCE      ( icore,dcore,rcore,dxl,rxl,ixl )
 
 mem          = nidx
 ncol         = mcb( 2 )
 ntwds        = 0
 ipos( 1 )    = ncol
 DO  i  = 1,20
   iblk(i)      = 0
 END DO
 iblk(1)      = mcb( 1 )
 iblk(9)      = 1
 iblk(10)     = 1
 CALL gopen  ( mcb, icore( ibuffi ), rdrew )
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
   ntwds        = ntwds + 4 + ntms*iprec
   IF ( ntwds > memtot ) GO TO 2000
   icore(mem)   = jcol
   icore(mem+1) = ntms
   IF ( iprec == 1 ) GO TO 160
   mindex         = mem/2+1
   DO  ii = 1,ntms
     dcore(mindex+ii) = dxl(INDEX+ii-1)
   END DO
   GO TO 180
   160 mindex     = mem + 1
   DO  ii  = 1,ntms
     rcore(mindex+ii) = rxl(INDEX+ii-1)
   END DO
   180 CONTINUE
   mem          = mem + 2 + ntms*iprec
   icore(mem  ) = jrow
   icore(mem+1) = ntms
   mem          = mem + 2
   185 CALL endget (iblk( 1 ) )
   GO TO 100
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
END SUBROUTINE ferrdm
