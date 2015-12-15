SUBROUTINE sofio ( isop, iblknm, buf )
     INCLUDE 'GINOX.COM'
 INCLUDE 'DSIOF.COM'
 
 INTEGER, INTENT(IN)                      :: isop
 INTEGER, INTENT(IN)                      :: iblknm
 INTEGER, INTENT(IN OUT)                  :: buf(10)
 COMMON / sofcom / nfiles, filnam( 10 ), filsiz( 10 )
 COMMON / sys    / blksiz, dirsiz, supsiz, avblks, hiblk
 COMMON / system / isysbf, iwr
 COMMON / sofdsn / sofdsn(10)
 CHARACTER (LEN=4) :: filnam
 CHARACTER (LEN=80) :: dsname
 CHARACTER (LEN=80) :: sofdsn
 INTEGER :: filsiz, hiblk
 
 IF ( lensof( 1 ) /= 0 ) GO TO 20
 numblk = 1
 IF ( lenwpb /= 0 ) numblk = isysbf / lenwpb
 DO  k = 1, nfiles
   lensof( k ) = 0
   dsname      = sofdsn(k)
   CALL dsinqr ( dsname, istat, isize)
   IF (istat == 0) CYCLE
   lensof( k ) = filsiz( k )
 END DO
 lasfil = 0
 20      CONTINUE
 IF ( isop == 7 ) GO TO 200
 iblk   = iblknm
 IF ( iblk <= 0 ) GO TO 700
 ifile  = 0
 DO  k = 1, nfiles
   IF ( iblk <= filsiz( k ) )GO TO 30
   iblk   = iblk - filsiz( k )
   CYCLE
   30      ifile  = k
   GO TO 100
 END DO
 WRITE( iwr, 9910 ) iblknm
 9910    FORMAT(' *** SUBSTRUCTURING ERROR - BLOCK NUMBER OUT OF RANGE ',  &
     ' -  BLOCK NUMBER IS ',i10)
 WRITE( iwr, 9915 )
 9915    FORMAT( //,' THE FOLLOWING SOF FILES WERE AVAILABLE',//)
 DO  k = 1, nfiles
   WRITE( iwr, 9920 ) sofdsn( k ), filsiz( k ), lensof( k )
   9920    FORMAT(' FILE ',a72' HAS ',i10, ' BLOCKS - BLOCKS USED ',i10)
 END DO
 CALL mesage (-61, 0, 0)
 100     IF ( lasfil == ifile ) GO TO 120
 IF ( lasfil /= 0 ) CALL dsclos ( 90 )
 ialloc = numblk * filsiz(ifile)
 dsname = sofdsn( ifile )
 iop = 0
 CALL dsopen (  dsname, 90, iop )
 lasfil = ifile
 120     IF ( isop == 1 ) GO TO 140
 IF ( ( iblk - lensof( ifile ) ) <= 1 ) GO TO 130
 num  = iblk - lensof( ifile ) - 1
 IF ( num == 0 ) GO TO 130
 DO  k = 1, num
   lensof( ifile ) = lensof( ifile ) + 1
   CALL dswrit ( 90, buf(4), nbuff, lensof( ifile ),icerr )
   IF ( icerr /= 0 ) GO TO 701
 END DO
 130     CONTINUE
 CALL dswrit ( 90, buf(4), nbuff, iblk, icerr )
 IF ( icerr /= 0 ) GO TO 701
 IF ( iblk   > lensof( ifile ) ) lensof( ifile ) = iblk
 IF ( iblknm > hiblk ) hiblk = iblknm
 GO TO 700
 140     CALL dsread ( 90, buf(4), nbuff, iblk )
 GO TO 700
 200     CALL dsclos( 90 )
 lasfil = 0
 700     GO TO 7000
 701     IF ( icerr == 28 ) WRITE ( iwr, 901 )
 IF ( icerr /= 28 ) WRITE ( iwr, 902 )
 CALL mesage (-61, 0, 0)
 901     FORMAT(///,' INSUFFICIENT SPACE FOR SOF FILE ON DEFAULT',  &
     ' DEVICE---JOB ABORTED.')
 902     FORMAT(///,' I/O ERROR OCCURRED ON SOF FILE, JOB ABORTED')
 7000    RETURN
END SUBROUTINE sofio
