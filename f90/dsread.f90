SUBROUTINE dsread ( iunit, buff, LEN, irec )
     
 INTEGER, INTENT(OUT)                     :: iunit
 INTEGER, INTENT(IN OUT)                  :: buff( LEN )
 INTEGER, INTENT(IN OUT)                  :: LEN
 INTEGER, INTENT(IN OUT)                  :: irec
 
 COMMON / system / sysbuf, iwr
 INCLUDE          'DSIOF.COM'
 
 IF ( irec  < 0 ) GO TO 701
!      PRINT *,' DSREAD,LEN,IREC,IUNIT=',LEN,IREC,IUNIT
 istat=0
 READ ( iunit, REC=irec, ERR=702, IOSTAT=istat ) buff
 IF ( istat == 0 ) GO TO 777
 ioerr = istat
 CALL dsmsg  ( 101 )
 CALL mesage ( -61, 0, 0 )
 701   WRITE ( iwr, 901 ) iunit, irec, mdsnam( iunit )
 901   FORMAT(//' ERROR IN DSREAD-BAD REC NO., UNIT=',i4,' REC=',i4  &
     ,      /,' FILE NAME=',a72)
 iccerr = 0
 CALL dsmsg  ( 101 )
 CALL mesage ( -61, 0, 0 )
 GO TO 777
 702   WRITE( iwr, 902 ) iunit, irec, istat, mdsnam( iunit )
 902   FORMAT(//', ERROR ENCOUNTERED IN DSREAD, UNIT=',i5,' RECORD='  &
     , i5,' STATUS=',i9,/' DSNAME=',a72 )
 iccerr = istat
 CALL dsmsg( 101 )
 CALL mesage( -61, 0, 0 )
 GO TO 777
 777   CONTINUE
 numrea = numrea + 1
 RETURN
END SUBROUTINE dsread
