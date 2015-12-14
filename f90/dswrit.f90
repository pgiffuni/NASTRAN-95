SUBROUTINE dswrit ( iunit, buff, LEN, irec, iccerr )
     
 INTEGER, INTENT(OUT)                     :: iunit
 INTEGER, INTENT(IN)                      :: buff( LEN )
 INTEGER, INTENT(IN OUT)                  :: LEN
 INTEGER, INTENT(IN OUT)                  :: irec
 INTEGER, INTENT(OUT)                     :: iccerr
 
 COMMON / system / sysbuf, iwr
 INCLUDE          'DSIOF.COM'
!      print *,' dswrit,len,IREC,UNIT=',len,irec,iunit
 IF ( irec <= 0 ) GO TO 701
 WRITE ( iunit, REC=irec, IOSTAT=istat, ERR=702 ) buff
 iccerr = 0
 GO TO 777
 701   WRITE( iwr, 901 ) iunit, irec, mdsnam( iunit )
 901   FORMAT(//' ERROR IN DSWRIT, BAD RECORD NO., UNIT=',i4,' REC=',i5  &
     ,      /,' FILE NAME=',a80 )
 iccerr = istat
 CALL dsmsg  ( 101 )
 CALL mesage ( -61, 0, 0 )
 702   WRITE( iwr, 902 ) iunit, irec, istat, mdsnam( iunit )
 902   FORMAT(//', ERROR ENCOUNTERED IN DSWRCC, UNIT=',i5,' RECORD='  &
     , i5,' STATUS=',i9,/' DSNAME=',a80 )
 iccerr = istat
 CALL dsmsg  ( 101 )
 CALL mesage ( -61, 0, 0 )
 777   CONTINUE
 numwri = numwri + 1
 RETURN
END SUBROUTINE dswrit


