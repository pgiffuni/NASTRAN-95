SUBROUTINE dsocff ( dsname, iunit, istatus )
     
 CHARACTER (LEN=80), INTENT(OUT)          :: dsname
 INTEGER, INTENT(OUT)                     :: iunit
 INTEGER, INTENT(OUT)                     :: istatus
 
 COMMON / system / sysbuf, iwr
 COMMON / machin / mac(3), lqro
 INCLUDE  'DSIOF.COM'
!  OPEN AND CLOSE FILE IN ORDER TO DELETE SPACE
 OPEN  ( UNIT=iunit, FILE=dsname    , IOSTAT=istatus, ERR=100  &
     ,       STATUS='UNKNOWN' )
 100   CLOSE ( UNIT=iunit, STATUS='DELETE', IOSTAT=istatus, ERR=701 )
! NOW, OPEN FILE AS NEW FOR NASTRAN
!      print *,' dsocff,nbuff=',nbuff
 nbuff4 = nbuff * ( MOD(lqro,100) / 10 )
 OPEN  ( UNIT=iunit, FILE=dsname, RECL=nbuff4, STATUS='NEW'  &
     ,       ACCESS='direct', FORM='unformatted',IOSTAT=istatus ,       ERR=702 )
 GO TO 777
 701   WRITE ( iwr, 901 ) iunit, istatus, dsname
 901   FORMAT(//,' FATAL ERROR IN DSOCFF, UNABLE TO CLOSE UNIT=',i4  &
     ,         ' STATUS='i4 ,       /,' FILE NAME=',a80 )
 iccerr = istatus
 CALL dsmsg  ( 101 )
 CALL mesage ( -61, 0, 0 )
 702   WRITE ( iwr, 902 ) iunit, istatus, dsname
 902   FORMAT(//,' FATAL ERROR IN DSOCFF, UNABLE TO OPEN UNIT=',i4  &
     ,         ' STATUS=',i4 ,       /,' FILE NAME=',a80 )
 iccerr = istatus
 CALL dsmsg  ( 101 )
 CALL mesage ( -61, 0, 0 )
 777   RETURN
END SUBROUTINE dsocff
