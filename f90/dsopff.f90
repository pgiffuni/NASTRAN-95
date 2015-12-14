SUBROUTINE dsopff ( dsname, iunit, istatus )
     
 CHARACTER (LEN=80), INTENT(OUT)          :: dsname
 INTEGER, INTENT(OUT)                     :: iunit
 INTEGER, INTENT(OUT)                     :: istatus
 
 COMMON / system / sysbuf, iwr
 COMMON / machin / mac(3), lqro
 INCLUDE          'DSIOF.COM'
 
 nbuff4 = nbuff * ( MOD(lqro,100) / 10 )
 OPEN  ( UNIT=iunit, FILE=dsname, RECL=nbuff4, FORM='UNFORMATTED'  &
     ,       ACCESS='DIRECT', IOSTAT=istatus, ERR=701 ,       STATUS='UNKNOWN' )
 GO TO 777
 701   WRITE ( iwr, 901 ) iunit, istatus, dsname
 901   FORMAT(//,' FATAL ERROR IN DSOPFF, UNABLE TO OPEN UNIT=',i4  &
     ,' IOSTAT=',i5 ,/,' FILE NAME=',a80 )
 iccerr = istatus
 CALL dsmsg  ( 101 )
 CALL mesage ( -61, 0, 0 )
 777   RETURN
END SUBROUTINE dsopff
