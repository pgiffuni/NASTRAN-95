SUBROUTINE sdr3
     INTEGER :: ofpfil(6)
 
 COMMON /system/ sysbuf, l
!*****
!  MAIN DRIVER FOR THE SDR3 MODULE...
!*****
 CALL sdr3a( ofpfil(1) )
!*****
!  IF ANY OF THE SIX DATA-BLOCKS DID NOT COMPLETE SORT2 CALL OFPDMP
!*****
 DO  i = 1,6
   IF( ofpfil(i) == 0 ) CYCLE
   WRITE(l,15)i,ofpfil(i)
   15 FORMAT(1H1,20(131(1H*)/),95H0DUE TO errors mentioned previously, s  &
       dr3 is calling the -ofp- TO output sdr3-INPUT-DATA-BLOCK-,i2,17H i  &
       n sort-i FORMAT/ 28H the sdr3 traceback NUMBER =,i3//20(131(1H*)/ ))
   ifile = i + 100
   CALL ofpdmp( ifile )
 END DO
 RETURN
END SUBROUTINE sdr3
