SUBROUTINE ofp1
     
!     THIS ROUTINE OUTPUTS A PAGE HEADING BASED ON PARAMETERS COMING
!     THROUGH COMMON.
!     THIS ROUTINE CALLS OPF1A, OFP1B OR OFP1C FOR ACTUAL PRINTING, SUCH
!     THAT OFP1A, OFP1B AND OFP1C CAN BE OVERLAYED IN PARALLEL.
 
 INTEGER :: l123(5),id(50),of(6)
 COMMON /system/ ksys(100)
 COMMON /zzzzzz/ core(1)
 EQUIVALENCE     (core(1),of(1),l123(1)), (id(1),of(6)),  &
     (nout,ksys(2)), (linet,ksys(12)), (iflag,ksys(33))
 
!     IFLAG IS WORD 33 OF /SYSTEM/ AND IS SET TO INCIDATE OFP PRINTED
!     LAST.
 
 CALL page1
 iflag = 1
 DO  i = 1,5
   line = l123(i)
   IF (line < 0) THEN
     GO TO  1000
   ELSE IF (line == 0) THEN
     GO TO   800
   END IF
   500  IF (line > 174) GO TO 600
   
! ... 1 THRU 174-
   
   CALL ofp1a (line)
   CYCLE
   600  IF (line > 380) GO TO 700
   
! ... 175 THRU 380 -
   
   CALL ofp1b (line)
   CYCLE
   
! ... 381 UP -
   
   700  CALL ofp1c (line)
   CYCLE
   
   800  WRITE  (nout,900)
   900  FORMAT (1H )
   
 END DO
 linet = linet + 4
 RETURN
END SUBROUTINE ofp1
