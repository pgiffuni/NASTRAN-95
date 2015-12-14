SUBROUTINE chkopn (name)
     
    !     CHECKS IF A CALL TO SOFOPN HAS BEEN MADE.
 
 
    INTEGER, INTENT(IN OUT)                  :: name(2)
    LOGICAL :: opnsof
 
    CHARACTER (LEN=23) :: ufm
    COMMON /xmssg / ufm
    COMMON /sofcom/ sofdum(25),opnsof
    COMMON /system/ nbuff,nout
 
    IF (opnsof) GO TO 20
    WRITE  (nout,10) ufm,name
10  FORMAT (a23,' 6204, SUBROUTINE ',2A4,' - THE SUBROUTINE SOFOPN ',  &
        'SHOULD BE CALLED PRIOR TO ANY OF THE SOF UTILITY ', 'SUBROUTINES.')
    CALL mesage (-61,0,0)

20  RETURN
END SUBROUTINE chkopn
