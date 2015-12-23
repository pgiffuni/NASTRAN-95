SUBROUTINE ifp4e (id)
     
    !     IFP4E, CALLED BY IFP4, CHECKS TO SEE THAT ID IS WITHIN PERMISSABLE
    !     RANGE OF FROM 1 TO 499999.
 
    INTEGER, INTENT(IN OUT)                  :: id
    LOGICAL :: nogo
    INTEGER :: output
    CHARACTER (LEN=23) :: ufm
    COMMON /xmssg / ufm
    COMMON /system/ sysbuf,output
 
    IF (id < 1) GO TO 100
    IF (id <= 499999) RETURN
 
    !     ERROR
 
100 nogo = .true.
    WRITE  (output,110) ufm,id
110 FORMAT (a23,' 4041, ID =',i12,' IS OUT OF PERMISSIBLE RANGE OF 1',  &
        ' TO 499999.')

    RETURN
END SUBROUTINE ifp4e
