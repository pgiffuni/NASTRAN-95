SUBROUTINE clstab (FILE,opt)
     
    INTEGER, INTENT(IN)                      :: FILE
    INTEGER, INTENT(IN OUT)                  :: opt
    INTEGER :: trlr(7)
    DATA trlr / 6*0,1 /
 
    trlr(1) = FILE
    CALL CLOSE (trlr,opt)
    CALL wrttrl (trlr)

    RETURN
END SUBROUTINE clstab
