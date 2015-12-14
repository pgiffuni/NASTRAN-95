SUBROUTINE ifp5a (num)
     
!     IFP5A PRINTS MESSAGE NUMBER LINE ONLY.
!     CALLING SUBROUTINE PRINTS THE MESSAGE.
 
 
 INTEGER, INTENT(IN)                      :: num
 LOGICAL :: nogo
 INTEGER :: output
 CHARACTER (LEN=23) :: ufm
 COMMON /xmssg / ufm
 COMMON /system/ sysbuf,output,nogo
 
 CALL page2 (4)
 i = num + 4080
 WRITE  (output,10) ufm,i
 10 FORMAT (a23,i15,1H.)
 nogo = .true.
 RETURN
END SUBROUTINE ifp5a
