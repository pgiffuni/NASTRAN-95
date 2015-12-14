SUBROUTINE ssghtp(order,z,lz)
!*****
!  SPECIAL IN-PLACE PARTITIONING ROUTINE USED ONLY BY SSGHT MODULE.
!*****
 
 
 INTEGER, INTENT(IN OUT)                  :: order(lz)
 INTEGER, INTENT(IN OUT)                  :: z(lz)
 INTEGER, INTENT(IN)                      :: lz
 INTEGER :: save1, save2, ptr
 i = 1
 isave = 1
 
!     CHECK TO SEE THAT POINTER TO NEXT SLOT HAS NOT BEEN USED YET.
 
 10 ptr = order(i)
 IF( ptr > 1000000 ) GO TO 40
 order(i) = ptr + 1000000
 
!     IF THE MOVE-TO LOCATION IS THE SAME, THEN DO NOTHING.
 
 IF( ptr == i ) GO TO 40
 
!     SAVE VALUE CURRENTLY IN SLOT WE ARE MOVING TO.
 
 save1 = z(ptr)
 
!     MOVE ITEM INTO SLOT
 
 z(ptr) = z(i)
 
!     SET POINTER TO WHERE -SAVE1- IS TO BE MOVED.
 
 20 jptr = order(ptr)
 IF( jptr > 1000000 ) GO TO 30
 order(ptr) = jptr + 1000000
 save2 = z(jptr)
 z(jptr) = save1
 save1 = save2
 ptr = jptr
 GO TO 20
 
!     END OF CHAIN.  GO BACK AND LOOK FOR ANOTHER.
 
 30 i = isave + 1
 GO TO 50
 40 i = i + 1
 50 isave = i
 IF( i <= lz ) GO TO 10
 
!     CLEAR OUT FLAGS AND RETURN.
 
 DO  i = 1,lz
   order(i) = order(i) - 1000000
 END DO
 RETURN
END SUBROUTINE ssghtp
