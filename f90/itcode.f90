FUNCTION itcode (itemx)
     
!     THE FUNCTION RETURNS AN INTEGER CODE NUMBER FOR ITEM.  THE CODE
!     NUMBER IS USED IN UPDATING THE MDI.  IF AN INCORRECT ITEM NAME IS
!     USED, THE VALUE RETURNED WILL BE -1.
 
 
 INTEGER, INTENT(IN OUT)                  :: itemx
 COMMON /itemdt/ nitem,item(7,1)
 COMMON /sys   / sys(5),ifrst
 
 DO  i = 1,nitem
   IF (itemx == item(1,i)) GO TO 20
 END DO
 
!     INVALID ITEM - RETURN -1
 
 itcode = -1
 RETURN
 
!     ITEM FOUND - RETURN MDI POSITION POINTER
 
 20 itcode = i + ifrst - 1
 RETURN
END FUNCTION itcode
