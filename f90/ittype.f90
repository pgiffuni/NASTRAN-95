FUNCTION ittype(itemx)
     
!*****
 
!     THIS FUNCTION RETURNS AN INTEGER CODE NUMBER TO INDICATE
!     WHETHER A PARTICULAR SOF ITEM IS A MATRIX OR TABLE.
!     THE RETURN CODES ARE
 
!          1 - MATRIX ITEM
!          0 - TABLE ITEM
!         -1 - ILLEGAL ITEM NAME
 
!*****
 
 
 
 INTEGER, INTENT(IN OUT)                  :: itemx
 COMMON / itemdt /       nitem    ,item(7,1)
 
 DO  i=1,nitem
   IF(itemx == item(1,i)) GO TO 20
 END DO
 
!     ILLIGAL ITEM - RETURN -1
 
 ittype = -1
 RETURN
 
!     ITEM FOUND - RETURN TYPE
 
 20 IF(item(2,i) <= 0) ittype = 0
 IF(item(2,i) > 0) ittype = 1
 RETURN
END FUNCTION ittype
