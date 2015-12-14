FUNCTION korsz ( i )
     
 INTEGER, INTENT(IN OUT)                  :: i
 COMMON / logout / lout
 COMMON / lstadd / lastad
 
 korsz = lastad - locfx(i) + 1
 CALL sswtch ( 13, l13 )
 IF ( l13 /= 0 ) WRITE ( lout, 2000 ) korsz
 RETURN
 2000  FORMAT(22X,' --- OPEN CORE =',i8,' WORDS (DECIMAL) ---')
END FUNCTION korsz
