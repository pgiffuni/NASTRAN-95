SUBROUTINE feerdd
!*******
 
!     SUBROUTINE TO INITIALIZE COMMON /FEERCX/
 
!*******
 INTEGER :: jfrcx(28)
 INTEGER :: kfrcx( 4)
 INTEGER :: lfrcx( 4)
 
 COMMON   /feercx/  ifrcx(37)
 
 DATA               jfrcx   / 101,6*0   ,102,6*0  ,201,6*0  ,202,6*0  /
 DATA               kfrcx   / 301       ,302      ,303      ,304      /
 DATA               lfrcx   / 305       ,306      ,307      ,308      /
 DATA               mfrcx   /   204   /
 
 DO  i = 1,28
   ifrcx(i) = jfrcx(i)
 END DO
 DO  i = 1,4
   ifrcx(i+28) = kfrcx(i)
   ifrcx(i+32) = lfrcx(i)
 END DO
 ifrcx(37) = mfrcx
 
 RETURN
END SUBROUTINE feerdd
