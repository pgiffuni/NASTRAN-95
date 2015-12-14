SUBROUTINE skpfrm (bframs)
     
 
 INTEGER, INTENT(IN OUT)                  :: bframs
 INTEGER :: bfrms,ploter,camera,a(10),adv10(3),con10
 REAL :: SAVE(2,4),xymax(2)
 COMMON /pltdat/ model,ploter,reg(2,2),axymax(2),edge(2),camera,  &
     skpplt(9),pxymax(7),origin(2)
 DATA    adv10 , con10 / 1,2,3, 3 /
 
 DO  i = 1,2
   SAVE(i,1) = reg(i,1)
   reg(i,1)  = 0.
   SAVE(i,2) = reg(i,2)
   reg(i,2)  = axymax(i)+2.*edge(i)
   SAVE(i,3) = origin(i)
   origin(i) = 0.
   SAVE(i,4) = edge(i)
   edge(i)   = 0.
 END DO
 xymax(1) = AMAX1(reg(1,2),reg(2,2))
 xymax(2) = AMIN1(reg(1,2),reg(2,2))
 reg(1,2) = xymax(1)
 reg(2,2) = xymax(2)
 bfrms    = MIN0(MAX0(bframs,1),5)
 
!     PLOTTER 1, 2
 
 a(1) = con10
 a(2) = adv10(camera)
 DO  i = 3,6
   a(i) = 0
 END DO
 DO  i = 1,bfrms
   CALL wplt10 (a,0)
 END DO
 GO TO 300
 
 300 DO  i = 1,2
   reg(i,1) = SAVE(i,1)
   reg(i,2) = SAVE(i,2)
   origin(i)= SAVE(i,3)
   edge(i)  = SAVE(i,4)
 END DO
 
 RETURN
END SUBROUTINE skpfrm
