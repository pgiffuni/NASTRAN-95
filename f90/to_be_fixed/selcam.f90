SUBROUTINE selcam (camera,pltnum,opt)
     
 
 INTEGER, INTENT(IN OUT)                  :: camera
 REAL, INTENT(IN OUT)                     :: pltnum
 INTEGER, INTENT(IN OUT)                  :: opt
 REAL :: edge(2)  ,origin(2),xymax(2) ,SAVE(2,3)
 INTEGER :: a(17)    , camnum   , ploter   ,pltnum ,       con10(2) ,cam10(3)
 
 COMMON /pltdat/  pd(20,2)
 
 EQUIVALENCE (   model,pd( 1,1)) , (  ploter,pd( 2,1))  &
     ,           (xymax(1),pd( 7,1)) , ( edge(1),pd( 9,1))  &
     ,           (  camnum,pd(11,1)) ,(origin(1),pd( 8,2))
 DATA con10,cam10 / 1,2, 1,2,3 /
 
 DO  i = 1,2
   SAVE(i,1) = edge(i)
   edge(i) = 0.
   SAVE(i,2) = origin(i)
   origin(i) = 0.
   SAVE(i,3) = xymax(i)
   xymax(i) = 0.
   a(i) = IABS(pltnum)
 END DO
 camnum = MIN0(MAX0(camera,1),3)
 IF(opt == 0.0) THEN
   GO TO   160
 ELSE
   GO TO   165
 END IF
 
!     PLOTTER 1, 2
 
 160 a(3) = a(1)
 a(1) = con10(1)
 a(2) = 0
 a(4) = SAVE(1,3) + 2.*SAVE(1,1) + .1
 a(5) = SAVE(2,3) + 2.*SAVE(2,1) + .1
 a(6) = 0
 CALL wplt10 (a,0)
 165 a(1) = con10(2)
 a(2) = cam10(camnum)
 DO  i = 3,6
   a(i) = 0
 END DO
 CALL wplt10 (a,0)
 
 DO  i=1,2
   edge(i) = SAVE(i,1)
   origin(i) = SAVE(i,2)
   xymax(i) = SAVE(i,3)
 END DO
 RETURN
END SUBROUTINE selcam
