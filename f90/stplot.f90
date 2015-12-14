SUBROUTINE stplot (pltnum)
     
 
 INTEGER, INTENT(IN OUT)                  :: pltnum
 INTEGER :: date(3),idte(8),CHAR,ploter,pltype,pltape, eof,camera,bframs
 REAL :: SAVE(2,2)
 COMMON /xxparm/ pbufsz,camera,bframs
 COMMON /system/ ksystm(65)
 COMMON /char94/ CHAR(60)
 COMMON /pltdat/ model,ploter,reg(2,2),xymax(13),chrscl,skpa1(3),  &
     cntx,skpa2(5),pltype,pltape,skpa3,eof
 EQUIVALENCE     (ksystm(15),date(1))
 DATA    idte  / 2*1H ,1H/, 2*1H , 1H/, 2*1H  /, lstplt, m / 0,0 /
 
 IF (pltnum < 0) GO TO 150
 
!     SELECT THE PROPER CAMERA
 
 CALL selcam (camera,pltnum,0)
 
!     GENERATE THE ID PLOT
 
 IF (ploter /= lstplt) CALL skpfrm (1)
 lstplt = ploter
 CALL idplot (id)
 IF (id == 0) GO TO 120
 CALL selcam (camera,pltnum,0)
 CALL skpfrm (1)
 
!     INSERT THE BLANK FRAMES ON FILM ONLY
 
 120 IF (camera == 2 .OR. IABS(pltype) /= 1) GO TO 130
 IF (bframs == 0) GO TO 130
 CALL selcam (1,0,1)
 CALL skpfrm (MAX0(bframs,1))
 130 CALL selcam (camera,0,1)
 
!     TYPE THE PLOT NUMBER IN UPPER LEFT AND RIGHT CORNERS OF THE PLOT
 
 IF (pltnum == 0) GO TO 135
 DO  i  = 1,2
   SAVE(i,1) = reg(i,1)
   reg (i,1) = 0.
   SAVE(i,2) = reg(i,2)
   reg (i,2) = xymax(i)
 END DO
 CALL typint (0,0,0,0,0,-1)
 CALL typint (reg(1,1)+chrscl,reg(2,2)-chrscl,+1,pltnum,0,0)
 
!     PRINT THE DATE
 
 IF (m /= 0) GO TO 1312
 DO  n = 1,7,3
   m = m + 1
   i = date(m)/10 + 1
   j = date(m) - (i-1)*10 + 1
   IF (i == 1) i = 48
   idte(n  ) = CHAR(i)
   idte(n+1) = CHAR(j)
 END DO
 
 1312 CALL tipe (8.*cntx,reg(2,2)-chrscl,1,idte(1),8,0)
 
 CALL typint (reg(1,2)-chrscl,reg(2,2)-chrscl,-1,pltnum,0,0)
 DO  i = 1,2
   reg(i,1) = SAVE(i,1)
   reg(i,2) = SAVE(i,2)
 END DO
 135 CALL typint (0,0,0,0,0,1)
 GO TO 200
 
!     TERMINATE A PLOT
 
 150 CALL skpfrm (1)
 CALL typint (0,0,0,0,0,1)
 IF (eof == 0) CALL seof (pltape)
 CALL sclose (pltape)
 
 200 RETURN
END SUBROUTINE stplot
