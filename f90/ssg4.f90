SUBROUTINE ssg4
     
!     DRIVER TO DO INERTIAL RELIEF PORTION OF SSG
 
!     DMAP SEQUENCE
 
!     SSG4  PL,QR,PO,MR,MLR,D,MLL,MOOB,MOAB,GO,USET/PLI,POI/V,N,IOMT $
 
 INTEGER :: GO,uset
 INTEGER :: pl,qr,po,d,pli,poi,scr1,scr2,scr3,scr4,scr5
 COMMON /bitpos/ um,uo,ur,usg,usb,ul,ua,uf,us,un,ug
 COMMON /BLANK/ iomt
 DATA pl,qr,po,mr,mlr,d,mll,moob,moab,pli,poi,scr1,scr2,scr3,scr4  &
     ,scr5,GO,uset  &
     / 101,102,103,104,105,106,107,108,109,201,202,301,302,303,304, 305,110,111/
 
!     COMPUTE  MR-1*QR=TEMP2
 
 CALL factor(mr,scr1,scr2,scr3,scr4,scr5)
 CALL ssg3a( mr, scr1, qr, scr3, scr4, scr5, -1, xxx )
 
!     COMPUTE  MLL*D+MLR=TEMP1
 
 CALL ssg2b(mll,d,mlr,scr4,0,2,1,scr1)
 
!     COMPUTE  TEMP1*TEMP2+PL=PLI
 
 CALL ssg2b(scr4,scr3,pl,pli,0,2,1,scr1)
 IF(iomt > 0) THEN
   GO TO    10
 ELSE
   GO TO    20
 END IF
 
!     COMPUTE  MOOB*GO+MOAB=SCR4
 
 10 CALL ssg2b(moob,GO,moab,scr4,0,2,1,scr1)
 
!     COMPUTE DI*TEMP2  =SCR2
 
 CALL ssg2b(d,scr3,0,scr2,0,2,1,scr1)
 CALL sdr1b(scr5,scr2,scr3,scr1,ua,ul,uo,uset,0,0)
 
!     COMPUTE  SCR4*SCR1+PO=POI
 
 CALL ssg2b(scr4,scr1,po,poi,0,2,1,scr3)
 20 RETURN
END SUBROUTINE ssg4
