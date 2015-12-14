SUBROUTINE apdcs
     INTEGER :: cp,acsid,cstma,iz(1)
 REAL :: rcp1(3),rcp4(3),rb1(3),rb4(3),rx1(3),rx4(3)  &
     ,              ra2(3),ra3(3),ra4(3),rb2(3),rb3(3),vx1(3),vx2(3)  &
     ,   vx3(3),acpl(3,3),v1(3),v2(3)
 COMMON /apd1c/ eid,pid,cp,nspan,nchord,lspan,lchord,igid  &
     ,              x1,y1,z1,x12,x4,y4,z4,x43,xop,x1p,alzo,mcstm  &
     ,              ncst1,ncst2,cidbx,acsid,iacs,silb,ncrd  &
     ,              scr1,scr2,scr3,scr4,scr5,ecta,bgpa,gpla,useta,sila  &
     ,              cstma,acpt,buf10,buf11,buf12,next,left,isiln
 COMMON /apd1d/ icpl(14),yp4,sg,cg,xp2,xp3,xp4,ra1(3)
 COMMON /zzzzzz/  z(1)
 EQUIVALENCE (z(1),iz(1))
 EQUIVALENCE (icpl(3),rb1(1)),(icpl(6),acpl(1,1))  &
     , (v1(1),rcp1(1)),(v2(1),rcp4(1))
 DATA degr/.017453293/
 icpl(2) = 1
! CREATE PANEL COORDINATE SYSTEM
! FIND CP TRANSFORMATION AND CONVERT POINT 1 AND 4 TO BASIC
 IF(cp == 0) GO TO 120
 IF(ncst1 == 0) GO TO 470
 DO  icp=ncst1,ncst2,14
   IF(iz(icp) == cp) GO TO 50
 END DO
 GO TO 470
 50 IF(iz(icp+1)-2 < 0) THEN
   GO TO    60
 ELSE IF (iz(icp+1)-2 == 0) THEN
   GO TO    70
 ELSE
   GO TO    80
 END IF
! CP RECTANGULAR
 60 rcp1(1)=x1
 rcp1(2)=y1
 rcp1(3)=z1
 rcp4(1)=x4
 rcp4(2)=y4
 rcp4(3)=z4
 GO TO 90
! CP CYLINDRICAL
 70 rcp1(1)=x1*COS(y1*degr)
 rcp1(2)=x1*SIN(y1*degr)
 rcp1(3)=z1
 rcp4(1)=x4*COS(y4*degr)
 rcp4(2) = x4*SIN(y4*degr)
 rcp4(3)=z4
 GO TO 90
! CP SPHERICAL
 80 rcp1(1)=x1*SIN(y1*degr)*COS(z1*degr)
 rcp1(2)=x1*SIN(y1*degr)*SIN(z1*degr)
 rcp1(3)=x1*COS(y1*degr)
 rcp4(1)=x4*SIN(y4*degr)*COS(z4*degr)
 rcp4(2)=x4*SIN(y4*degr)*SIN(z4*degr)
 rcp4(3)=x4*COS(y4*degr)
! CONVERT TO BASIC
 90 CALL gmmats(z(icp+5),3,3,0,rcp1,3,1,0,rb1)
 CALL gmmats(z(icp+5),3,3,0,rcp4,3,1,0,rb4)
 j=icp+1
 DO  i=1,3
   k=j+i
   rb1(i)=rb1(i)+z(k)
 END DO
 DO  i=1,3
   k=j+i
   rb4(i)=rb4(i)+z(k)
 END DO
 GO TO 130
! COORDS ARE IN BASIC
 120 rb1(1)=x1
 rb1(2)=y1
 rb1(3)=z1
 rb4(1)=x4
 rb4(2)=y4
 rb4(3)=z4
! FIND R1 THRU IN R4 AERO CS
 130 IF(acsid == 0) GO TO 150
 j=iacs+1
 DO  i=1,3
   k=j+i
   rx1(i)=rb1(i)-z(k)
   rx4(i)=rb4(i)-z(k)
 END DO
 CALL gmmats(z(iacs+5),3,3,1,rx1,3,1,0,ra1)
 CALL gmmats (z(iacs+5),3,3,1,rx4,3,1,0,ra4)
 GO TO 170
 150 DO  i=1,3
   ra1(i)=rb1(i)
   ra4(i)=rb4(i)
 END DO
 
!     STOP IF BODY
 
 IF(igid < 0) GO TO 1000
! CALCULATE R2 AND R3 IN AC CS
 170 DO  i=2,3
   ra2(i)=ra1(i)
   ra3(i)=ra4(i)
 END DO
 ra2(1)=ra1(1)+x12
 ra3(1)=ra4(1)+x43
 ee=SQRT((ra4(3)-ra1(3))**2 + (ra4(2)-ra1(2))**2)
 sg=(ra4(3)-ra1(3))/ee
 cg=(ra4(2)-ra1(2))/ee
! LOCATE POINTS 2,3,4 IN PANEL CORDINATE SYSTEM
 xp2=x12
 xp4=ra4(1)-ra1(1)
 xp3=ra3(1)-ra1(1)
 yp4=ee
! TRANSFORM R2 AND R3 INTO BASIC
 IF(acsid == 0) GO TO 200
 CALL gmmats(z(iacs+5),3,3,0,ra2,3,1,0,rb2)
 CALL gmmats(z(iacs+5),3,3,0,ra3,3,1,0,rb3)
 j=iacs+1
 DO  i=1,3
   k=j+i
   rb2(i) = rb2(i) + z(k)
   rb3(i) = rb3(i) + z(k)
 END DO
 GO TO 220
 200 DO  i=1,3
   rb2(i)=ra2(i)
   rb3(i)=ra3(i)
 END DO
! FIND PANEL COORDINATE SYSTEM
 220 DO  i=1,3
   vx1(i)=rb2(i)-rb1(i)
   vx2(i)=rb4(i)-rb1(i)
   vx3(i) = rb3(i) - rb1(i)
   IF ( x12. EQ. 0.0 ) vx1(i) = vx3(i)
 END DO
 CALL saxb(vx1,vx2,v1)
 sx1=sadotb(v1,v1)
 CALL saxb(vx1,vx3,v2)
 sx2=sadotb(v2,v2)
 IF(sx1 < sx2) GO TO 250
 sx1=1.0/SQRT(sx1)
 DO  i=1,3
   vx3(i)=v1(i)*sx1
 END DO
 GO TO 270
 250 sx2=1.0/SQRT(sx2)
 DO  i=1,3
   vx3(i)=v2(i)*sx2
 END DO
 270 IF(acsid /= 0) GO TO 275
 vx1(1) = 1.0
 vx1(2) = 0.0
 vx1(3) = 0.0
 GO TO 285
 275 j=iacs+5
 DO  i=1,3
   k=j+3*(i-1)
   vx1(i)=z(k)
 END DO
 285 CONTINUE
 CALL saxb(vx3,vx1,vx2)
 DO  i=1,3
   acpl(1,i)=vx1(i)
   acpl(2,i)=vx2(i)
   acpl(3,i)=vx3(i)
 END DO
! WRITE TRANSFORMATION ON CSTMA
 icpl(1)=mcstm
 CALL WRITE(cstma,icpl(1),14,0)
 1000 RETURN
 470 CALL mesage(-30,25,cp)
 GO TO 1000
END SUBROUTINE apdcs
